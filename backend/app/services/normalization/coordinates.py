"""Genome build detection via landmark-gene coordinate evidence (never guessed from filename)."""

from __future__ import annotations
from collections import Counter
from typing import Iterable
from app.models.variant import GenomeBuild, Variant

LANDMARK_GENES: dict[str, tuple[str, int, int, int, int]] = {
    "TP53":  ("17", 7_565_000, 7_595_000, 7_655_000, 7_692_000),
    "EGFR":  ("7",  55_080_000, 55_330_000, 55_015_000, 55_275_000),
    "BRCA1": ("17", 41_190_000, 41_280_000, 43_040_000, 43_130_000),
    "BRCA2": ("13", 32_889_000, 32_974_000, 32_315_000, 32_400_000),
    "KRAS":  ("12", 25_355_000, 25_405_000, 25_205_000, 25_255_000),
}

def _strip_chr(chrom: str | None) -> str | None:
    if chrom is None:
        return None
    return chrom[3:] if chrom.lower().startswith("chr") else chrom

def detect_genome_build(variants: Iterable[Variant]) -> tuple[GenomeBuild, dict]:
    votes: Counter[str] = Counter()
    landmark_hits: dict[str, list] = {}
    for v in variants:
        if v.gene is None or v.gene not in LANDMARK_GENES or v.position is None:
            continue
        chrom_expected, hg19_start, hg19_end, hg38_start, hg38_end = LANDMARK_GENES[v.gene]
        if _strip_chr(v.chromosome) != chrom_expected:
            continue
        landmark_hits.setdefault(v.gene, []).append((v.variant_id, v.position))
        in_hg19 = hg19_start <= v.position <= hg19_end
        in_hg38 = hg38_start <= v.position <= hg38_end
        if in_hg19 and not in_hg38:
            votes["GRCh37"] += 1
        elif in_hg38 and not in_hg19:
            votes["GRCh38"] += 1

    evidence = {"landmark_variants_found": {g: len(h) for g, h in landmark_hits.items()}, "votes": dict(votes)}
    if not votes:
        return GenomeBuild.UNKNOWN, evidence
    top_build, top_count = votes.most_common(1)[0]
    if top_count == sum(votes.values()):
        return (GenomeBuild.GRCH37 if top_build == "GRCh37" else GenomeBuild.GRCH38), evidence
    evidence["conflict"] = True
    return GenomeBuild.UNKNOWN, evidence

def normalize_ref_alt(ref: str | None, alt: str | None) -> tuple[str | None, str | None]:
    r = ref.strip().upper() if ref else ref
    a = alt.strip().upper() if alt else alt
    return r, a

def apply_genome_build(variants: list[Variant], forced_build: str | None = None) -> tuple[list[Variant], dict]:
    import os
    build, evidence = detect_genome_build(variants)
    forced = (forced_build or os.environ.get("GENOME_BUILD") or "").strip().lower()
    if forced:
        mapping = {"hg19": GenomeBuild.GRCH37, "grch37": GenomeBuild.GRCH37,
                   "hg38": GenomeBuild.GRCH38, "grch38": GenomeBuild.GRCH38}
        if forced not in mapping:
            raise ValueError(f"unknown genome build {forced!r}; use hg19 or hg38")
        evidence["detected_build"] = build.name
        evidence["forced_by_user"] = forced
        if build != GenomeBuild.UNKNOWN and build != mapping[forced]:
            evidence["warning"] = "forced build disagrees with landmark-gene detection"
        build = mapping[forced]
    for v in variants:
        v.genome_build = build
        v.ref, v.alt = normalize_ref_alt(v.ref, v.alt)
    return variants, evidence
