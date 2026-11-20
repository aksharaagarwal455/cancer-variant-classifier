"""VCF-style REF/ALT -> HGVS genomic notation, indels included (no reference FASTA needed)."""

from __future__ import annotations
from app.models.variant import Variant

def _trim_shared_prefix_suffix(ref: str, alt: str) -> tuple[str, str, int, int]:
    prefix_len = 0
    max_prefix = min(len(ref), len(alt))
    while prefix_len < max_prefix and ref[prefix_len] == alt[prefix_len]:
        prefix_len += 1
    suffix_len = 0
    max_suffix = min(len(ref), len(alt)) - prefix_len
    while suffix_len < max_suffix and ref[len(ref)-1-suffix_len] == alt[len(alt)-1-suffix_len]:
        suffix_len += 1
    return ref[prefix_len:len(ref)-suffix_len], alt[prefix_len:len(alt)-suffix_len], prefix_len, suffix_len

def build_indel_hgvs(variant: Variant) -> tuple[str | None, str | None]:
    ref, alt, pos = variant.ref, variant.alt, variant.position
    if not (ref and alt and pos):
        return None, "missing ref/alt/position"
    if ref == alt:
        return None, "ref equals alt, not a variant"
    if len(ref) == 1 and len(alt) == 1:
        return None, "single-base substitution, not an indel"
    chrom = variant.chromosome if variant.chromosome and variant.chromosome.lower().startswith("chr") else f"chr{variant.chromosome}"
    t_ref, t_alt, prefix_len, suffix_len = _trim_shared_prefix_suffix(ref, alt)
    if t_ref == "" and t_alt != "":
        left = pos + prefix_len - 1
        return f"{chrom}:g.{left}_{left+1}ins{t_alt}", None
    if t_alt == "" and t_ref != "":
        start = pos + prefix_len
        end = start + len(t_ref) - 1
        return (f"{chrom}:g.{start}del" if start == end else f"{chrom}:g.{start}_{end}del"), None
    if t_ref != "" and t_alt != "":
        start = pos + prefix_len
        end = start + len(t_ref) - 1
        return (f"{chrom}:g.{start}delins{t_alt}" if start == end else f"{chrom}:g.{start}_{end}delins{t_alt}"), None
    return None, "could not determine indel shape after trimming"
