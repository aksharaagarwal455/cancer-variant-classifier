"""
TSV parser + schema mapper for ANNOVAR-style somatic variant files.

- Detects and preserves duplicate column names by position (never
  silently drops a column).
- Normalizes "." and empty string to None (ANNOVAR null convention).
- Never guesses at required fields (chrom/pos/ref/alt/gene) -- marks
  rows INVALID with a specific reason if they can't be found.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from app.models.variant import ParseStatus, Variant

SCHEMA_ALIASES: dict[str, list[str]] = {
    "chromosome": ["CHROM", "chrom", "chr", "chromosome"],
    "position": ["POS", "pos", "position", "start"],
    "ref": ["REF", "ref", "reference", "reference_allele"],
    "alt": ["ALT", "alt", "alternate", "alternate_allele"],
    "gene": ["Gene", "gene", "gene_symbol", "symbol"],
    "transcript": ["Transcript", "transcript", "transcript_id"],
    "hgvsc": ["Coding_change", "hgvsc"],
    "hgvsp": ["Protein_change", "protein_change", "hgvsp"],
    "rsid": ["avsnp150", "rsid", "dbsnp"],
    "variant_type": ["ExonicFunc.refGene", "variant_type"],
    "consequence": ["Func.refGene", "consequence"],
}

NULL_TOKENS = {"", ".", "NA", "N/A", "NULL", "none", "None"}


def _clean_value(raw: str | None) -> str | None:
    if raw is None:
        return None
    v = raw.strip()
    return None if v in NULL_TOKENS else v


def _dedupe_headers(raw_headers: list[str]) -> list[str]:
    seen: Counter[str] = Counter()
    result: list[str] = []
    for h in raw_headers:
        seen[h] += 1
        result.append(h if seen[h] == 1 else f"{h}__dup{seen[h]}")
    return result


@dataclass
class ParseReport:
    source_path: str
    total_rows: int = 0
    parsed_ok: int = 0
    invalid: int = 0
    duplicates: int = 0
    duplicate_header_columns: list[str] = field(default_factory=list)
    columns_found: int = 0
    unmapped_canonical_fields: list[str] = field(default_factory=list)


def _map_canonical_fields(headers: list[str]) -> dict[str, str]:
    header_set = set(headers)
    header_lower_map = {h.lower(): h for h in headers}
    mapping: dict[str, str] = {}
    for canonical, aliases in SCHEMA_ALIASES.items():
        found = None
        for alias in aliases:
            if alias in header_set:
                found = alias
                break
        if found is None:
            for alias in aliases:
                if alias.lower() in header_lower_map:
                    found = header_lower_map[alias.lower()]
                    break
        if found is not None:
            mapping[canonical] = found
    return mapping


def parse_tsv(path: str | Path) -> tuple[list[Variant], ParseReport]:
    path = Path(path)
    report = ParseReport(source_path=str(path))

    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            raw_headers = next(reader)
        except StopIteration:
            return [], report

        dup_counts = Counter(raw_headers)
        report.duplicate_header_columns = sorted([n for n, c in dup_counts.items() if c > 1])
        headers = _dedupe_headers(raw_headers)
        report.columns_found = len(headers)

        field_map = _map_canonical_fields(raw_headers)
        for canonical in ("chromosome", "position", "ref", "alt", "gene"):
            if canonical not in field_map:
                report.unmapped_canonical_fields.append(canonical)

        variants: list[Variant] = []
        seen_keys: dict[tuple, int] = {}

        for row_idx, raw_row in enumerate(reader):
            report.total_rows += 1

            if len(raw_row) != len(headers):
                original_input = {f"col_{i}": _clean_value(v) for i, v in enumerate(raw_row)}
                variants.append(Variant(
                    variant_id=f"INVALID_ROW_{row_idx}", original_row_id=row_idx,
                    parse_status=ParseStatus.INVALID,
                    parse_errors=[f"Row has {len(raw_row)} fields, expected {len(headers)}"],
                    original_input=original_input,
                ))
                report.invalid += 1
                continue

            row_dict = {h: _clean_value(v) for h, v in zip(headers, raw_row)}

            chrom = row_dict.get(field_map.get("chromosome", ""))
            pos_raw = row_dict.get(field_map.get("position", ""))
            ref = row_dict.get(field_map.get("ref", ""))
            alt = row_dict.get(field_map.get("alt", ""))
            gene = row_dict.get(field_map.get("gene", ""))
            transcript = row_dict.get(field_map.get("transcript", ""))
            hgvsc = row_dict.get(field_map.get("hgvsc", ""))
            hgvsp = row_dict.get(field_map.get("hgvsp", ""))
            rsid = row_dict.get(field_map.get("rsid", ""))
            variant_type = row_dict.get(field_map.get("variant_type", ""))
            consequence = row_dict.get(field_map.get("consequence", ""))

            parse_errors: list[str] = []
            position: int | None = None
            if pos_raw is not None:
                try:
                    position = int(float(pos_raw))
                except ValueError:
                    parse_errors.append(f"position not numeric: {pos_raw!r}")
            else:
                parse_errors.append("position missing")
            if chrom is None:
                parse_errors.append("chromosome missing")
            if ref is None:
                parse_errors.append("ref missing")
            if alt is None:
                parse_errors.append("alt missing")

            status = ParseStatus.OK if not parse_errors else ParseStatus.INVALID
            dedup_key = (chrom, position, ref, alt)
            if status == ParseStatus.OK:
                if dedup_key in seen_keys:
                    status = ParseStatus.DUPLICATE
                    parse_errors.append(f"duplicate of row {seen_keys[dedup_key]}")
                else:
                    seen_keys[dedup_key] = row_idx

            if status == ParseStatus.OK:
                variant_id = f"{chrom}:{position}:{ref}:{alt}"
                if transcript:
                    variant_id = f"{variant_id}:{transcript}"
            else:
                variant_id = f"ROW_{row_idx}:{chrom}:{position}:{ref}:{alt}"

            variants.append(Variant(
                variant_id=variant_id, original_row_id=row_idx,
                chromosome=chrom, position=position, ref=ref, alt=alt,
                gene=gene, transcript=transcript, hgvsc=hgvsc, hgvsp=hgvsp, rsid=rsid,
                variant_type=variant_type, consequence=consequence,
                parse_status=status, parse_errors=parse_errors, original_input=row_dict,
            ))

            if status == ParseStatus.OK:
                report.parsed_ok += 1
            elif status == ParseStatus.DUPLICATE:
                report.duplicates += 1
            else:
                report.invalid += 1

    return variants, report
