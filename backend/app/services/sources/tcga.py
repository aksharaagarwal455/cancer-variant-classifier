"""
TCGA evidence via the NCI GDC API (open access, somatic mutations).
DB-first like the other sources. GDC is GRCh38, so matching is by gene +
protein change (build-independent), not by genomic coordinate.
Counts are approximate: gene and aa_change filters are matched across
transcript consequences.
"""

from __future__ import annotations

import json
import re
from typing import Callable

import httpx

from app.models.variant import Variant
from app.services.sources.base import (
    CircuitBreaker, EvidenceResult, SourceStats, make_timeout, not_assessable, request_with_retry,
    run_lookups, skipped, source_enabled,
)

GDC_URL = "https://api.gdc.cancer.gov/ssm_occurrences"
TABLE = "tcga_evidence"
WORKERS = 4
_AA3 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}
_PC = re.compile(r"p\.\(?([A-Za-z]{1,3})(\d+)([A-Za-z\*]{1,3})\)?")


def _one(aa: str) -> str | None:
    if len(aa) == 1:
        return aa.upper() if aa != "*" else "*"
    return _AA3.get(aa.capitalize())


def parse_protein_change(hgvsp: str | None) -> str | None:
    """'p.V600E' / 'p.Val600Glu' -> 'V600E'. Frameshifts and unparseable -> None."""
    if not hgvsp:
        return None
    m = _PC.search(hgvsp)
    if not m:
        return None
    ref, pos, alt = m.groups()
    if "fs" in alt.lower():
        return None
    r, a = _one(ref), _one(alt)
    if not r or not a:
        return None
    return f"{r}{pos}{a}"


def fetch_tcga_cases(gene: str, aa: str, client: httpx.Client) -> tuple[str, dict]:
    filters = {"op": "and", "content": [
        {"op": "=", "content": {"field": "ssm.consequence.transcript.gene.symbol", "value": gene}},
        {"op": "=", "content": {"field": "ssm.consequence.transcript.aa_change", "value": aa}},
    ]}
    resp = request_with_retry(client, "GET", GDC_URL, max_retries=1, base_backoff_seconds=0.5,
                              params={"filters": json.dumps(filters), "size": 1, "fields": "ssm.ssm_id"})
    if resp is None or resp.status_code != 200:
        return "SOURCE_UNAVAILABLE", {"error": f"HTTP {resp.status_code if resp is not None else 'none'}"}
    try:
        total = resp.json()["data"]["pagination"]["total"]
    except (ValueError, KeyError, TypeError):
        return "SOURCE_UNAVAILABLE", {"error": "unexpected GDC response shape"}
    return "AVAILABLE", {"gene": gene, "aa_change": aa, "cases": int(total)}


class TcgaClient:
    def __init__(self, db_path: str = "data/variant_classifier.db", timeout: float = 10.0):
        self.db_path = db_path
        self.timeout = timeout
        self.stats = SourceStats()
        self.max_seconds: float | None = None

    def annotate_batch(self, variants: list[Variant], progress: Callable[[int], None] | None = None) -> dict[str, EvidenceResult]:
        self.stats = SourceStats()
        results: dict[str, EvidenceResult] = {}
        if not source_enabled("tcga"):
            for v in variants:
                results[v.variant_id] = skipped("tcga", "disabled (TCGA_ENABLED=0)")
            self.stats.skip(len(variants))
            if progress:
                progress(len(variants))
            return results

        breaker = CircuitBreaker(threshold=4, deadline_seconds=self.max_seconds)
        groups: dict[str, list[Variant]] = {}
        for v in variants:
            aa = parse_protein_change(v.hgvsp)
            if not v.gene or not aa:
                results[v.variant_id] = not_assessable("tcga", "no gene / substitution-type protein change to look up")
                if progress:
                    progress(1)
            else:
                groups.setdefault(f"{v.gene}:{aa}", []).append(v)

        with httpx.Client(timeout=make_timeout(self.timeout)) as client:
            results.update(run_lookups(
                source="tcga", table=TABLE, groups=groups,
                fetch_for=lambda v: fetch_tcga_cases(v.gene, parse_protein_change(v.hgvsp), client),
                db_path=self.db_path, breaker=breaker, stats=self.stats, workers=WORKERS, progress=progress))
        return results
