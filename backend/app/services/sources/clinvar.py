"""
ClinVar client (NCBI E-utilities) -- cache-first, rate-limited, parallel, time-bounded.

IMPORTANT: this search is by genomic POSITION, so it can return a different allele at the same
position. The scoring engines therefore use the exact-allele ClinVar data that arrives through
MyVariant (see evidence_service.clinvar_from_myvariant). The NCBI result is kept only as
`clinvar_position` (information / a review flag) and, to save time, is cache-only unless
CLINVAR_NCBI_ENABLED=1.

Without an NCBI_API_KEY NCBI allows 3 requests/s, with a key 10/s. Each uncached variant needs up to
two requests (esearch + esummary), so ClinVar is the rate-limit-bound source: put a key in .env.
"""

from __future__ import annotations

import os
from typing import Callable

import httpx

from app.models.variant import Variant
from app.services.sources.base import (
    CircuitBreaker, EvidenceResult, RateLimiter, SourceStats, make_timeout,
    not_assessable, request_with_retry, run_lookups,
)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TABLE = "clinvar_evidence"
WORKERS = 3


def _variant_key(variant: Variant) -> str | None:
    if not (variant.chromosome and variant.position and variant.ref and variant.alt):
        return None
    chrom = variant.chromosome[3:] if variant.chromosome.lower().startswith("chr") else variant.chromosome
    return f"{chrom}:{variant.position}:{variant.ref}:{variant.alt}"


class ClinVarClient:
    def __init__(self, db_path: str = "data/variant_classifier.db", timeout: float = 15.0):
        self.db_path = db_path
        self.timeout = timeout
        self.api_key = os.environ.get("NCBI_API_KEY", "").strip() or None
        self._limiter = RateLimiter(0.11 if self.api_key else 0.35)
        self.stats = SourceStats()
        self.max_seconds: float | None = None

    def _params(self, extra: dict) -> dict:
        p = {"db": "clinvar", "retmode": "json", **extra}
        if self.api_key:
            p["api_key"] = self.api_key
        return p

    def _get(self, client: httpx.Client, endpoint: str, params: dict):
        self._limiter.wait()
        return request_with_retry(client, "GET", f"{EUTILS_BASE}/{endpoint}", params=self._params(params))

    def _live_fetch(self, variant: Variant, client: httpx.Client) -> tuple[str, dict]:
        """Returns (status, data) from the live API -- never writes to the cache itself."""
        chrom = variant.chromosome[3:] if variant.chromosome.lower().startswith("chr") else variant.chromosome
        resp = self._get(client, "esearch.fcgi", {"term": f"{chrom}[chr] AND {variant.position}[chrpos37]", "retmax": "5"})
        if resp is None or resp.status_code != 200:
            return "SOURCE_UNAVAILABLE", {"error": f"esearch failed: HTTP {resp.status_code if resp is not None else 'none'}"}
        try:
            id_list = resp.json().get("esearchresult", {}).get("idlist", [])
        except ValueError:
            return "SOURCE_UNAVAILABLE", {"error": "invalid esearch JSON"}
        if not id_list:
            return "AVAILABLE", {"found": False}

        resp2 = self._get(client, "esummary.fcgi", {"id": ",".join(id_list)})
        if resp2 is None or resp2.status_code != 200:
            return "SOURCE_UNAVAILABLE", {"error": f"esummary failed: HTTP {resp2.status_code if resp2 is not None else 'none'}"}
        try:
            result = resp2.json().get("result", {})
            records = [result[uid] for uid in id_list if uid in result]
        except ValueError:
            return "SOURCE_UNAVAILABLE", {"error": "invalid esummary JSON"}

        data = {
            "found": True,
            "variation_ids": id_list,
            "records": [
                {
                    "clinical_significance": (r.get("germline_classification") or {}).get("description"),
                    "review_status": (r.get("germline_classification") or {}).get("review_status"),
                    "condition": [t.get("trait_name") for t in (r.get("germline_classification") or {}).get("trait_set", [])],
                    "accession": r.get("accession"),
                    "last_evaluated": (r.get("germline_classification") or {}).get("last_evaluated"),
                }
                for r in records
            ],
        }
        return "AVAILABLE", data

    def annotate_batch(self, variants: list[Variant], progress: Callable[[int], None] | None = None) -> dict[str, EvidenceResult]:
        self.stats = SourceStats()
        breaker = CircuitBreaker(threshold=5, deadline_seconds=self.max_seconds)
        results: dict[str, EvidenceResult] = {}
        groups: dict[str, list[Variant]] = {}
        for v in variants:
            key = _variant_key(v)
            if key is None:
                results[v.variant_id] = not_assessable("clinvar", "missing coordinates")
                if progress:
                    progress(1)
            else:
                groups.setdefault(key, []).append(v)
        enabled = os.environ.get("CLINVAR_NCBI_ENABLED", "0") == "1"
        reason = None if enabled else "NCBI position search off (exact-allele ClinVar comes from MyVariant)"
        with httpx.Client(timeout=make_timeout(self.timeout)) as client:
            results.update(run_lookups(source="clinvar", table=TABLE, groups=groups,
                                       fetch_for=lambda v: self._live_fetch(v, client), db_path=self.db_path,
                                       breaker=breaker, stats=self.stats, workers=WORKERS, progress=progress,
                                       cache_only_reason=reason))
        return results
