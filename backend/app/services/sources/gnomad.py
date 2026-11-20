"""gnomAD client (GraphQL) -- cache-first. gnomAD rate-limits hard, so calls are spaced >= 1 s apart."""

from __future__ import annotations

import os
from typing import Callable

import httpx

from app.models.variant import GenomeBuild, Variant
from app.services.sources.base import (
    CircuitBreaker, EvidenceResult, RateLimiter, SourceStats, make_timeout, not_assessable,
    request_with_retry, run_lookups, skipped, source_enabled,
)

GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
TABLE = "gnomad_evidence"
WORKERS = 2

_BUILD_TO_DATASET = {GenomeBuild.GRCH37: "gnomad_r2_1", GenomeBuild.GRCH38: "gnomad_r4"}

_QUERY = """
query VariantAF($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    genome { ac an af }
    exome { ac an af }
  }
}
"""


def build_gnomad_variant_id(variant: Variant) -> str | None:
    if not (variant.chromosome and variant.position and variant.ref and variant.alt):
        return None
    chrom = variant.chromosome[3:] if variant.chromosome.lower().startswith("chr") else variant.chromosome
    return f"{chrom}-{variant.position}-{variant.ref}-{variant.alt}"


class GnomadClient:
    def __init__(self, db_path: str = "data/variant_classifier.db", timeout: float = 12.0):
        self.db_path = db_path
        self.timeout = timeout
        self._limiter = RateLimiter(float(os.environ.get("GNOMAD_MIN_INTERVAL", "1.0")))
        self.stats = SourceStats()
        self.max_seconds: float | None = None

    def _live_fetch(self, variant: Variant, gid: str, dataset: str, client: httpx.Client) -> tuple[str, dict]:
        self._limiter.wait()
        resp = request_with_retry(client, "POST", GNOMAD_API_URL, max_retries=1,
                                  json={"query": _QUERY, "variables": {"variantId": gid, "dataset": dataset}},
                                  headers={"Content-Type": "application/json"})
        if resp is None or resp.status_code != 200:
            return "SOURCE_UNAVAILABLE", {"error": f"HTTP {resp.status_code if resp is not None else 'none'}"}
        try:
            payload = resp.json()
        except ValueError:
            return "SOURCE_UNAVAILABLE", {"error": "invalid JSON"}
        if payload.get("errors"):
            return "AVAILABLE", {"found": False, "errors": payload["errors"]}
        variant_data = (payload.get("data") or {}).get("variant")
        return "AVAILABLE", {"found": variant_data is not None, "variant": variant_data}

    def annotate_batch(self, variants: list[Variant], progress: Callable[[int], None] | None = None) -> dict[str, EvidenceResult]:
        self.stats = SourceStats()
        results: dict[str, EvidenceResult] = {}
        if not source_enabled("gnomad"):
            for v in variants:
                results[v.variant_id] = skipped("gnomad", "disabled (GNOMAD_ENABLED=0)")
            self.stats.skip(len(variants))
            if progress:
                progress(len(variants))
            return results

        breaker = CircuitBreaker(threshold=3, deadline_seconds=self.max_seconds)
        groups: dict[str, list[Variant]] = {}
        meta: dict[str, tuple[str, str]] = {}
        for v in variants:
            if v.genome_build not in _BUILD_TO_DATASET:
                results[v.variant_id] = not_assessable("gnomad", "genome build unknown/unsupported")
            else:
                gid = build_gnomad_variant_id(v)
                if gid is None:
                    results[v.variant_id] = not_assessable("gnomad", "could not construct gnomAD id")
                else:
                    dataset = _BUILD_TO_DATASET[v.genome_build]
                    key = f"{dataset}:{gid}"
                    groups.setdefault(key, []).append(v)
                    meta[key] = (gid, dataset)
            if v.variant_id in results and progress:
                progress(1)

        with httpx.Client(timeout=make_timeout(self.timeout)) as client:
            def fetch_for(v: Variant):
                key = f"{_BUILD_TO_DATASET[v.genome_build]}:{build_gnomad_variant_id(v)}"
                gid, dataset = meta[key]
                return self._live_fetch(v, gid, dataset, client)

            results.update(run_lookups(source="gnomad", table=TABLE, groups=groups, fetch_for=fetch_for,
                                       db_path=self.db_path, breaker=breaker, stats=self.stats,
                                       workers=WORKERS, progress=progress))
        return results
