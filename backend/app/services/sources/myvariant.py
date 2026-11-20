"""MyVariant.info client -- cache-first, 50 variants per HTTP request (one POST for a whole chunk)."""

from __future__ import annotations

import time
from datetime import datetime, timezone
from typing import Callable

import httpx

from app.database.db import get_local, save_local
from app.models.variant import GenomeBuild, Variant
from app.services.normalization.hgvs import build_indel_hgvs
from app.services.sources.base import (
    CircuitBreaker, EvidenceResult, SourceStats, SourceStatus, is_offline, make_timeout, not_assessable,
    request_with_retry, run_parallel, skipped,
)

MYVARIANT_BASE_URL = "https://myvariant.info/v1"
TABLE = "myvariant_evidence"
CHUNK = 50
WORKERS = 4
_BUILD_TO_ASSEMBLY = {GenomeBuild.GRCH37: "hg19", GenomeBuild.GRCH38: "hg38"}


def build_hgvs_genomic_id(variant: Variant) -> str | None:
    if not (variant.chromosome and variant.position and variant.ref and variant.alt):
        return None
    if len(variant.ref) != 1 or len(variant.alt) != 1:
        return None
    if variant.ref not in "ACGT" or variant.alt not in "ACGT":
        return None
    chrom = variant.chromosome if variant.chromosome.lower().startswith("chr") else f"chr{variant.chromosome}"
    return f"{chrom}:g.{variant.position}{variant.ref}>{variant.alt}"


def resolve_hgvs_genomic_id(variant: Variant) -> tuple[str | None, str | None]:
    sub_id = build_hgvs_genomic_id(variant)
    if sub_id is not None:
        return sub_id, None
    if len(variant.ref or "") == 1 and len(variant.alt or "") == 1:
        return None, "substitution failed validation"
    return build_indel_hgvs(variant)


class MyVariantClient:
    def __init__(self, db_path: str = "data/variant_classifier.db", timeout: float = 20.0):
        self.db_path = db_path
        self.timeout = timeout
        self.stats = SourceStats()
        self.max_seconds: float | None = None

    def _live_fetch(self, hgvs_id: str, assembly: str, client: httpx.Client) -> tuple[str, dict]:
        """Single-id lookup (used by scripts/benchmark.py)."""
        status, by_id = self._live_fetch_many([hgvs_id], assembly, client)
        if status != "AVAILABLE":
            return status, by_id
        return "AVAILABLE", by_id.get(hgvs_id, {"notfound": True})

    def _live_fetch_many(self, ids: list[str], assembly: str, client: httpx.Client) -> tuple[str, dict]:
        resp = request_with_retry(client, "POST", f"{MYVARIANT_BASE_URL}/variant",
                                  data={"ids": ",".join(ids), "assembly": assembly})
        if resp is None or resp.status_code != 200:
            return "SOURCE_UNAVAILABLE", {"error": f"HTTP {resp.status_code if resp is not None else 'none'}"}
        try:
            payload = resp.json()
        except ValueError:
            return "SOURCE_UNAVAILABLE", {"error": "invalid JSON"}
        if not isinstance(payload, list):
            return "SOURCE_UNAVAILABLE", {"error": "unexpected response shape"}
        by_id: dict[str, dict] = {}
        for entry in payload:
            q = entry.get("query")
            if q is not None and q not in by_id:  # keep the first hit per query
                by_id[q] = {"notfound": True} if entry.get("notfound") else entry
        return "AVAILABLE", by_id

    def annotate_batch(self, variants: list[Variant], progress: Callable[[int], None] | None = None) -> dict[str, EvidenceResult]:
        self.stats = SourceStats()
        breaker = CircuitBreaker(threshold=3, deadline_seconds=self.max_seconds)
        results: dict[str, EvidenceResult] = {}
        pending: dict[tuple[str, str], list[Variant]] = {}   # (assembly, hgvs_id) -> variants needing a live call

        for v in variants:
            if v.genome_build not in _BUILD_TO_ASSEMBLY:
                results[v.variant_id] = not_assessable("myvariant", "genome build unknown/unsupported")
                if progress:
                    progress(1)
                continue
            hgvs_id, err = resolve_hgvs_genomic_id(v)
            if hgvs_id is None:
                results[v.variant_id] = not_assessable("myvariant", err or "could not construct HGVS id")
                if progress:
                    progress(1)
                continue
            assembly = _BUILD_TO_ASSEMBLY[v.genome_build]
            local = get_local(TABLE, f"{assembly}:{hgvs_id}", db_path=self.db_path)
            if local is not None:
                self.stats.hit()
                results[v.variant_id] = EvidenceResult("myvariant", SourceStatus(local["status"]), local["data"],
                                                       retrieved_at=local["retrieved_at"], source_type="local_db")
                if progress:
                    progress(1)
            else:
                pending.setdefault((assembly, hgvs_id), []).append(v)

        if not pending:
            return results

        if is_offline():
            for vs in pending.values():
                for v in vs:
                    results[v.variant_id] = skipped("myvariant", "offline mode: not in cache")
                    self.stats.skip()
                    if progress:
                        progress(1)
            return results

        by_assembly: dict[str, list[str]] = {}
        for assembly, hgvs_id in pending:
            by_assembly.setdefault(assembly, []).append(hgvs_id)
        chunks = [(a, ids[i:i + CHUNK]) for a, ids in by_assembly.items() for i in range(0, len(ids), CHUNK)]

        with httpx.Client(timeout=make_timeout(self.timeout)) as client:
            def do_chunk(chunk: tuple[str, list[str]]) -> dict[str, EvidenceResult]:
                assembly, ids = chunk
                out: dict[str, EvidenceResult] = {}
                n_variants = sum(len(pending[(assembly, i)]) for i in ids)
                if not breaker.allow():
                    reason = breaker.reason() or "skipped"
                    for i in ids:
                        for v in pending[(assembly, i)]:
                            out[v.variant_id] = skipped("myvariant", reason)
                    self.stats.skip(n_variants)
                    return out
                t0 = time.monotonic()
                try:
                    status, payload = self._live_fetch_many(ids, assembly, client)
                except Exception as e:  # noqa: BLE001
                    status, payload = "SOURCE_UNAVAILABLE", {"error": type(e).__name__}
                ok = status == "AVAILABLE"
                self.stats.live(ok, time.monotonic() - t0, n=n_variants)
                breaker.record(ok)
                now = datetime.now(timezone.utc).isoformat()
                for i in ids:
                    if ok:
                        data = payload.get(i, {"notfound": True})
                        save_local(TABLE, f"{assembly}:{i}", "AVAILABLE", data, source_type="live_api", db_path=self.db_path)
                        res = EvidenceResult("myvariant", SourceStatus.AVAILABLE, data, retrieved_at=now, source_type="live_api")
                    else:
                        res = EvidenceResult("myvariant", SourceStatus.UNAVAILABLE, payload, error=payload.get("error"),
                                             retrieved_at=now, source_type="live_api")
                    for v in pending[(assembly, i)]:
                        out[v.variant_id] = res
                if progress:
                    progress(n_variants)
                return out

            for part in run_parallel(chunks, do_chunk, workers=WORKERS):
                results.update(part)
        return results
