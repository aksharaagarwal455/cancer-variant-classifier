"""
Evidence aggregation: merges local ANNOVAR annotation + the four sources per variant
(MyVariant, ClinVar, gnomAD, TCGA/GDC), keeping each distinct (never overwritten).

The four sources run CONCURRENTLY (one thread each; inside each source lookups are also parallel),
so wall time is about the slowest single source, not the sum. A source that crashes is reported
and the run continues with that source missing.
"""

from __future__ import annotations

import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import Any

from app import config  # noqa: F401  (loads .env: NCBI_API_KEY etc.)
from app.database.db import get_engine, purge_failures
from app.models.variant import Variant
from app.services.sources.base import EvidenceResult, SourceStatus
from app.services.sources.clinvar import ClinVarClient
from app.services.sources.gnomad import GnomadClient
from app.services.sources.myvariant import MyVariantClient
from app.services.sources.tcga import TcgaClient

LOCAL_ANNOTATION_FIELDS = [
    "CLNSIG", "variant_classification", "AF", "AF_popmax", "AF_raw", "ExAC_ALL",
    "SIFT_pred", "Polyphen2_HVAR_pred", "Polyphen2_HDIV_pred", "MutationTaster_pred",
    "FATHMM_pred", "MetaSVM_pred", "MetaLR_pred", "REVEL_score", "CADD_phred",
    "ExonicFunc.refGene", "Func.refGene", "cosmic103",
]


@dataclass
class VariantEvidence:
    variant_id: str
    local_annotation: dict[str, Any] = field(default_factory=dict)
    myvariant: EvidenceResult | None = None
    clinvar: EvidenceResult | None = None          # EXACT-allele ClinVar (derived from MyVariant) -- used for scoring
    clinvar_position: EvidenceResult | None = None  # NCBI search by position -- information only, never scored
    gnomad: EvidenceResult | None = None
    tcga: EvidenceResult | None = None


def clinvar_from_myvariant(mv: EvidenceResult | None) -> EvidenceResult:
    """Exact-allele ClinVar: MyVariant looks the variant up by its full genomic HGVS id (chrom, pos, ref, alt)."""
    if mv is None:
        return EvidenceResult("clinvar", SourceStatus.UNAVAILABLE, {"error": "MyVariant not queried"},
                              error="MyVariant not queried", source_type="n/a")
    if mv.status != SourceStatus.AVAILABLE:
        return EvidenceResult("clinvar", mv.status, {}, error=mv.error or "MyVariant unavailable",
                              source_type=mv.source_type, retrieved_at=mv.retrieved_at)
    cv = (mv.data or {}).get("clinvar")
    rcv = cv.get("rcv") if isinstance(cv, dict) else None
    rcv = rcv if isinstance(rcv, list) else ([rcv] if isinstance(rcv, dict) else [])
    records = []
    for r in rcv:
        cond = r.get("conditions")
        cond = cond.get("name") if isinstance(cond, dict) else None
        records.append({"clinical_significance": r.get("clinical_significance"), "review_status": r.get("review_status"),
                        "condition": [cond] if cond else [], "accession": r.get("accession"),
                        "last_evaluated": r.get("last_evaluated"), "origin": r.get("origin")})
    data = {"found": bool(records), "records": records, "match": "exact_allele",
            "variation_id": cv.get("variant_id") or cv.get("allele_id") if isinstance(cv, dict) else None}
    return EvidenceResult("clinvar", SourceStatus.AVAILABLE, data, retrieved_at=mv.retrieved_at, source_type=mv.source_type)


def extract_local_annotation(variant: Variant) -> dict[str, Any]:
    return {f: variant.original_input.get(f) for f in LOCAL_ANNOTATION_FIELDS if f in variant.original_input}


class EvidenceService:
    def __init__(self, db_path: str = "data/variant_classifier.db"):
        self.db_path = db_path
        get_engine(db_path)
        self.purged = purge_failures(db_path)   # old versions stored failures; never trust them
        self.clients: dict[str, Any] = {
            "myvariant": MyVariantClient(db_path=db_path),
            "clinvar": ClinVarClient(db_path=db_path),
            "gnomad": GnomadClient(db_path=db_path),
            "tcga": TcgaClient(db_path=db_path),
        }
        self.timings: dict[str, dict] = {}
        self.errors: dict[str, str] = {}

    def gather(self, variants: list[Variant], show_progress: bool = False,
               max_seconds: float | None = None) -> dict[str, VariantEvidence]:
        self.timings, self.errors = {}, {}
        for c in self.clients.values():
            if hasattr(c, "max_seconds"):
                c.max_seconds = max_seconds

        progress_ctx = None
        callbacks: dict[str, Any] = {name: None for name in self.clients}
        if show_progress:
            from rich.progress import BarColumn, Progress, TextColumn, TimeElapsedColumn
            progress_ctx = Progress(TextColumn("{task.description:<10}"), BarColumn(),
                                    TextColumn("{task.completed}/{task.total}"), TimeElapsedColumn())
            progress_ctx.start()
            for name in self.clients:
                tid = progress_ctx.add_task(name, total=len(variants))
                callbacks[name] = (lambda n=1, t=tid: progress_ctx.advance(t, n))

        def run(name: str):
            t0 = time.monotonic()
            try:
                res = self.clients[name].annotate_batch(variants, progress=callbacks[name])
            except Exception as e:  # noqa: BLE001  (one broken source must never break the run)
                res = {}
                self.errors[name] = f"{type(e).__name__}: {e}"
            return name, res, time.monotonic() - t0

        try:
            with ThreadPoolExecutor(max_workers=len(self.clients)) as ex:
                outs = list(ex.map(run, list(self.clients)))
        finally:
            if progress_ctx:
                progress_ctx.stop()

        by_source: dict[str, dict[str, EvidenceResult]] = {}
        for name, res, wall in outs:
            by_source[name] = res
            stats = getattr(self.clients[name], "stats", None)
            self.timings[name] = {**(stats.summary() if stats else {}), "wall_s": round(wall, 2)}

        evidence: dict[str, VariantEvidence] = {}
        for v in variants:
            evidence[v.variant_id] = VariantEvidence(
                variant_id=v.variant_id,
                local_annotation=extract_local_annotation(v),
                myvariant=by_source["myvariant"].get(v.variant_id),
                clinvar=clinvar_from_myvariant(by_source["myvariant"].get(v.variant_id)),
                clinvar_position=by_source["clinvar"].get(v.variant_id),
                gnomad=by_source["gnomad"].get(v.variant_id),
                tcga=by_source["tcga"].get(v.variant_id),
            )
        return evidence
