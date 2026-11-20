"""FastAPI app: upload a TSV, classify it, return JSON for the browser UI. DB-first evidence."""

from __future__ import annotations

import io
import os
import tempfile
import threading
import time
from pathlib import Path

import pandas as pd
from fastapi import Body, FastAPI, File, HTTPException, UploadFile
from fastapi.responses import FileResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from app import config
from app.services.classification.acmg import evaluate_acmg
from app.services.classification.amp import evaluate_amp
from app.services.classification.evidence_utils import population_af
from app.services.classification.quality import TIER_LABELS, assess_quality
from app.services.evidence.evidence_service import EvidenceService
from app.services.normalization.coordinates import apply_genome_build
from app.services.parser.tsv_parser import ParseStatus, parse_tsv

app = FastAPI(title="Variant Classifier")
STATIC_DIR = Path(__file__).resolve().parent.parent / "static"
_DEFAULT_DB = Path(__file__).resolve().parent.parent.parent.parent / "data" / "variant_classifier.db"
if not _DEFAULT_DB.exists():
    _ALT_DB = Path("data/variant_classifier.db").resolve()
    if _ALT_DB.exists():
        _DEFAULT_DB = _ALT_DB
DB_PATH = os.environ.get("DB_PATH") or str(_DEFAULT_DB)


@app.get("/")
def index():
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/health")
def health():
    return {"status": "ok"}


_RUN_LOCK = threading.Lock()


@app.post("/api/classify")
def classify_endpoint(file: UploadFile = File(...), limit: int = 30, offline: bool = False, max_seconds: int = 120,
                      genome_build: str | None = None):
    """Plain `def` on purpose: FastAPI runs it in a worker thread, so slow lookups never freeze the server."""
    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tmp:
        tmp.write(file.file.read())
        tmp_path = tmp.name
    try:
        with _RUN_LOCK:
            previous = os.environ.get("EVIDENCE_OFFLINE")
            if offline:
                os.environ["EVIDENCE_OFFLINE"] = "1"
            try:
                return _classify(tmp_path, limit, max_seconds, genome_build)
            finally:
                if previous is None:
                    os.environ.pop("EVIDENCE_OFFLINE", None)
                else:
                    os.environ["EVIDENCE_OFFLINE"] = previous
    finally:
        Path(tmp_path).unlink(missing_ok=True)


def _classify(tmp_path: str, limit: int, max_seconds: int, genome_build: str | None = None) -> dict:
    variants, parse_report = parse_tsv(tmp_path)
    ok_variants = [v for v in variants if v.parse_status == ParseStatus.OK]
    if not ok_variants:
        raise HTTPException(status_code=422, detail="No valid variants found. The file needs CHROM, POS, REF, ALT and Gene columns (ANNOVAR-style TSV).")

    try:
        ok_variants, build_evidence = apply_genome_build(ok_variants, forced_build=genome_build)
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    warnings: list[str] = []
    if ok_variants[0].genome_build.name == "UNKNOWN":
        warnings.append("Genome build could not be detected from this file, so MyVariant and gnomAD lookups were skipped and "
                        "results rely on the annotation columns in the file. Choose hg19 or hg38 in 'Genome build' and run again.")
    if build_evidence.get("warning"):
        warnings.append(str(build_evidence["warning"]))
    if limit > 0:
        ok_variants = ok_variants[:limit]

    svc = EvidenceService(db_path=DB_PATH)
    t_ev = time.monotonic()
    evidence_map = svc.gather(ok_variants, max_seconds=max_seconds)
    evidence_seconds = round(time.monotonic() - t_ev, 2)

    results = []
    tier_counts: dict[str, int] = {}
    completeness_counts: dict[str, int] = {}
    for v in ok_variants:
        ev = evidence_map[v.variant_id]
        acmg_result = evaluate_acmg(v, ev)
        amp_result = evaluate_amp(v, ev)
        pvs1 = any(c.code == "PVS1" and c.status == "SUPPORTED" for c in acmg_result.criteria)
        completeness, flags = assess_quality(ev, acmg_result.classification.value, amp_result.tier.value, pvs1=pvs1, gene=v.gene)
        af, af_src = population_af(ev)
        cv_recs = ev.clinvar.data.get("records", []) if ev.clinvar and ev.clinvar.status.value == "AVAILABLE" else []
        tier_counts[amp_result.tier.value] = tier_counts.get(amp_result.tier.value, 0) + 1
        completeness_counts[completeness] = completeness_counts.get(completeness, 0) + 1
        results.append({
            "variant_id": v.variant_id,
            "gene": v.gene,
            "chromosome": v.chromosome,
            "position": v.position,
            "ref": v.ref,
            "alt": v.alt,
            "genome_build": v.genome_build.value,
            "acmg_classification": acmg_result.classification.value,
            "acmg_score": acmg_result.score,
            "acmg_criteria": "; ".join(f"{c.code}={c.status}" for c in acmg_result.criteria if c.status == "SUPPORTED"),
            "acmg_explanation": acmg_result.explanation,
            "amp_tier": amp_result.tier.value,
            "amp_tier_label": TIER_LABELS[amp_result.tier.value],
            "amp_explanation": amp_result.explanation,
            "evidence_completeness": completeness,
            "review_flags": flags,
            "population_af": af,
            "population_af_source": af_src,
            "clinvar_exact_allele": sorted({r.get("clinical_significance") for r in cv_recs if r.get("clinical_significance")}),
            "clinvar_source": ev.clinvar.source_type if ev.clinvar else "N/A",
            "gnomad_source": ev.gnomad.source_type if ev.gnomad else "N/A",
            "myvariant_source": ev.myvariant.source_type if ev.myvariant else "N/A",
            "tcga_source": ev.tcga.source_type if ev.tcga else "N/A",
            "tcga_cases": (ev.tcga.data or {}).get("cases") if ev.tcga and ev.tcga.status.value == "AVAILABLE" else None,
        })

    for name, err in svc.errors.items():
        warnings.append(f"Source '{name}' crashed and was skipped: {err}")
    return {
        "warnings": warnings,
        "total_variants": len(ok_variants),
        "total_parsed": parse_report.total_rows,
        "invalid": parse_report.invalid,
        "genome_build": ok_variants[0].genome_build.value,
        "tier_counts": tier_counts,
        "completeness_counts": completeness_counts,
        "evidence_seconds": evidence_seconds,
        "source_timings": svc.timings,
        "results": results,
    }


class ExportPayload(BaseModel):
    results: list[dict]
    format: str = "xlsx"


@app.post("/api/export")
def export_results(payload: ExportPayload):
    rows = payload.results
    if not rows:
        raise HTTPException(status_code=400, detail="No results provided to export.")

    export_rows = []
    for r in rows:
        cv_sig = r.get("clinvar_exact_allele", [])
        if isinstance(cv_sig, list):
            cv_sig = "; ".join(cv_sig)
        flags = r.get("review_flags", [])
        if isinstance(flags, list):
            flags = " | ".join(flags)
        export_rows.append({
            "Gene": r.get("gene", ""),
            "Variant_ID": r.get("variant_id", ""),
            "Chromosome": r.get("chromosome", ""),
            "Position": r.get("position", ""),
            "Ref": r.get("ref", ""),
            "Alt": r.get("alt", ""),
            "Genome_Build": r.get("genome_build", ""),
            "AMP_Tier": r.get("amp_tier", ""),
            "AMP_Tier_Label": r.get("amp_tier_label", ""),
            "AMP_Rationale": r.get("amp_explanation", ""),
            "ACMG_Classification": r.get("acmg_classification", ""),
            "ACMG_Score": r.get("acmg_score", ""),
            "ACMG_Criteria": r.get("acmg_criteria", ""),
            "ACMG_Explanation": r.get("acmg_explanation", ""),
            "ClinVar": cv_sig,
            "gnomAD_AF": r.get("population_af", ""),
            "TCGA_Cases": r.get("tcga_cases", ""),
            "Evidence_Completeness": r.get("evidence_completeness", ""),
            "Review_Flags": flags,
        })

    df = pd.DataFrame(export_rows)
    fmt = payload.format.lower().strip()
    if fmt == "xlsx":
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as writer:
            df.to_excel(writer, index=False, sheet_name="Classified_Variants")
        buf.seek(0)
        return StreamingResponse(
            buf,
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            headers={"Content-Disposition": "attachment; filename=classified_variants.xlsx"},
        )
    else:
        buf = io.StringIO()
        df.to_csv(buf, sep="\t", index=False)
        return StreamingResponse(
            io.BytesIO(buf.getvalue().encode("utf-8")),
            media_type="text/tab-separated-values",
            headers={"Content-Disposition": "attachment; filename=classified_variants.tsv"},
        )


app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
