"""Tests for the correctness fixes: exact-allele ClinVar, exome AF, absent-from-gnomAD, calibrated REVEL, flags."""
from __future__ import annotations

import os

from app.models.variant import GenomeBuild, Variant
from app.services.classification.acmg import evaluate_acmg
from app.services.classification.amp import evaluate_amp
from app.services.classification.evidence_utils import population_af, revel_strength
from app.services.classification.quality import TIER_LABELS, assess_quality
from app.services.evidence.evidence_service import VariantEvidence, clinvar_from_myvariant
from app.services.sources.base import EvidenceResult, SourceStatus


def V(gene="GENE", hgvsp=None):
    return Variant(variant_id="v", original_row_id=1, chromosome="1", position=100, ref="A", alt="T",
                   genome_build=GenomeBuild.GRCH37, gene=gene, hgvsp=hgvsp)


def mv(data, status=SourceStatus.AVAILABLE, st="live_api"):
    return EvidenceResult("myvariant", status, data, source_type=st)


def test_clinvar_is_taken_from_exact_allele_myvariant_data():
    res = clinvar_from_myvariant(mv({"clinvar": {"variant_id": 1, "rcv": [
        {"accession": "RCV1", "clinical_significance": "Pathogenic", "review_status": "criteria provided, single submitter",
         "conditions": {"name": "X"}}]}}))
    assert res.status == SourceStatus.AVAILABLE and res.data["match"] == "exact_allele"
    assert res.data["records"][0]["clinical_significance"] == "Pathogenic"


def test_variant_not_in_clinvar_is_found_false_not_unavailable():
    res = clinvar_from_myvariant(mv({"notfound": True}))
    assert res.status == SourceStatus.AVAILABLE and res.data["found"] is False and res.data["records"] == []


def test_clinvar_propagates_skipped_when_myvariant_skipped():
    res = clinvar_from_myvariant(mv({}, SourceStatus.UNAVAILABLE, "skipped"))
    assert res.status == SourceStatus.UNAVAILABLE and res.source_type == "skipped"


def test_population_af_uses_exome_when_genome_missing():
    gn = EvidenceResult("gnomad", SourceStatus.AVAILABLE, {"found": True, "variant": {"genome": None, "exome": {"af": 0.07}}})
    ev = VariantEvidence("v", {}, gnomad=gn)
    assert population_af(ev) == (0.07, "gnomAD exome")


def test_population_af_takes_the_highest_source():
    gn = EvidenceResult("gnomad", SourceStatus.AVAILABLE, {"found": True, "variant": {"genome": {"af": 0.0001}, "exome": {"af": 0.02}}})
    af, src = population_af(VariantEvidence("v", {"AF": "0.001"}, gnomad=gn))
    assert af == 0.02 and "exome" in src


def test_gnomad_not_found_without_file_frequency_is_absent():
    gn = EvidenceResult("gnomad", SourceStatus.AVAILABLE, {"found": False, "errors": [{"message": "Variant not found"}]})
    assert population_af(VariantEvidence("v", {}, gnomad=gn))[0] == 0.0


def test_gnomad_invalid_id_error_is_not_treated_as_absent():
    gn = EvidenceResult("gnomad", SourceStatus.AVAILABLE, {"found": False, "errors": [{"message": "Invalid variant ID"}]})
    assert population_af(VariantEvidence("v", {}, gnomad=gn))[0] is None


def test_revel_calibration_thresholds():
    assert revel_strength("0.95") == ("pathogenic", "STRONG")
    assert revel_strength("0.80") == ("pathogenic", "MODERATE")
    assert revel_strength("0.70") == ("pathogenic", "SUPPORTING")
    assert revel_strength("0.50") is None
    assert revel_strength("0.25") == ("benign", "SUPPORTING")
    assert revel_strength(".") is None


def test_high_revel_plus_rare_variant_reaches_likely_pathogenic(monkeypatch):
    monkeypatch.delenv("ACMG_REVEL_PP3", raising=False)
    ev = VariantEvidence("v", {"ExonicFunc.refGene": "nonsynonymous SNV", "AF": "0", "REVEL_score": "0.95"})
    r = evaluate_acmg(V(), ev)
    assert r.score == 6 and r.classification.value == "LIKELY_PATHOGENIC"   # PM2 2 + PP3 strong 4


def test_revel_switch_turns_calibration_off(monkeypatch):
    monkeypatch.setenv("ACMG_REVEL_PP3", "0")
    ev = VariantEvidence("v", {"ExonicFunc.refGene": "nonsynonymous SNV", "AF": "0", "REVEL_score": "0.95"})
    assert evaluate_acmg(V(), ev).score == 2                               # PM2 only


def test_ignore_clinvar_also_hides_it_from_amp(monkeypatch):
    cv = EvidenceResult("clinvar", SourceStatus.AVAILABLE, {"records": [{"clinical_significance": "Pathogenic"}]})
    ev = VariantEvidence("v", {"ExonicFunc.refGene": "nonsynonymous SNV", "AF": "0"}, clinvar=cv)
    monkeypatch.delenv("ACMG_IGNORE_CLINVAR", raising=False)
    assert evaluate_amp(V("BRCA1"), ev).tier.value == "TIER_I"
    monkeypatch.setenv("ACMG_IGNORE_CLINVAR", "1")
    assert evaluate_amp(V("BRCA1"), ev).tier.value != "TIER_I"


def test_quality_flags_missing_sources():
    ev = VariantEvidence("v", {"AF": "0.001"})
    comp, flags = assess_quality(ev, "VUS", "TIER_III")
    assert comp == "LOW" and any("gnomAD" in f for f in flags) and any("MyVariant" in f for f in flags)


def test_quality_flags_position_only_clinvar_entry():
    ok = EvidenceResult("x", SourceStatus.AVAILABLE, {})
    cv = EvidenceResult("clinvar", SourceStatus.AVAILABLE, {"found": False, "records": []})
    pos = EvidenceResult("clinvar", SourceStatus.AVAILABLE, {"found": True, "records": [{"clinical_significance": "Pathogenic"}]})
    ev = VariantEvidence("v", {"REVEL_score": "0.5"}, myvariant=ok, clinvar=cv, clinvar_position=pos, gnomad=ok)
    comp, flags = assess_quality(ev, "VUS", "TIER_III")
    assert comp == "HIGH" and any("not for this exact allele" in f for f in flags)


def test_quality_flags_acmg_amp_disagreement():
    _, flags = assess_quality(VariantEvidence("v", {}), "PATHOGENIC", "TIER_IV")
    assert any("disagree" in f for f in flags)


def test_every_tier_has_a_label():
    assert set(TIER_LABELS) == {"TIER_I", "TIER_II", "TIER_III", "TIER_IV"}


def test_pvs1_outside_known_lof_genes_is_flagged():
    _, flags = assess_quality(VariantEvidence("v", {}), "PATHOGENIC", "TIER_III", pvs1=True, gene="RAD54L")
    assert any("loss-of-function" in f for f in flags)
    _, flags = assess_quality(VariantEvidence("v", {}), "PATHOGENIC", "TIER_I", pvs1=True, gene="BRCA2")
    assert not any("loss-of-function" in f for f in flags)


# ---------------------------------------------------------------- API
def _client():
    from fastapi.testclient import TestClient
    from app.api.main import app
    return TestClient(app)


def test_api_offline_flag_does_not_leak_into_the_server_process(monkeypatch):
    monkeypatch.delenv("EVIDENCE_OFFLINE", raising=False)
    with open("data/input/S11_labeled_subset.tsv", "rb") as f:
        r = _client().post("/api/classify?limit=3&offline=true", files={"file": ("x.tsv", f, "text/plain")})
    assert r.status_code == 200 and r.json()["total_variants"] == 3
    assert os.environ.get("EVIDENCE_OFFLINE") is None
    row = r.json()["results"][0]
    assert row["amp_tier_label"] and row["evidence_completeness"] in ("HIGH", "MEDIUM", "LOW") and "review_flags" in row


def test_api_returns_a_clear_error_for_a_file_with_no_variants():
    import io
    r = _client().post("/api/classify", files={"file": ("bad.tsv", io.BytesIO(b"foo\tbar\n1\t2\n"), "text/plain")})
    assert r.status_code == 422 and "CHROM" in r.json()["detail"]
