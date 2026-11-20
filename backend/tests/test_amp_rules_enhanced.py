"""Tests for enhanced AMP somatic tiering rules and expanded gene panels."""

from __future__ import annotations

import os

from app.models.variant import GenomeBuild, Variant
from app.services.classification.amp import AMPTier, evaluate_amp
from app.services.evidence.evidence_service import VariantEvidence
from app.services.sources.base import EvidenceResult, SourceStatus


def make_var(gene="TP53", **local):
    v = Variant(
        variant_id=f"chr17:100:A:T:{gene}",
        original_row_id=0,
        chromosome="17",
        position=100,
        ref="A",
        alt="T",
        genome_build=GenomeBuild.GRCH37,
        gene=gene,
    )
    return v, VariantEvidence(variant_id=v.variant_id, local_annotation=local)


def test_amp_tp53_pathogenic_null_is_tier_i():
    v, e = make_var("TP53", **{"ExonicFunc.refGene": "stopgain", "AF": "0.00001"})
    res = evaluate_amp(v, e)
    assert res.tier == AMPTier.TIER_I
    assert any("Tier I candidate" in ev.detail for ev in res.evidence)


def test_amp_tp53_revel_damaging_is_tier_i():
    v, e = make_var("TP53", **{"ExonicFunc.refGene": "nonsynonymous_SNV", "REVEL_score": "0.85", "AF": "0.00001"})
    res = evaluate_amp(v, e)
    assert res.tier == AMPTier.TIER_I


def test_amp_tp53_vus_is_tier_ii():
    v, e = make_var("TP53", **{"ExonicFunc.refGene": "nonsynonymous_SNV", "REVEL_score": "0.30", "AF": "0.00001"})
    res = evaluate_amp(v, e)
    assert res.tier == AMPTier.TIER_II
    assert "recognized cancer gene" in res.explanation


def test_amp_hrr_panel_genes_pathogenic_are_tier_i():
    for g in ["PTEN", "ATM", "RAD51C", "CHEK2", "FANCA", "BRIP1", "BARD1", "FANCL", "CHEK1"]:
        v, e = make_var(g, **{"variant_classification": "Pathogenic", "AF": "0.00001"})
        assert evaluate_amp(v, e).tier == AMPTier.TIER_I, f"{g} should be Tier I when pathogenic"


def test_amp_tier2_genes_null_variant_are_tier_ii():
    for g in ["POLE", "CDK12", "XRCC2"]:
        v, e = make_var(g, **{"ExonicFunc.refGene": "frameshift_deletion", "AF": "0.00001"})
        assert evaluate_amp(v, e).tier == AMPTier.TIER_II, f"{g} should be Tier II with null variant"


def test_amp_pik3ca_recurrent_hotspot_is_tier_ii_when_weak():
    v, e = make_var("PIK3CA", **{"ExonicFunc.refGene": "nonsynonymous_SNV", "AF": "0.00001"})
    e.tcga = EvidenceResult("tcga", SourceStatus.AVAILABLE, {"cases": 4})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_II


def test_amp_local_annotated_benign_is_tier_iv():
    v, e = make_var("BRCA1", **{"variant_classification": "Benign", "AF": "0.00001"})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_IV


def test_amp_non_panel_gene_pathogenic_is_tier_iii():
    # RAD50 and FANCD2 are not in Tier I or II panels; even if null/pathogenic, they remain Tier III (VUS in clinical guidelines)
    v, e = make_var("RAD50", **{"ExonicFunc.refGene": "stopgain", "AF": "0.00001"})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_III
