from app.models.variant import GenomeBuild, Variant
from app.services.classification.acmg import ACMGClassification, evaluate_acmg
from app.services.classification.amp import AMPTier, evaluate_amp
from app.services.evidence.evidence_service import VariantEvidence
from app.services.sources.base import EvidenceResult, SourceStatus


def make(gene="BRCA2", **local):
    v = Variant(variant_id="chr13:1:A:T", original_row_id=0, chromosome="13", position=1, ref="A", alt="T",
                genome_build=GenomeBuild.GRCH37, gene=gene)
    return v, VariantEvidence(variant_id=v.variant_id, local_annotation=local)


def clinvar(sig):
    return EvidenceResult("clinvar", SourceStatus.AVAILABLE, {"records": [{"clinical_significance": sig}]})


def test_frameshift_rare_is_pathogenic():
    v, e = make(**{"ExonicFunc.refGene": "frameshift_deletion", "AF": "0.00001"})
    r = evaluate_acmg(v, e)
    assert r.classification == ACMGClassification.PATHOGENIC and r.score == 10


def test_splice_site_counts_as_null_by_default(monkeypatch):
    monkeypatch.delenv("ACMG_SPLICE_PVS1", raising=False)
    v, e = make(**{"Func.refGene": "splicing", "ExonicFunc.refGene": "", "AF": "0.00001"})
    assert evaluate_acmg(v, e).classification == ACMGClassification.PATHOGENIC


def test_splice_rule_can_be_switched_off(monkeypatch):
    monkeypatch.setenv("ACMG_SPLICE_PVS1", "0")
    v, e = make(**{"Func.refGene": "splicing", "ExonicFunc.refGene": "", "AF": "0.00001"})
    r = evaluate_acmg(v, e)
    assert r.classification == ACMGClassification.VUS and r.score == 2


def test_exonic_splicing_missense_is_not_null():
    v, e = make(**{"Func.refGene": "exonic;splicing", "ExonicFunc.refGene": "nonsynonymous_SNV", "AF": "0.00001"})
    assert evaluate_acmg(v, e).score == 2


def test_clinvar_benign_blocks_pvs1():
    v, e = make(**{"ExonicFunc.refGene": "stopgain", "AF": "0.00001"})
    e.clinvar = clinvar("Benign")
    r = evaluate_acmg(v, e)
    assert not any(c.code == "PVS1" and c.status == "SUPPORTED" for c in r.criteria)


def test_clinvar_hidden_when_ignored(monkeypatch):
    monkeypatch.setenv("ACMG_IGNORE_CLINVAR", "1")
    v, e = make(**{"AF": "0.00001"})
    e.clinvar = clinvar("Pathogenic")
    assert not any(c.code == "PP5" and c.status == "SUPPORTED" for c in evaluate_acmg(v, e).criteria)


def test_ba1_common_variant_is_benign():
    v, e = make(**{"AF": "0.2", "ExonicFunc.refGene": "stopgain"})
    assert evaluate_acmg(v, e).classification == ACMGClassification.BENIGN


def test_amp_tier_i_for_established_gene_with_null_variant():
    v, e = make("BRCA2", **{"ExonicFunc.refGene": "frameshift_insertion", "AF": "0.00001"})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_I


def test_amp_common_variant_is_tier_iv():
    v, e = make("BRCA2", **{"AF": "0.3"})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_IV


def test_amp_unknown_gene_without_evidence_is_tier_iii():
    v, e = make("FOO1", **{"AF": "0.00001"})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_III


def test_amp_tcga_recurrence_promotes_tier_ii_gene():
    v, e = make("TP53", **{"AF": "0.00001"})
    e.tcga = EvidenceResult("tcga", SourceStatus.AVAILABLE, {"cases": 50})
    assert evaluate_amp(v, e).tier == AMPTier.TIER_II
