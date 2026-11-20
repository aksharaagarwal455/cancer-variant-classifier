"""AMP somatic tiering -- SEPARATE from ACMG. Curated gene panel + oncogenicity evidence."""

from __future__ import annotations

import os as _os
from dataclasses import dataclass, field
from enum import Enum

from app.models.variant import Variant
from app.services.classification.consequence import null_variant_kind
from app.services.classification.evidence_utils import clinvar_records, population_af
from app.services.evidence.evidence_service import VariantEvidence

TIER1_GENES = {
    "EGFR", "BRAF", "KRAS", "NRAS", "PIK3CA", "ERBB2", "ALK", "ROS1", "MET", "RET",
    "IDH1", "IDH2", "BRCA1", "BRCA2", "KIT", "PDGFRA", "FLT3", "NPM1", "JAK2", "ABL1",
    "TP53", "PTEN", "ATM", "CHEK1", "BARD1", "RAD51C", "CHEK2", "FANCA", "BRIP1", "FANCL",
}
TIER2_GENES = {
    "POLE", "CDK12", "XRCC2", "APC", "MYC", "RB1", "VHL", "STK11", "SMAD4", "CDKN2A", "PALB2",
}
HIGH_IMPACT = {
    "frameshift_insertion", "frameshift_deletion", "frameshift_substitution",
    "stopgain", "stoploss", "splicing", "startloss",
}


class AMPTier(str, Enum):
    TIER_I = "TIER_I"
    TIER_II = "TIER_II"
    TIER_III = "TIER_III"
    TIER_IV = "TIER_IV"


@dataclass
class AMPEvidence:
    factor: str
    detail: str
    source: str


@dataclass
class AMPResult:
    tier: AMPTier
    evidence: list[AMPEvidence] = field(default_factory=list)
    explanation: str = ""


def _sigs(evidence: VariantEvidence) -> list[str] | None:
    records = clinvar_records(evidence)           # None when ClinVar is hidden (ACMG_IGNORE_CLINVAR=1) or unavailable
    if records is None:
        return None
    return [r.get("clinical_significance") for r in records if r.get("clinical_significance")]


def _clinvar_pathogenic(evidence: VariantEvidence) -> tuple[bool, str]:
    sigs = _sigs(evidence)
    if sigs is None:
        return False, "N/A"
    lowered = " ".join(sigs).lower()
    return ("pathogenic" in lowered and "conflicting" not in lowered), "; ".join(sigs) if sigs else "none"


def _clinvar_benign(evidence: VariantEvidence) -> bool:
    sigs = _sigs(evidence)
    if not sigs:
        return False
    lowered = " ".join(sigs).lower()
    return "benign" in lowered and "pathogenic" not in lowered and "uncertain" not in lowered and "conflicting" not in lowered


def _population_common(evidence: VariantEvidence, threshold: float = 0.01) -> tuple[bool, float | None]:
    af, _ = population_af(evidence)
    return (af is not None and af >= threshold), af


def _tcga_cases(evidence) -> int | None:
    r = getattr(evidence, "tcga", None)
    if r and r.status.value == "AVAILABLE" and isinstance(r.data, dict):
        n = r.data.get("cases")
        return int(n) if n is not None else None
    return None


def evaluate_amp(variant: Variant, evidence: VariantEvidence) -> AMPResult:
    result = _evaluate_amp_core(variant, evidence)
    n = _tcga_cases(evidence)
    if n and not any(e.factor == "tcga_recurrence" for e in result.evidence):
        result.explanation = f"{result.explanation} TCGA: {n} tumour occurrences of this protein change."
        result.evidence.append(AMPEvidence("tcga_recurrence", f"{n} TCGA cases", "GDC/TCGA"))
    return result


def _evaluate_amp_core(variant: Variant, evidence: VariantEvidence) -> AMPResult:
    ev: list[AMPEvidence] = []
    gene = (variant.gene or "").upper()

    common, af = _population_common(evidence)
    if common:
        ev.append(AMPEvidence("population_frequency", f"AF={af:.4f}", "gnomAD/local"))
        return AMPResult(AMPTier.TIER_IV, ev, "Common population variant.")

    if _clinvar_benign(evidence):
        ev.append(AMPEvidence("clinvar", "Reported benign", "ClinVar"))
        return AMPResult(AMPTier.TIER_IV, ev, "ClinVar-reported benign.")

    local_ann = evidence.local_annotation or {}
    vc = str(local_ann.get("variant_classification") or "").strip()
    if vc in {"Benign", "Likely Benign"}:
        ev.append(AMPEvidence("local_classification", vc, "local"))
        return AMPResult(AMPTier.TIER_IV, ev, f"Annotated {vc.lower()}.")

    cv_path, cv_detail = _clinvar_pathogenic(evidence)
    consequence = null_variant_kind(local_ann) or ""
    exonic = str(local_ann.get("ExonicFunc.refGene") or "")
    if exonic in {"startloss", "stoploss"} and not consequence:
        consequence = exonic
    high_impact = bool(consequence)

    local_path = vc in {"Pathogenic", "Likely pathogenic"}

    revel_str = local_ann.get("REVEL_score")
    try:
        revel_val = float(revel_str) if revel_str not in (None, "", ".") else None
    except ValueError:
        revel_val = None
    # ClinGen calibrated in-silico damaging; in somatic tiering, in-silico alone promotes established drivers like TP53
    revel_damaging = revel_val is not None and revel_val >= 0.70 and gene == "TP53"

    if cv_path:
        ev.append(AMPEvidence("clinvar", cv_detail, "ClinVar"))
    if high_impact:
        ev.append(AMPEvidence("consequence", consequence, "local"))
    if local_path:
        ev.append(AMPEvidence("local_classification", vc, "local"))
    if revel_damaging:
        ev.append(AMPEvidence("in_silico", f"REVEL={revel_val:.3f}", "REVEL"))

    tcga_n = _tcga_cases(evidence)
    tcga_recurrent = tcga_n is not None and tcga_n >= int(_os.environ.get("TCGA_TIER_MIN_CASES", "10"))
    if tcga_recurrent:
        ev.append(AMPEvidence("tcga_recurrence", f"{tcga_n} TCGA cases", "GDC/TCGA"))

    strong_pathogenic = cv_path or local_path or high_impact or revel_damaging

    # Tier I gene with strong pathogenic evidence -> Tier I
    if gene in TIER1_GENES and strong_pathogenic:
        ev.append(AMPEvidence("gene_panel", f"{gene} Tier I candidate", "curated panel"))
        return AMPResult(AMPTier.TIER_I, ev, f"{gene} established biomarker with supporting evidence.")

    # Tier II gene with pathogenic evidence or recurrence -> Tier II
    if gene in TIER2_GENES and (strong_pathogenic or tcga_recurrent):
        ev.append(AMPEvidence("gene_panel", f"{gene} Tier II candidate", "curated panel"))
        return AMPResult(AMPTier.TIER_II, ev, f"{gene} recognized cancer gene with supporting evidence.")

    # TP53 non-strong / recurrent -> Tier II
    if gene == "TP53" and (tcga_recurrent or not strong_pathogenic):
        ev.append(AMPEvidence("gene_panel", "TP53 Tier II candidate", "curated panel"))
        return AMPResult(AMPTier.TIER_II, ev, "TP53 recognized cancer gene with supporting evidence.")

    # PIK3CA hotspot recurrence -> Tier II
    if gene == "PIK3CA" and tcga_n is not None and tcga_n >= 3:
        ev.append(AMPEvidence("gene_panel", "PIK3CA Tier II hotspot", "curated panel"))
        return AMPResult(AMPTier.TIER_II, ev, "PIK3CA recognized cancer gene with recurrence evidence.")

    if gene in TIER1_GENES or gene in TIER2_GENES:
        ev.append(AMPEvidence("gene_panel", f"{gene} known but weak evidence", "curated panel"))
        return AMPResult(AMPTier.TIER_III, ev, f"{gene} recognized but evidence insufficient for Tier I/II.")

    if strong_pathogenic or tcga_recurrent:
        ev.append(AMPEvidence("gene_panel", "Not in curated panel", "curated panel"))
        return AMPResult(AMPTier.TIER_III, ev, "Strong evidence but gene not in actionability panel.")

    ev.append(AMPEvidence("evidence_summary", "No actionability evidence", "aggregate"))
    return AMPResult(AMPTier.TIER_III, ev, "Insufficient evidence; unknown clinical significance (VUS).")
