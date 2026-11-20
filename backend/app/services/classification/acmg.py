"""
ACMG classification -- documented subset (PVS1, PM2, PP3, PP5, BA1, BS1, BP4, BP6)
combined via ClinGen/Tavtigian 2018 Bayesian points (PMID 29300386).
NOT the full 28-criteria ACMG guideline -- disclosed, not hidden.

PVS1 is not applied when ClinVar reports the variant benign/likely benign
(curated evidence contradicts the null-variant assumption).
Set ACMG_IGNORE_CLINVAR=1 to skip ClinVar entirely (for independent validation).
"""

from __future__ import annotations
import os
from dataclasses import dataclass, field
from enum import Enum
from app.models.variant import Variant
from app.services.classification.consequence import null_variant_kind
from app.services.classification.evidence_utils import clinvar_ignored, clinvar_records, population_af, revel_strength
from app.services.evidence.evidence_service import VariantEvidence

DAMAGING_PRED_VALUES = {"D", "P", "A", "H"}
IN_SILICO_COLUMNS = ["SIFT_pred", "Polyphen2_HVAR_pred", "Polyphen2_HDIV_pred", "MutationTaster_pred", "FATHMM_pred", "MetaSVM_pred", "MetaLR_pred"]
BA1_THRESHOLD, BS1_THRESHOLD, PM2_THRESHOLD = 0.05, 0.01, 0.0001
POINTS = {"VERY_STRONG": 8, "STRONG": 4, "MODERATE": 2, "SUPPORTING": 1}


class ACMGClassification(str, Enum):
    PATHOGENIC = "PATHOGENIC"
    LIKELY_PATHOGENIC = "LIKELY_PATHOGENIC"
    VUS = "VUS"
    LIKELY_BENIGN = "LIKELY_BENIGN"
    BENIGN = "BENIGN"


@dataclass
class Criterion:
    code: str
    status: str
    strength: str | None
    direction: str
    evidence: str
    source: str


@dataclass
class ACMGResult:
    classification: ACMGClassification
    score: int
    criteria: list[Criterion] = field(default_factory=list)
    explanation: str = ""


def _get_population_af(evidence: VariantEvidence) -> tuple[float | None, str]:
    return population_af(evidence)


def _count_in_silico(evidence: VariantEvidence) -> tuple[int, int, list[str]]:
    damaging, benign, used = 0, 0, []
    for col in IN_SILICO_COLUMNS:
        val = evidence.local_annotation.get(col)
        if not val:
            continue
        v = val.strip().upper()
        used.append(f"{col}={v}")
        if v in DAMAGING_PRED_VALUES:
            damaging += 1
        elif v in ("T", "B", "N"):
            benign += 1
    cadd_raw = evidence.local_annotation.get("CADD_phred")
    if cadd_raw not in (None, ""):
        try:
            if float(cadd_raw) >= 20:
                damaging += 1
                used.append(f"CADD_phred={cadd_raw}")
        except (TypeError, ValueError):
            pass
    return damaging, benign, used


def _clinvar_signal(evidence: VariantEvidence) -> tuple[str | None, str]:
    if clinvar_ignored():
        return None, "ClinVar ignored (ACMG_IGNORE_CLINVAR=1)"
    records = clinvar_records(evidence)
    if records is None:
        return None, "ClinVar not available"
    sigs = [r.get("clinical_significance") for r in records if r.get("clinical_significance")]
    if not sigs:
        return None, "no ClinVar significance recorded"
    lowered = " ".join(sigs).lower()
    if "pathogenic" in lowered and "benign" not in lowered and "conflicting" not in lowered:
        return "pathogenic", "; ".join(sigs)
    if "benign" in lowered and "pathogenic" not in lowered and "conflicting" not in lowered:
        return "benign", "; ".join(sigs)
    return None, f"ambiguous: {'; '.join(sigs)}"


def evaluate_acmg(variant: Variant, evidence: VariantEvidence) -> ACMGResult:
    criteria: list[Criterion] = []
    score = 0
    ba1 = False

    cv_signal, cv_detail = _clinvar_signal(evidence)

    null_kind = null_variant_kind(evidence.local_annotation)   # exonic null OR canonical splice site (switchable)
    if null_kind and cv_signal == "benign":
        criteria.append(Criterion("PVS1", "NOT_SUPPORTED", None, "pathogenic",
                                  f"Null-type consequence ({null_kind}) but ClinVar reports benign/likely benign; PVS1 not applied",
                                  "local+ClinVar"))
    elif null_kind:
        criteria.append(Criterion("PVS1", "SUPPORTED", "VERY_STRONG", "pathogenic", f"Null variant: {null_kind}", "local"))
        score += POINTS["VERY_STRONG"]
    else:
        consequence = (evidence.local_annotation.get("ExonicFunc.refGene") or "").strip().lower()
        criteria.append(Criterion("PVS1", "NOT_SUPPORTED", None, "pathogenic", f"Not null: {consequence or 'unknown'}", "local"))

    af, af_source = _get_population_af(evidence)
    if af is None:
        criteria.append(Criterion("PM2", "NOT_ASSESSABLE", None, "pathogenic", "No AF data", af_source))
    elif af >= BA1_THRESHOLD:
        criteria.append(Criterion("BA1", "SUPPORTED", "STAND_ALONE", "benign", f"AF={af:.5f}", af_source))
        ba1 = True
    elif af >= BS1_THRESHOLD:
        criteria.append(Criterion("BS1", "SUPPORTED", "STRONG", "benign", f"AF={af:.5f}", af_source))
        score -= POINTS["STRONG"]
    elif af < PM2_THRESHOLD:
        criteria.append(Criterion("PM2", "SUPPORTED", "MODERATE", "pathogenic", f"AF={af:.6f} rare", af_source))
        score += POINTS["MODERATE"]
    else:
        criteria.append(Criterion("PM2", "NOT_SUPPORTED", None, "pathogenic", f"AF={af:.6f}", af_source))

    revel = revel_strength(evidence.local_annotation.get("REVEL_score")) if os.environ.get("ACMG_REVEL_PP3", "1") != "0" else None
    dmg, ben, used = _count_in_silico(evidence)
    if revel and revel[0] == "pathogenic":
        criteria.append(Criterion("PP3", "SUPPORTED", revel[1], "pathogenic",
                                  f"REVEL={evidence.local_annotation.get('REVEL_score')} (ClinGen-calibrated {revel[1].lower()})", "local REVEL"))
        score += POINTS[revel[1]]
    elif revel and revel[0] == "benign":
        criteria.append(Criterion("BP4", "SUPPORTED", revel[1], "benign",
                                  f"REVEL={evidence.local_annotation.get('REVEL_score')} (ClinGen-calibrated {revel[1].lower()})", "local REVEL"))
        score -= POINTS[revel[1]]
    elif dmg >= 3:
        criteria.append(Criterion("PP3", "SUPPORTED", "SUPPORTING", "pathogenic", f"{dmg} damaging", "local"))
        score += POINTS["SUPPORTING"]
    elif ben >= 3:
        criteria.append(Criterion("BP4", "SUPPORTED", "SUPPORTING", "benign", f"{ben} benign", "local"))
        score -= POINTS["SUPPORTING"]
    else:
        criteria.append(Criterion("PP3/BP4", "NOT_SUPPORTED", None, "pathogenic", f"{dmg}D/{ben}B of {len(used)}", "local"))

    if cv_signal == "pathogenic":
        criteria.append(Criterion("PP5", "SUPPORTED", "SUPPORTING", "pathogenic", cv_detail, "ClinVar"))
        score += POINTS["SUPPORTING"]
    elif cv_signal == "benign":
        criteria.append(Criterion("BP6", "SUPPORTED", "SUPPORTING", "benign", cv_detail, "ClinVar"))
        score -= POINTS["SUPPORTING"]
    else:
        criteria.append(Criterion("PP5/BP6", "NOT_ASSESSABLE", None, "pathogenic", cv_detail, "ClinVar"))

    if ba1:
        classification, explanation = ACMGClassification.BENIGN, "BA1 stand-alone rule."
    elif score >= 10:
        classification, explanation = ACMGClassification.PATHOGENIC, f"Score {score} >= 10."
    elif score >= 6:
        classification, explanation = ACMGClassification.LIKELY_PATHOGENIC, f"Score {score} in [6,9]."
    elif score >= 0:
        classification, explanation = ACMGClassification.VUS, f"Score {score} in [0,5]."
    elif score >= -6:
        classification, explanation = ACMGClassification.LIKELY_BENIGN, f"Score {score} in [-6,-1]."
    else:
        classification, explanation = ACMGClassification.BENIGN, f"Score {score} <= -7."

    return ACMGResult(classification=classification, score=score, criteria=criteria, explanation=explanation)
