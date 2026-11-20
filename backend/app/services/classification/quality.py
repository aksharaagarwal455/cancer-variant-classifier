"""Tier labels and an honest 'how complete was the evidence?' indicator for every variant.

evidence_completeness is NOT a probability that the call is right. It only says how many of the
four evidence streams actually reached the engines for this variant (a call made without gnomAD
or ClinVar is weaker than one made with them), and review_flags says exactly what was missing.
"""

from __future__ import annotations

from app.services.classification.acmg import IN_SILICO_COLUMNS
from app.services.classification.evidence_utils import clinvar_ignored, population_af

TIER_LABELS = {
    "TIER_I": "Tier I - strong clinical significance (candidate; confirm therapy evidence for the tumour type)",
    "TIER_II": "Tier II - potential clinical significance",
    "TIER_III": "Tier III - unknown significance (VUS)",
    "TIER_IV": "Tier IV - benign / likely benign",
}


# Genes where loss of function is an established disease mechanism (PVS1 is only valid there). Curated, short, editable.
LOF_MECHANISM_GENES = {"BRCA1", "BRCA2", "TP53", "PTEN", "APC", "RB1", "ATM", "PALB2", "CHEK2", "BARD1", "RAD51C",
                       "RAD51D", "STK11", "SMAD4", "CDKN2A", "VHL", "MLH1", "MSH2", "MSH6", "PMS2", "NF1", "NF2", "CDH1"}


def assess_quality(evidence, acmg_class: str, amp_tier: str, pvs1: bool = False,
                   gene: str | None = None) -> tuple[str, list[str]]:
    flags: list[str] = []
    points, possible = 0, 0

    # 1) population frequency from a real gnomAD answer
    possible += 1
    gn = evidence.gnomad
    if gn is not None and gn.status.value == "AVAILABLE":
        points += 1
    else:
        af, src = population_af(evidence)
        flags.append("gnomAD unavailable/skipped; frequency from input-file columns" if af is not None
                     else "no population frequency available")

    # 2) MyVariant (carries exact-allele ClinVar, CADD, dbNSFP)
    possible += 1
    mv = evidence.myvariant
    if mv is not None and mv.status.value == "AVAILABLE":
        points += 1
    else:
        flags.append("MyVariant unavailable/not assessable (no exact-allele ClinVar)")

    # 3) exact-allele ClinVar (not counted when deliberately hidden for validation)
    if not clinvar_ignored():
        possible += 1
        cv = evidence.clinvar
        if cv is not None and cv.status.value == "AVAILABLE":
            points += 1
        else:
            flags.append("ClinVar exact-allele data unavailable")
        pos = getattr(evidence, "clinvar_position", None)
        if (cv is not None and cv.status.value == "AVAILABLE" and not cv.data.get("found")
                and pos is not None and pos.status.value == "AVAILABLE" and pos.data.get("found")):
            flags.append("ClinVar has an entry at this position but not for this exact allele (ignored)")

    # 4) in-silico evidence
    possible += 1
    la = evidence.local_annotation
    n_pred = sum(1 for c in IN_SILICO_COLUMNS if (la.get(c) or "").strip() not in ("", ".", "NA"))
    has_revel = (la.get("REVEL_score") or "").strip() not in ("", ".", "NA")
    if n_pred >= 3 or has_revel:
        points += 1
    else:
        flags.append("few in-silico predictors available")

    # PVS1 is only valid when loss of function causes disease in that gene
    if pvs1 and (gene or "").upper() not in LOF_MECHANISM_GENES:
        flags.append("PVS1 applied but loss-of-function as the disease mechanism is not verified for this gene")

    # 5) cross-framework disagreement
    pathogenic = acmg_class in ("PATHOGENIC", "LIKELY_PATHOGENIC")
    benign = acmg_class in ("BENIGN", "LIKELY_BENIGN")
    if (pathogenic and amp_tier == "TIER_IV") or (benign and amp_tier in ("TIER_I", "TIER_II")):
        flags.append("ACMG and AMP disagree: manual review needed")

    ratio = points / possible
    completeness = "HIGH" if ratio == 1 else "MEDIUM" if ratio >= 0.75 else "LOW"
    return completeness, flags
