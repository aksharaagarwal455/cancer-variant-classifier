"""Write rules/acmg_rules.json and rules/amp_rules.json from the constants in the code.

  python scripts/export_rules.py          # (re)write the files
  python scripts/export_rules.py --check  # exit 1 if the files no longer match the code

A unit test runs the --check logic, so the rules folder can never silently drift from the engines.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app.services.classification import acmg, amp
from app.services.classification import evidence_utils as eu
from app.services.classification.consequence import EXONIC_NULL
from app.services.classification.quality import LOF_MECHANISM_GENES


def current_rules() -> dict[str, dict]:
    return {
        "acmg_rules.json": {
            "framework": "ClinGen / Tavtigian 2018 Bayesian points (PMID 29300386)",
            "implemented_criteria": ["PVS1", "PM2", "PP3", "PP5", "BA1", "BS1", "BP4", "BP6"],
            "not_implemented": "about 20 of the 28 ACMG/AMP criteria (e.g. PS1, PS3, PM1, PM5, PP1, PP2, PP4, BS2-BS4, BP1-BP3, BP5, BP7)",
            "points": acmg.POINTS,
            "allele_frequency_thresholds": {"BA1": acmg.BA1_THRESHOLD, "BS1": acmg.BS1_THRESHOLD, "PM2": acmg.PM2_THRESHOLD},
            "pvs1_exonic_null_consequences": sorted(EXONIC_NULL),
            "pvs1_canonical_splice": {"annovar_func_refgene": "splicing", "switch": "ACMG_SPLICE_PVS1=0 turns it off",
                                       "note": "added after inspecting two validation misses; post-hoc"},
            "pvs1_review_flag_unless_gene_in": sorted(LOF_MECHANISM_GENES),
            "in_silico_columns": acmg.IN_SILICO_COLUMNS,
            "damaging_prediction_values": sorted(acmg.DAMAGING_PRED_VALUES),
            "cadd_phred_damaging_at_or_above": 20,
            "pp3_revel_calibration": {"source": "ClinGen SVI, Pejaver et al. 2022 (PMID 36413997)", "switch": "ACMG_REVEL_PP3=0 turns it off",
                                     "pathogenic_at_or_above": {s: c for c, s in eu.REVEL_PATHOGENIC},
                                     "benign_at_or_below": {s: c for c, s in eu.REVEL_BENIGN},
                                     "note": "used first; if REVEL is missing the >=3-predictor vote below is used"},
            "pp3_needs_damaging_at_least": 3,
            "bp4_needs_benign_at_least": 3,
            "class_cutoffs": {"PATHOGENIC": ">= 10", "LIKELY_PATHOGENIC": "6 to 9", "VUS": "0 to 5",
                              "LIKELY_BENIGN": "-6 to -1", "BENIGN": "<= -7, or BA1"},
            "population_frequency": "highest of gnomAD exome, gnomAD genome and the file AF columns; gnomAD 'Variant not found' with no file AF counts as absent (AF 0)",
            "clinvar": "PP5/BP6 from EXACT-allele ClinVar P/LP or B/LB (taken from MyVariant, which matches on chrom/pos/ref/alt); the NCBI position search is informational only; hidden from BOTH engines with ACMG_IGNORE_CLINVAR=1; benign ClinVar blocks PVS1",
        },
        "amp_rules.json": {
            "framework": "AMP/ASCO/CAP-style somatic tiering, evaluated independently of ACMG",
            "tier1_genes": sorted(amp.TIER1_GENES),
            "tier2_genes": sorted(amp.TIER2_GENES),
            "common_variant_allele_frequency": 0.01,
            "clinvar_and_frequency_rules": "same helpers as ACMG (exact-allele ClinVar, highest-AF across gnomAD exome/genome/file); ACMG_IGNORE_CLINVAR=1 also hides ClinVar from AMP",
            "tcga_recurrence_min_cases_default": 10,
            "tcga_env_override": "TCGA_TIER_MIN_CASES",
            "strong_evidence_means": [
                "ClinVar pathogenic (not conflicting)",
                "local annotation pathogenic / likely pathogenic",
                "null variant (frameshift, stopgain, stoploss, splicing, startloss)",
                "in-silico pathogenic (REVEL >= 0.70)",
                "TCGA recurrence >= threshold",
            ],
            "order_of_decision": [
                "common population variant -> TIER_IV",
                "ClinVar benign (unanimous / non-conflicting) -> TIER_IV",
                "annotated benign / likely benign -> TIER_IV",
                "Tier I gene + strong evidence -> TIER_I",
                "Tier II gene + strong evidence or recurrence -> TIER_II",
                "TP53 or PIK3CA hotspot/recurrence with weak evidence -> TIER_II",
                "panel gene without strong evidence -> TIER_III",
                "strong evidence, gene not in panel -> TIER_III",
                "otherwise -> TIER_III (VUS)",
            ],
        },
    }


def main() -> int:
    out_dir = ROOT / "rules"
    out_dir.mkdir(exist_ok=True)
    bad = False
    for name, content in current_rules().items():
        text = json.dumps(content, indent=2, sort_keys=False) + "\n"
        p = out_dir / name
        if "--check" in sys.argv:
            if not p.exists() or p.read_text() != text:
                print(f"OUT OF DATE: {p}")
                bad = True
        else:
            p.write_text(text)
            print(f"wrote {p}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
