# Validation

Reproduce: `python scripts/ablation.py` (cache-only, seconds) and `python scripts/make_report.py`.

## Method

* Ground truth: ClinVar significance in the input file (`CLNSIG`). Positive = Pathogenic / Likely pathogenic; conflicting and blank excluded.
* ClinVar is hidden from BOTH engines (`ACMG_IGNORE_CLINVAR=1`), so the label is not used to make the prediction.
* Prediction = ACMG class Pathogenic or Likely pathogenic.
* The hg19 file has 2251 rows, 473 unique variants, 297 with a usable label (168 positive).

## Rule ablation (ClinVar hidden, cache-only)

| Dataset | Rules | TP | FN | FP | TN | Sensitivity | Specificity | Accuracy | AUC |
|---|---|---|---|---|---|---|---|---|---|
| S11 (42) | original | 4 | 2 | 0 | 36 | 0.667 | 1.000 | 0.952 | 0.824 |
| S11 (42) | + splice PVS1 | 6 | 0 | 0 | 36 | 1.000 | 1.000 | 1.000 | 1.000 |
| S11 (42) | + calibrated REVEL only | 4 | 2 | 0 | 36 | 0.667 | 1.000 | 0.952 | 0.907 |
| S11 (42) | both (default) | 6 | 0 | 0 | 36 | 1.000 | 1.000 | 1.000 | 1.000 |
| hg19 (297) | original | 71 | 97 | 4 | 125 | 0.423 | 0.969 | 0.660 | 0.750 |
| hg19 (297) | + splice PVS1 | 92 | 76 | 5 | 124 | 0.548 | 0.961 | 0.727 | 0.827 |
| hg19 (297) | + calibrated REVEL only | 100 | 68 | 12 | 117 | 0.595 | 0.907 | 0.731 | 0.859 |
| hg19 (297) | both (default) | 121 | 47 | 13 | 116 | 0.720 | 0.899 | 0.798 | 0.920 |

"Original" also benefits from the gnomAD exome-frequency fix, which cannot be switched off (hg19 AUC 0.726 -> 0.750).

## What the numbers do and do not say

* Detected: every null-type pathogenic variant in the hg19 set (frameshift 27/27, stop-gain 44/44, splice 21/21).
* Missense: 29 of 75 ClinVar-pathogenic missense variants are now called P/LP (was 0 of 75). 46 are still VUS, and 8 ClinVar-VUS or
  benign missense variants are over-called. The single in-frame deletion is missed.
* The 5 other false positives are null-type variants that ClinVar lists as VUS or not provided; PVS1 over-calls these.

## Discussion & Considerations

1. **Splice rule sensitivity:** Inclusion of canonical splice site variants within null consequences accounts for significant sensitivity recovery on loss-of-function genes.
2. **REVEL Calibration:** The ClinGen SVI calibrated thresholds (Pejaver et al. 2022) provide standardized odds of pathogenicity for in-silico scoring without overfitting.
3. **Independent Evaluation:** Hiding ClinVar (`ACMG_IGNORE_CLINVAR=1`) demonstrates predictive performance driven strictly by molecular consequence, population frequencies, and calibrated predictors.
4. **Implementation Scope:** 8 primary ACMG criteria are currently automated. Somatic tiers incorporate curated actionability panels, recurrence, and variant consequences.
