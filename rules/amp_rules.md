# AMP somatic tiering rules implemented

Source of truth: `backend/app/services/classification/amp.py`. Machine-readable copy: `amp_rules.json`.
AMP is evaluated **independently** of ACMG: it uses its own inputs and is never copied from the ACMG class.

## Decision order

1. Common population variant (frequency >= 0.01) -> **Tier IV**
2. ClinVar reports benign (unanimous / non-conflicting) -> **Tier IV**
3. Annotated benign / likely benign (`variant_classification`) -> **Tier IV**
4. Gene in the Tier I panel **and** strong evidence -> **Tier I**
5. Gene in the Tier II panel **and** strong evidence or recurrence -> **Tier II**
6. TP53 or PIK3CA hotspot/recurrence with weak evidence -> **Tier II**
7. Panel gene without strong evidence -> **Tier III**
8. Strong evidence but gene not in panel -> **Tier III**
9. Otherwise -> **Tier III** (unknown significance / VUS)

"Strong evidence" means any of:
- ClinVar pathogenic (not conflicting)
- Local annotation pathogenic or likely pathogenic (`variant_classification`)
- Null variant (same definition as PVS1: frameshift, stopgain, stoploss, splicing, startloss)
- In-silico pathogenic prediction (`REVEL_score >= 0.70`)
- TCGA recurrence of at least 10 cases (`TCGA_TIER_MIN_CASES` changes the threshold)

## Gene panels

The gene lists are in `amp_rules.json` (`tier1_genes`, `tier2_genes`):
- **Tier I genes**: Curated Level A/B cancer biomarkers and core homologous recombination / DNA damage repair genes with therapeutic or established clinical significance:
  `ABL1`, `ALK`, `ATM`, `BARD1`, `BRAF`, `BRCA1`, `BRCA2`, `BRIP1`, `CHEK1`, `CHEK2`, `EGFR`, `ERBB2`, `FANCA`, `FANCL`, `FLT3`, `IDH1`, `IDH2`, `JAK2`, `KIT`, `KRAS`, `MET`, `NPM1`, `NRAS`, `PDGFRA`, `PIK3CA`, `PTEN`, `RAD51C`, `RET`, `ROS1`, `TP53`.
- **Tier II genes**: Emerging or secondary actionable cancer genes:
  `APC`, `CDK12`, `CDKN2A`, `MYC`, `PALB2`, `POLE`, `RB1`, `SMAD4`, `STK11`, `VHL`, `XRCC2`.

These are curated research panels for this platform (no live OncoKB/CIViC therapy lookups are performed).

## TCGA evidence

Counts come from the NCI GDC `ssm_occurrences` endpoint, matched by gene and protein change (not by genomic
coordinate, because GDC is GRCh38). Counts are approximate and are used only for AMP, not for ACMG scoring.

ClinVar, population frequency and the null-variant rule are read through the same helpers as ACMG, so both engines see the
same evidence. `ACMG_IGNORE_CLINVAR=1` hides ClinVar from AMP as well (earlier versions hid it from ACMG only).
