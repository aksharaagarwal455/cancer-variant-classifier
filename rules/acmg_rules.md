# ACMG rules implemented

Source of truth: `backend/app/services/classification/acmg.py` and `consequence.py`.
The machine-readable copy is `acmg_rules.json` (regenerate with `python scripts/export_rules.py`; a test fails if it drifts).

This is a **documented subset**: 8 of about 28 ACMG/AMP criteria. Missense variants usually need criteria that are not
implemented (PS1, PS3, PM1, PM5 ...), so many ClinVar-pathogenic missense variants stay VUS. That is a known limitation,
not a hidden one.

## Criteria and points (ClinGen / Tavtigian 2018)

| Criterion | Strength | Points | Rule in this tool |
|---|---|---|---|
| PVS1 | Very strong, pathogenic | +8 | Frameshift, stop-gain, stop-loss, or canonical splice site (see below). Not applied if ClinVar says benign. |
| PM2 | Moderate, pathogenic | +2 | Population frequency below 0.0001 |
| PP3 | Supporting / moderate / strong, pathogenic | +1 / +2 / +4 | REVEL, ClinGen-calibrated (Pejaver 2022): >= 0.644 / 0.773 / 0.932. If REVEL is missing: at least 3 damaging predictors = supporting (CADD phred >= 20 counts as one). Switch: `ACMG_REVEL_PP3=0` |
| PP5 | Supporting, pathogenic | +1 | EXACT-allele ClinVar pathogenic / likely pathogenic (not conflicting), taken from MyVariant. Hidden from both engines when `ACMG_IGNORE_CLINVAR=1` |
| BA1 | Stand-alone, benign | forces Benign | Frequency >= 0.05 |
| BS1 | Strong, benign | -4 | Frequency >= 0.01 |
| BP4 | Supporting to very strong, benign | -1 to -8 | REVEL <= 0.290 / 0.183 / 0.016 / 0.003, else at least 3 benign predictors = supporting |
| BP6 | Supporting, benign | -1 | ClinVar benign / likely benign |

## Class cut-offs (sum of points)

| Total points | Class |
|---|---|
| >= 10 | Pathogenic |
| 6 to 9 | Likely pathogenic |
| 0 to 5 | VUS |
| -6 to -1 | Likely benign |
| <= -7 (or BA1) | Benign |

## Population frequency source

The highest of: gnomAD exome AF, gnomAD genome AF, and `AF`, `AF_popmax`, `AF_raw`, `ExAC_ALL` from the input file (the
conservative choice; the old code looked at gnomAD genome only and ignored exome-only variants). If gnomAD answers
"Variant not found" and the file has no frequency, the variant counts as absent (AF 0), which is what PM2 asks for. A gnomAD
"Invalid variant ID" error is not treated as absence.

## Splice-site PVS1 (post-hoc rule)

ANNOVAR `Func.refGene == "splicing"` means within 2 bp of an exon/intron boundary, i.e. canonical splice sites.
Such variants count as null variants for PVS1. `exonic;splicing` variants are **not** treated as null (they are
exonic missense/synonymous changes near a boundary).

This rule was added after a validation run showed that the only two missed ClinVar-pathogenic variants in the small
labelled set were canonical splice-site variants. It is clinically reasonable, but it was found by looking at those
misses, so the improvement on that set is **not independent evidence**. Switch it off with `ACMG_SPLICE_PVS1=0`.
It does not model nonsense-mediated decay escape or last-exon exceptions.
