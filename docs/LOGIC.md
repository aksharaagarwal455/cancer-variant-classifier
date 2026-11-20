# How the tool decides (full logic)

## 1. Pipeline

```
TSV -> parse (duplicate headers, "." = null, merge duplicate chrom:pos:ref:alt)
    -> detect genome build (landmark-gene coordinates) -> normalise alleles -> genomic HGVS id
    -> evidence (4 sources in parallel; cache first: Redis -> SQLite -> live API)
         MyVariant.info  exact-allele ClinVar, CADD, dbNSFP        gnomAD  allele frequency (exome + genome)
         TCGA / GDC      tumour recurrence of the protein change   NCBI ClinVar by position (information only)
    -> ACMG engine (points)   AMP engine (tiers)   -> quality flags -> TSV / JSON / summary, CLI table or GUI
```

## 2. Evidence rules (what each source contributes)

* **ClinVar** is read for the EXACT allele (MyVariant matches chrom, position, ref, alt). The NCBI search is by position and can
  return a different allele, so it is never scored; it only raises a review flag.
* **Population frequency** = the highest of gnomAD exome, gnomAD genome and the file's AF columns. "Variant not found" in gnomAD with
  no file frequency counts as absent.
* **In-silico**: REVEL first (ClinGen-calibrated strengths), otherwise at least 3 agreeing predictors (SIFT, PolyPhen x2,
  MutationTaster, FATHMM, MetaSVM, MetaLR; CADD phred >= 20 counts as one).
* **Consequence**: frameshift, stop-gain, stop-loss, and canonical splice sites count as null variants.
* **TCGA**: number of tumour cases with the same protein change (substitutions only).

## 3. ACMG classification (points, ClinGen / Tavtigian 2018)

| Criterion | Points | Trigger |
|---|---|---|
| PVS1 | +8 | null variant (not applied if ClinVar says benign) |
| PM2 | +2 | frequency < 0.0001 or absent |
| PP3 | +1 / +2 / +4 | REVEL >= 0.644 / 0.773 / 0.932 (or 3+ damaging predictors = +1) |
| PP5 | +1 | exact-allele ClinVar pathogenic / likely pathogenic |
| BA1 | Benign | frequency >= 0.05 |
| BS1 | -4 | frequency >= 0.01 |
| BP4 | -1 to -8 | REVEL <= 0.290 / 0.183 / 0.016 / 0.003 (or 3+ benign predictors = -1) |
| BP6 | -1 | exact-allele ClinVar benign / likely benign |

Total >= 10 Pathogenic; 6-9 Likely pathogenic; 0-5 VUS; -6 to -1 Likely benign; <= -7 or BA1 Benign.
Example: frameshift (8) + rare (2) = 10 = Pathogenic. Missense with REVEL 0.95 (4) + rare (2) = 6 = Likely pathogenic.

## 4. AMP tier (evaluated separately from ACMG)

Order of decision: common variant (>= 1%) -> Tier IV; ClinVar benign -> Tier IV; Tier I gene + strong evidence -> Tier I;
Tier II gene + strong evidence -> Tier II; everything else -> Tier III. "Strong evidence" = exact-allele ClinVar pathogenic, or a
null variant, or TCGA recurrence >= 10 cases. Tier I genes: EGFR, BRAF, KRAS, NRAS, PIK3CA, ERBB2, ALK, ROS1, MET, RET, IDH1/2,
BRCA1/2, KIT, PDGFRA, FLT3, NPM1, JAK2, ABL1. Tier II genes: TP53, PTEN, APC, MYC, RB1, VHL, STK11, SMAD4, CDKN2A, ATM, PALB2,
CHEK1, BARD1, POLE.

These are **candidate** tiers. Real AMP Tier I needs therapy or guideline evidence for that variant in that tumour type, which
needs a drug-evidence source (OncoKB, CIViC) and the tumour type; neither is in this version.

## 5. Review columns written for every variant

* `amp_tier`, `amp_tier_label`, `acmg_classification`, `acmg_score`, `acmg_criteria`, explanations.
* `evidence_completeness` HIGH / MEDIUM / LOW = how many of gnomAD, MyVariant, exact-allele ClinVar and in-silico evidence actually
  reached the engines. It is NOT a probability of being right.
* `review_flags`: what was missing or suspicious (gnomAD skipped, position-only ClinVar entry ignored, PVS1 used in a gene where
  loss of function is not verified, ACMG and AMP disagree).
* `population_af`, `population_af_source`, `clinvar_exact_allele`.

## 6. Time limits

Every HTTP call has a hard timeout (5 s connect, 8-20 s read), at most 2 quick retries, no wait over 4 s, failures are never cached,
and each source stops after `--max-seconds` or repeated failures. `--offline` uses the cache only.

## 7. Known deviations from current ClinGen guidance (be ready to explain these)

| Topic | What this tool does | Current guidance / consequence |
|---|---|---|
| PM2 strength | Moderate (+2), as in the original 2015 ACMG/AMP rules and the Tavtigian points | ClinGen SVI (2020) recommends Supporting (+1). Our calls are slightly more generous for rare variants. |
| PP5 / BP6 | +1 / -1 from exact-allele ClinVar | ClinGen SVI (2018) recommends retiring them (they can double count other evidence). Kept because ClinVar is the main label source; hidden in validation. |
| PVS1 | Applied to any null variant, with a review flag outside a short list of loss-of-function genes | The full ClinGen PVS1 decision tree (last exon, NMD escape, alternative transcripts, gene mechanism) is not implemented. |
| Germline rules on tumour data | ACMG points are germline-style; AMP tiers are the somatic view | A somatic variant should be judged by AMP/ASCO/CAP; ACMG here is a pathogenicity score, not a germline report. |
| Tier I | Gene panel + strong evidence | Real Tier I needs therapy evidence for the tumour type (OncoKB / CIViC), which this version does not query. |
| REVEL | ClinGen-calibrated thresholds | REVEL training data overlaps ClinVar, so accuracy measured against ClinVar labels is optimistic. |
