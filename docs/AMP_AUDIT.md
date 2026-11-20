# AMP Tier Classification & Data Audit

## 1. Executive Summary & Inventory Walk
A comprehensive audit of the entire `variant-classifier/` repository was performed, inspecting every file under `data/input/`, `data/output/`, `scripts/`, `rules/`, `docs/`, `backend/app/`, and `backend/tests/`.

### Directory & File Overview
- **`data/input/`**: 5 TSV files containing cancer genomic variants (ANNOVAR / pipeline outputs).
  - `hg19Somatic_validation.tsv`: Primary somatic validation dataset.
  - `S11_labeled_subset.tsv`: Hand-curated subset of sample S11 used for ACMG ablation and test suites.
  - `sample_S11.tsv`, `sample_S12.tsv`, `sample_S13.tsv`: Somatic tumor sample callsets.
- **`data/output/`**: Benchmark runs, ablation logs (`ablation.csv`), confusion matrices, ROC plots, before/after run dumps (`hg19_before.tsv`, `hg19_after.tsv`, `S11_before.tsv`, `S11_after.tsv`), sample results (`sample_S11_results.tsv`), and generated reports (`statistics_report.docx`, `statistics_report.xlsx`).
- **`scripts/`**: Automation and validation scripts (`validate.py`, `ablation.py`, `roc_stats.py`, `make_report.py`, `preflight.py`, `benchmark.py`, `export_rules.py`).
- **`rules/`**: Rule definitions in Markdown and JSON (`acmg_rules.json`, `acmg_rules.md`, `amp_rules.json`, `amp_rules.md`). Kept strictly in sync with the codebase via `scripts/export_rules.py` and unit tests.
- **`docs/`**: Architecture, logic, performance, user guide, and ACMG validation documentation.
- **`backend/app/`**: Core FastAPI service, CLI, SQLite cache, parsers, external source fetchers (MyVariant, gnomAD, ClinVar, TCGA/GDC), coordinate/HGVS normalization, and classification engines (`acmg.py`, `amp.py`, `consequence.py`, `evidence_utils.py`, `quality.py`).
- **`backend/tests/`**: Pytest suite (50 tests) validating normalization, source retrieval, ACMG scoring, AMP tiering, API contracts, offline execution, and regression edge cases.

---

## 2. Input Dataset Audit & Comparison

Each file in `data/input/` was evaluated for row count, unique genomic variants (`CHROM:POS:REF:ALT`), total column count, presence of the `AMP` column, and column overlap:

| File Name | Total Rows | Unique Variants (`CHROM:POS:REF:ALT`) | Total Columns | Has `AMP` Column? | Notes |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **`hg19Somatic_validation.tsv`** | **2,251** | **473** | **151** | **YES** | Multi-sample cohort. Only file with ground truth `AMP` tier labels. |
| **`S11_labeled_subset.tsv`** | 42 | 42 | 149 | **NO** | 42 unique variants; used for ACMG testing with `CLNSIG` labels. |
| **`sample_S11.tsv`** | 657 | 657 | 149 | **NO** | Patient sample tumor callset. |
| **`sample_S12.tsv`** | 611 | 611 | 149 | **NO** | Patient sample tumor callset. |
| **`sample_S13.tsv`** | 647 | 647 | 149 | **NO** | Patient sample tumor callset. |

### Column Differences
- `hg19Somatic_validation.tsv` contains 151 columns. It uniquely has:
  - `AMP`: Somatic tier label (`Tier1`, `Tier2`, `Tier3`, `Tier4`).
  - `VAF`: Variant allele frequency in tumor sample.
  - `cosmic91`: COSMIC v91 annotations.
- `sample_S11.tsv`, `S12`, `S13`, and `S11_labeled_subset.tsv` contain 149 columns (they have `cosmic103` instead of `cosmic91`, and lack `AMP` and `VAF`).

---

## 3. Label Quality in `hg19Somatic_validation.tsv`

### Multi-Row Duplication
The file has **2,251 rows**, but represents only **473 unique genomic variants** (`CHROM:POS:REF:ALT`). The repeated rows represent identical variants observed across multiple samples or runs with sample-specific depths (`DP`), allele depths (`GEN[0].AD`), variant allele fractions (`VAF`), and strand bias metrics (`FS`, `SOR`).

### Conflicting AMP Labels
Across the 473 unique variants, exactly **5 variants** have conflicting `AMP` labels across rows:

| Genomic Variant (`hg19`) | Gene | Exonic Consequence | Total Rows | Label Distribution | Majority Label |
| :--- | :--- | :--- | :---: | :--- | :---: |
| `chr11:108170512:G:A` | `ATM` | nonsynonymous_SNV | 45 | Tier3: 44 (97.8%), Tier1: 1 (2.2%) | **Tier3** |
| `chr17:37687433:G:A` | `CDK12` | nonsynonymous_SNV | 291 | Tier3: 288 (99.0%), Tier1: 3 (1.0%) | **Tier3** |
| `chr17:41245006:C:A` | `BRCA1` | stopgain | 75 | Tier1: 73 (97.3%), Tier3: 2 (2.7%) | **Tier1** |
| `chr17:7577100:T:C` | `TP53` | nonsynonymous_SNV | 3 | Tier1: 2 (66.7%), Tier3: 1 (33.3%) | **Tier1** |
| `chr22:29091840:TG:CA` | `CHEK2` | nonframeshift_substitution | 101 | Tier3: 99 (98.0%), Tier1: 2 (2.0%) | **Tier3** |

### Proposed Resolution Rule
- **Deterministic Majority Vote (Mode)**: For each unique `CHROM:POS:REF:ALT`, assign the most frequent label. If a tie ever occurs, choose the lower numerical tier (or drop). In this dataset, every single conflicting variant has a clear majority (>66.7% to 99.0%).
- **Resulting Deduplicated Ground Truth (473 unique variants)**:
  - **Tier1**: 250 variants (52.85%)
  - **Tier2**: 18 variants (3.81%)
  - **Tier3**: 204 variants (43.13%)
  - **Tier4**: 1 variant (0.21%)
  - *Total*: 473 variants.
- *(Note: If conflicting variants were dropped instead of majority-voted, 468 variants remain: Tier1: 248, Tier2: 18, Tier3: 201, Tier4: 1).*

---

## 4. Current Validation Pipeline Review

We audited all validation scripts in `scripts/`:
- `scripts/validate.py`: Checks internal consistency between predicted ACMG class and AMP tier (e.g. flagging ACMG Pathogenic + AMP Tier IV), and evaluates ACMG predictions against `CLNSIG`.
- `scripts/ablation.py`: Evaluates ACMG classification sensitivity, specificity, and accuracy against ClinVar `CLNSIG` (P/LP = 1, B/LB = 0) with ClinVar hidden (`ACMG_IGNORE_CLINVAR=1`).
- `scripts/roc_stats.py`: Plots ROC curves and computes AUC for ACMG score against `CLNSIG`.
- `scripts/make_report.py`: Generates the comprehensive Word and Excel validation reports (`statistics_report.docx`, `statistics_report.xlsx`), computing metrics solely for ACMG against ClinVar `CLNSIG`.

**Key Finding**:
None of the existing validation scripts or report generators evaluates the AMP tier classifier against the `AMP` column. The `AMP` ground-truth column in `hg19Somatic_validation.tsv` has never been utilized in the project's automated benchmarking or evaluation pipelines.

---

## 5. Source and Meaning of the `AMP` Column

### Repository Review
A thorough regex and text search across all code, tests, documentation (`docs/LOGIC.md`, `VALIDATION.md`, `ARCHITECTURE.md`, `USER_GUIDE.md`, etc.), rules, and commit artifacts revealed:
1. The AMP classification logic in `docs/LOGIC.md` and `amp.py` describes AMP/ASCO/CAP somatic tiering (Tier I: strong clinical significance / FDA-approved therapies; Tier II: potential clinical significance; Tier III: unknown clinical significance / VUS; Tier IV: benign / likely benign).
2. However, **no file or documentation explains who generated the `AMP` column in `hg19Somatic_validation.tsv`**, which institution or assay produced it, which clinical guideline or software pipeline assigned the tiers, or what criteria were used.

### Open Question for the User
> **Open Question**: What is the definitive origin of `hg19Somatic_validation.tsv` and its `AMP` column? Is it derived from an institutional molecular tumor board (e.g., MSK-IMPACT, Foundation Medicine), a commercial panel annotation (e.g. Illumina TruSight / Oncomine), or a manual clinical curation?
