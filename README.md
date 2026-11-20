# Cancer Variant Classifier

An automated variant interpretation platform for cancer genomics, implementing the **AMP/ASCO/CAP (2017)** somatic tiering guidelines and **ACMG/AMP (2015)** germline pathogenicity scoring.

The platform takes ANNOVAR-annotated somatic variant files (TSV), integrates biological evidence from external databases (ClinVar, gnomAD, and TCGA recurrence), computes pathogenicity and clinical actionability tiers, and provides results through a Web GUI and an interactive CLI.

---

## Features

- **AMP Somatic Tiering (Tiers I - IV)**: Evaluates variants according to clinical actionability using curated biomarker panels (e.g., TP53, PTEN, ATM, BRCA1/2, EGFR, KRAS, PIK3CA), mutation consequence, in-silico prediction scores (REVEL), and TCGA tumor cohort recurrence.
- **ACMG Pathogenicity Scoring**: Computes Bayesian points based on standard ACMG criteria (PVS1, PM2, PP3, PP5, BA1, BS1, BP4, BP6).
- **Evidence Sources**: Integrates allele frequencies from gnomAD, exact-allele clinical classifications from ClinVar (via MyVariant.info), and somatic case recurrence counts from TCGA/GDC.
- **Local SQLite Caching**: External API responses are cached in a local SQLite database (`data/variant_classifier.db`), enabling fast, reproducible runs without repetitive network latency.
- **Web Interface**: Real-time text search, multi-column filters (AMP Tier, ACMG Class, Quality, Evidence), expandable rationale cards, and 1-click export to TSV and Excel (`.xlsx`).
- **Interactive CLI**: Batch processing of TSV cohorts with progress tracking and runtime summaries.

---

## Installation

### Requirements
- Python 3.11 or 3.12

### Setup
```bash
git clone https://github.com/aksharaagarwal455/cancer-variant-classifier.git
cd cancer-variant-classifier

python -m venv venv
# On Windows:
venv\Scripts\activate
# On Linux / macOS:
source venv/bin/activate

pip install -r requirements.txt
```

---

## Usage

### 1. Web Graphical User Interface
Start the local server:
```bash
python -m uvicorn app.api.main:app --app-dir backend --port 8000 --reload
```
Open your browser at:
```
http://127.0.0.1:8000
```
On Windows, you can alternatively double-click `run_gui.bat`.

### 2. Interactive CLI
Launch the interactive command-line interface:
```bash
# On Windows:
run_cli.bat

# Or using Python directly:
python -m app.cli
```

To run non-interactively on a specific file:
```bash
python -m app.cli --input data/input/sample_S11.tsv --limit 50 --offline
```

### 3. Export Options
From the Web GUI, classified results can be exported as:
- **TSV**: Clean tab-separated table of classified variants.
- **Excel (.xlsx)**: Spreadsheet including clinical rationales, ACMG scores, and evidence breakdown.

---

## Tests

Execute the unit and integration test suite:
```bash
python -m pytest backend/tests/
```

---

## Repository Structure

```
variant-classifier/
├── backend/
│   ├── app/
│   │   ├── api/             # FastAPI web endpoints and export services
│   │   ├── database/        # SQLite database connection and cache tables
│   │   ├── models/          # Variant data models
│   │   ├── services/
│   │   │   ├── classification/ # ACMG scoring, AMP tiering, quality checks
│   │   │   ├── evidence/       # Evidence aggregation service
│   │   │   ├── normalization/  # Genomic coordinate mapping
│   │   │   ├── parser/         # TSV parser for ANNOVAR outputs
│   │   │   └── sources/        # API clients (ClinVar, gnomAD, MyVariant, TCGA)
│   │   └── static/          # Web GUI front-end assets
│   └── tests/               # Test suite
├── data/
│   ├── input/               # Sample input TSV files
│   ├── output/              # Classification results and evaluation files
│   └── variant_classifier.db # Cached evidence database
├── docs/                    # Technical documentation and evaluation reports
├── rules/                   # Rule definitions (JSON and Markdown)
├── scripts/                 # Benchmarking and evaluation utilities
├── Dockerfile               # Container deployment configuration
├── render.yaml              # Cloud deployment specification
├── requirements.txt         # Dependencies
└── run_gui.bat              # GUI launcher for Windows
```

---

## References

1. Richards, S. et al. (2015). Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. *Genetics in Medicine*, 17(5), 405-424.
2. Li, M. M. et al. (2017). Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. *The Journal of Molecular Diagnostics*, 19(1), 4-23.
3. Pejaver, V. et al. (2022). Calibration of computational tools for genomic variant pathogenicity classification using ClinGen guidelines. *Genome Medicine*, 14(1), 143.
