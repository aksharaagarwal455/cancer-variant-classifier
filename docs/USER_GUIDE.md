# User guide

All commands are run from the project root, with the virtual environment active.

## 1. One-time setup (Python 3.12)

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r backend/requirements.txt
cp .env.example .env
```
Optional: put `NCBI_API_KEY=...` in `.env` (free from NCBI; only used by the optional NCBI position search). Redis is optional (`brew install redis`, `brew services start redis`).

## 2. Check before you run anything long

```bash
python scripts/preflight.py --input data/input/hg19Somatic_validation.tsv
```
Shows what is already cached, which APIs answer (with latency) and a rough time estimate. `--no-network` skips the pings.

## 3. Measure the time per variant

```bash
python scripts/benchmark.py --input data/input/sample_S12.tsv --n 10
```
Live calls, one at a time, cache bypassed: min / median / max seconds per source. Nothing is written to the cache.

## 4. Classify

```bash
PYTHONPATH=backend python -m app.cli --input data/input/sample_S11.tsv --limit 30
PYTHONPATH=backend python -m app.cli --input data/input/sample_S11.tsv --limit 0 --offline
```
Options: `--limit N` (0 = all), `--show N` (print the first N variants with tiers; -1 = all), `--offline` (cache only), `--max-seconds N` (limit for live lookups), `--genome-build hg19|hg38`,
`--output PATH`, `--db-path PATH`. Without `--input` a file picker appears (numbered list with variant counts, then the variant limit and a cache-only question; if the genome build cannot be detected it asks you for hg19 or hg38). Progress bars show each source; a timing table
and "seconds per variant" are printed at the end.

Fast-run switches: `GNOMAD_ENABLED=0`, `TCGA_ENABLED=0` (e.g. `GNOMAD_ENABLED=0 PYTHONPATH=backend python -m app.cli ...`).

## 4b. Reading the output

Every variant gets `acmg_classification`, `amp_tier` (+ `amp_tier_label`), `evidence_completeness` (HIGH/MEDIUM/LOW) and `review_flags`.
`evidence_completeness` says how much evidence reached the engines; it is not a probability of being right. Read `docs/LOGIC.md`.

## 5. GUI

```bash
python -m uvicorn app.api.main:app --app-dir backend --port 8000
```
Open http://127.0.0.1:8000, choose a TSV, set the limit (30 is quick), pick the genome build (or leave Auto-detect), tick *Cache only* for an instant run, click Classify. A yellow warning appears if the build could not be detected. Hover the flag symbol in the last column to see the review flags. Stop with Ctrl+C.

## 6. Reproduce the validation statistics (and the Word + Excel report)

```bash
ACMG_IGNORE_CLINVAR=1 PYTHONPATH=backend python -m app.cli --input data/input/S11_labeled_subset.tsv --limit 0 --output data/output/S11_after.tsv
ACMG_SPLICE_PVS1=0 ACMG_REVEL_PP3=0 ACMG_IGNORE_CLINVAR=1 PYTHONPATH=backend python -m app.cli --input data/input/S11_labeled_subset.tsv --limit 0 --output data/output/S11_before.tsv
ACMG_IGNORE_CLINVAR=1 PYTHONPATH=backend python -m app.cli --input data/input/hg19Somatic_validation.tsv --limit 0 --output data/output/hg19_after.tsv
ACMG_SPLICE_PVS1=0 ACMG_REVEL_PP3=0 ACMG_IGNORE_CLINVAR=1 PYTHONPATH=backend python -m app.cli --input data/input/hg19Somatic_validation.tsv --limit 0 --output data/output/hg19_before.tsv
python scripts/roc_stats.py data/output/S11_after.tsv
python scripts/ablation.py
python scripts/make_report.py
```
Add `--offline` to the CLI commands for a cache-only run. `make_report.py` writes `data/output/statistics_report.xlsx`
(formula-driven), `statistics_report.docx` and `report_figures/`.
Menu alternative: `python3 choose_and_run.py`.

## 7. Tests

```bash
python -m pytest
```
50 tests cover the ACMG/AMP rules, exact-allele ClinVar, exome frequency, calibrated REVEL, review flags, the timeout/retry/circuit-breaker
logic, caching of failures, parallelism, offline mode, MyVariant batching, the API, the parser, the rules folder staying in sync with the
code, and the confidence-interval formula.

## 8. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| Looks frozen | Check the progress bars; run `preflight.py`; use `--offline` or `GNOMAD_ENABLED=0 TCGA_ENABLED=0` |
| Many `skipped` lookups | A source was failing or hit `--max-seconds`; results use input-file columns for those variants |
| gnomAD shows HTTP 429 | Rate limit; the breaker skips it after repeated 429s. Run again later (nothing bad was cached) |
| Want the NCBI position search too | `CLINVAR_NCBI_ENABLED=1` (information only; scoring uses exact-allele ClinVar from MyVariant) |
| `ModuleNotFoundError` when starting the GUI | You ran a different Python's `uvicorn`. Always use `python -m uvicorn ...` with the virtual environment active |
| `python` is not 3.12 | `/opt/homebrew/bin/python3.12 -m venv .venv` |
| ROC/confusion PNG missing | `roc_stats.py` and `make_report.py` need a labelled output file (one with a CLNSIG column) |
