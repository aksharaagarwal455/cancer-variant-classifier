# Cancer Variant Classification Platform

A research tool that classifies annotated cancer variants with two frameworks, **ACMG/AMP pathogenicity points** and
**AMP/ASCO/CAP-style tiers (I-IV)**, using four evidence sources: MyVariant.info (exact-allele ClinVar, CADD, dbNSFP),
gnomAD (allele frequency), TCGA/GDC (tumour recurrence) and an NCBI ClinVar position check.
Every result is cached (Redis if present, then SQLite), every network call has a hard timeout, and every variant carries
review flags that say what evidence was missing.

> Research / evidence-aggregation tool. Not a clinical diagnostic device. Requires review by a qualified clinical genomics
> professional before any clinical use.

## Quick start (Python 3.12)

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -r backend/requirements.txt
cp .env.example .env
PYTHONPATH=backend python -m pytest            # 50 tests, no internet needed
PYTHONPATH=backend python -m app.cli           # interactive: pick a file, a limit, cache-only yes/no
python -m uvicorn app.api.main:app --app-dir backend --port 8000     # browser GUI at http://127.0.0.1:8000
```

Non-interactive: `PYTHONPATH=backend python -m app.cli --input data/input/sample_S11.tsv --limit 30`
(`--offline` = cache only, `--max-seconds N` = limit for live lookups, `--show N` = print N variants with tiers,
`--genome-build hg19|hg38` = override detection).

## What you get for every variant

`acmg_classification` and `acmg_score`, `amp_tier` with a plain-English `amp_tier_label`, `evidence_completeness`
(HIGH / MEDIUM / LOW), `review_flags`, `population_af`, `clinvar_exact_allele`, source badges and explanations. A run also writes a
summary JSON with seconds per source and per variant.

## Validation (ClinVar hidden from the engines)

hg19 set, 297 labelled variants (168 ClinVar-pathogenic): sensitivity 0.720, specificity 0.899, AUC 0.920 (original rules:
0.423 / 0.969 / 0.750). Small S11 set: 6 of 6 recovered, but that gain is post-hoc. REVEL's gain is optimistic. Proof of concept,
not a clinical benchmark. Details: `docs/VALIDATION.md`; reproduce with `python scripts/ablation.py` and
`python scripts/make_report.py` (writes `data/output/statistics_report.docx` and `.xlsx`).

## Where things are

| Path | Contents |
|---|---|
| `docs/LOGIC.md` | how every decision is made, and known deviations from ClinGen guidance |
| `docs/USER_GUIDE.md`, `DEMO_SCRIPT.md` | commands, demo order, troubleshooting |
| `docs/ARCHITECTURE.md`, `PERFORMANCE.md`, `VALIDATION.md` | design, timings, statistics |
| `rules/` | the ACMG and AMP rules as `.md` and `.json` (a test keeps them identical to the code) |
| `backend/app/` | parser, genome-build detection, sources, evidence service, ACMG/AMP engines, CLI, API + GUI |
| `scripts/` | `preflight.py`, `benchmark.py`, `ablation.py`, `roc_stats.py`, `validate.py`, `make_report.py`, `export_rules.py` |
| `data/input/`, `data/output/` | sample inputs, current results, statistics reports |

## Limits

Only 8 of about 28 ACMG criteria are implemented (missense variants often stay VUS); AMP tiers are candidates (no OncoKB / CIViC
therapy lookup); ACMG points are germline-style. Read `docs/LOGIC.md` section 7 before presenting.
