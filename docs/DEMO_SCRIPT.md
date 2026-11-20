# Demo script (about 10 minutes)

Run `python scripts/preflight.py` and the first two commands below **before** the meeting.

1. Open with one line: "A research tool that classifies cancer variants with ACMG and AMP standards using four evidence
   sources including TCGA. Validation numbers are proof-of-concept, not a clinical benchmark."
2. Pick a file in the CLI and run from the cache (about a second):
   ```bash
   PYTHONPATH=backend python -m app.cli
   ```
   Choose `S11_labeled_subset.tsv` (number 1), limit 0, cache only = y. Point at the progress bars, the timing table (seconds per
   source and per variant) and the tier table with the first 15 variants.
3. Statistics with ClinVar hidden, rule by rule:
   ```bash
   python scripts/ablation.py
   ```
   Say it plainly: the splice rule was found by looking at the two S11 misses (not independent), and REVEL's gain is optimistic
   because REVEL's training data overlaps ClinVar. hg19 (297 variants): sensitivity 0.423 -> 0.720, specificity 0.969 -> 0.899.
4. Show `data/output/statistics_report.xlsx` (Summary, Ablation, Misses_by_consequence) and the Word report.
5. Show the limitation honestly: 46 of 75 ClinVar-pathogenic missense variants are still VUS; only 8 of about 28 ACMG criteria
   exist; AMP tiers are candidates (no OncoKB/CIViC therapy evidence).
6. Show the design choice: ClinVar is read for the exact allele (via MyVariant). The old position search gave pathogenic evidence to
   55 variants the file labels VUS or empty. Open `docs/LOGIC.md`.
7. GUI (optional): `python -m uvicorn app.api.main:app --app-dir backend --port 8000`, choose `sample_S11.tsv`, limit 30,
   tick *Cache only*. Hover the flag symbol for review flags.
8. Rules: `rules/acmg_rules.md`, `rules/amp_rules.md`; tests: `python -m pytest`.

Do not run uncached files live in a meeting (`sample_S13.tsv`, all of `hg19Somatic_validation.tsv`) unless preflight shows few live
lookups; use cache only instead.
