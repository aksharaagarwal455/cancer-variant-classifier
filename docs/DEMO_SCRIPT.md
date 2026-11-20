# Demonstration Guide

A walkthrough for presenting the cancer variant classification platform.

### Step 1: Pre-run Verification
Ensure the environment is ready:
```bash
python scripts/preflight.py
```

### Step 2: Interactive CLI Demonstration
Run the CLI in cached mode:
```bash
python -m app.cli
```
- Select sample cohort: `data/input/sample_S11.tsv`
- Set variant limit (e.g., 30)
- View output summary table displaying AMP tiers and ACMG classifications.

### Step 3: Web GUI Interface
Launch the browser interface:
```bash
python -m uvicorn app.api.main:app --app-dir backend --port 8000 --reload
```
Open `http://localhost:8000`:
- Upload a variant TSV file.
- Demonstrate real-time search across genes and variant IDs.
- Demonstrate filtering by AMP Tier (Tier I, II, III, IV) and ACMG pathogenicity.
- Click a variant row to display the expandable clinical rationale card.
- Click **Download Excel (.xlsx)** and **Download TSV** to demonstrate instant exporting.

### Step 4: Unit Testing
Run the test suite to verify classification rule integrity:
```bash
python -m pytest backend/tests/
```
All 59 unit tests validate deterministic rule execution and database caching.
