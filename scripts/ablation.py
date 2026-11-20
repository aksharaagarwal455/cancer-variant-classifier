"""Rule ablation: run each labelled set with the rule switches off/on and write data/output/ablation.csv.

  python scripts/ablation.py          # cache-only (offline), takes a few seconds

ClinVar is hidden from the engines (ACMG_IGNORE_CLINVAR=1) so the label is not used for the prediction.
Switches: ACMG_SPLICE_PVS1 (canonical splice site counts as null for PVS1),
          ACMG_REVEL_PP3  (ClinGen-calibrated REVEL strength for PP3/BP4).
"""
from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DATASETS = [("S11 labelled subset", "S11_labeled_subset"), ("hg19 somatic validation", "hg19Somatic_validation")]
CONFIGS = [("original rules", "0", "0"), ("+ splice-site PVS1", "1", "0"),
           ("+ calibrated REVEL only", "0", "1"), ("both (current default)", "1", "1")]


def truth(c: str):
    c = c.lower()
    if "conflicting" in c or c == "":
        return None
    return 1 if "pathogenic" in c else 0


def main() -> None:
    rows = []
    for title, stem in DATASETS:
        for name, splice, revel in CONFIGS:
            out = Path(tempfile.mkdtemp()) / f"{stem}.tsv"
            env = {**os.environ, "PYTHONPATH": str(ROOT / "backend"), "ACMG_IGNORE_CLINVAR": "1",
                   "ACMG_SPLICE_PVS1": splice, "ACMG_REVEL_PP3": revel, "EVIDENCE_OFFLINE": "1"}
            subprocess.run([sys.executable, "-m", "app.cli", "--input", str(ROOT / "data/input" / f"{stem}.tsv"),
                            "--limit", "0", "--offline", "--output", str(out)], cwd=ROOT, env=env,
                           check=True, capture_output=True)
            df = pd.read_csv(out, sep="\t", dtype=str).fillna("")
            df["y"] = df["CLNSIG"].apply(truth)
            df["score"] = pd.to_numeric(df["acmg_score"], errors="coerce")
            df = df[df["y"].notna() & df["score"].notna()]
            y = df["y"].astype(int).values
            pred = df["acmg_classification"].isin(["PATHOGENIC", "LIKELY_PATHOGENIC"]).astype(int).values
            tn, fp, fn, tp = confusion_matrix(y, pred, labels=[0, 1]).ravel()
            rows.append({"dataset": title, "rules": name, "n": len(y), "positives": int(y.sum()), "TP": int(tp), "FN": int(fn),
                         "FP": int(fp), "TN": int(tn), "sensitivity": round(tp / (tp + fn), 3),
                         "specificity": round(tn / (tn + fp), 3), "accuracy": round((tp + tn) / len(y), 3),
                         "AUC": round(roc_auc_score(y, df["score"]), 3)})
    out_csv = ROOT / "data/output/ablation.csv"
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    main()
