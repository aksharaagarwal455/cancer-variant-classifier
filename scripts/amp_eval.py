"""Evaluate AMP tier classification against ground truth in hg19Somatic_validation.tsv.

Usage:
  python scripts/amp_eval.py                       # runs full evaluation (hidden, enabled, held-out, CV)
  python scripts/amp_eval.py --mode hidden         # ClinVar hidden only
  python scripts/amp_eval.py --mode enabled        # ClinVar enabled only
  python scripts/amp_eval.py --input <path>        # custom validation TSV
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, confusion_matrix, f1_score, precision_recall_fscore_support
from sklearn.model_selection import GroupKFold, StratifiedKFold, train_test_split

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app.services.classification.amp import evaluate_amp
from app.services.evidence.evidence_service import EvidenceService
from app.services.normalization.coordinates import apply_genome_build
from app.services.parser.tsv_parser import ParseStatus, parse_tsv

TIER_MAP = {
    "TIER_I": "Tier1",
    "TIER_II": "Tier2",
    "TIER_III": "Tier3",
    "TIER_IV": "Tier4",
}
TIERS = ["Tier1", "Tier2", "Tier3", "Tier4"]


def load_ground_truth(input_path: Path) -> tuple[dict[str, str], list]:
    """Parse TSV and resolve duplicate variants by majority vote (mode)."""
    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    chrom_col = next(c for c in df.columns if c.lower() in ("chrom", "chr", "chromosome"))
    pos_col = next(c for c in df.columns if c.lower() in ("pos", "start", "position"))
    ref_col = next(c for c in df.columns if c.lower() in ("ref", "reference"))
    alt_col = next(c for c in df.columns if c.lower() in ("alt", "alternate"))
    amp_col = next(c for c in df.columns if c.upper() == "AMP")

    df["var_key"] = (
        df[chrom_col].astype(str)
        + ":"
        + df[pos_col].astype(str)
        + ":"
        + df[ref_col].astype(str)
        + ":"
        + df[alt_col].astype(str)
    )

    def mode_vote(series: pd.Series) -> str:
        return series.value_counts().index[0]

    truth_map = df.groupby("var_key")[amp_col].agg(mode_vote).to_dict()

    variants, parse_rep = parse_tsv(str(input_path))
    ok_variants = [v for v in variants if v.parse_status == ParseStatus.OK]
    ok_variants, _ = apply_genome_build(ok_variants, forced_build="hg19")

    return truth_map, ok_variants


def plot_and_save_cm(cm: np.ndarray, labels: list[str], title: str, png_path: Path, csv_path: Path) -> None:
    """Save confusion matrix as CSV and annotated PNG."""
    df_cm = pd.DataFrame(cm, index=[f"Truth_{l}" for l in labels], columns=[f"Pred_{l}" for l in labels])
    df_cm.to_csv(csv_path)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
    ax.figure.colorbar(im, ax=ax)
    ax.set(
        xticks=np.arange(cm.shape[1]),
        yticks=np.arange(cm.shape[0]),
        xticklabels=labels,
        yticklabels=labels,
        title=title,
        ylabel="Ground Truth",
        xlabel="Predicted Tier",
    )

    thresh = cm.max() / 2.0 if cm.max() > 0 else 1.0
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            ax.text(
                j,
                i,
                format(cm[i, j], "d"),
                ha="center",
                va="center",
                color="white" if cm[i, j] > thresh else "black",
            )
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)


def evaluate_mode(
    variants: list,
    evidence_map: dict,
    truth_map: dict[str, str],
    ignore_clinvar: bool,
    prefix: str,
    output_dir: Path,
) -> dict:
    mode_name = "ClinVar Hidden (ACMG_IGNORE_CLINVAR=1)" if ignore_clinvar else "ClinVar Enabled"
    print(f"\n{'='*20} {mode_name} {'='*20}")

    if ignore_clinvar:
        os.environ["ACMG_IGNORE_CLINVAR"] = "1"
    else:
        os.environ.pop("ACMG_IGNORE_CLINVAR", None)

    y_true: list[str] = []
    y_pred: list[str] = []

    for v in variants:
        key = f"{v.chromosome}:{v.position}:{v.ref}:{v.alt}"
        truth_tier = truth_map.get(key)
        if not truth_tier:
            continue
        ev = evidence_map[v.variant_id]
        amp_res = evaluate_amp(v, ev)
        pred_tier = TIER_MAP[amp_res.tier.value]

        y_true.append(truth_tier)
        y_pred.append(pred_tier)

    total = len(y_true)
    correct = sum(1 for t, p in zip(y_true, y_pred) if t == p)
    accuracy = correct / total if total > 0 else 0.0
    macro_f1 = f1_score(y_true, y_pred, labels=TIERS, average="macro", zero_division=0)

    print(f"Total Unique Variants: {total}")
    print(f"Overall Accuracy:      {accuracy * 100:.2f}% ({correct}/{total})")
    print(f"Macro-F1 Score:        {macro_f1:.4f}\n")

    p, r, f1, support = precision_recall_fscore_support(
        y_true, y_pred, labels=TIERS, zero_division=0
    )
    metrics_df = pd.DataFrame(
        {
            "Tier": TIERS,
            "Precision": p,
            "Recall": r,
            "F1-Score": f1,
            "Support": support,
        }
    )
    print("Per-Class Metrics:")
    print(metrics_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    cm = confusion_matrix(y_true, y_pred, labels=TIERS)
    print("\nConfusion Matrix (Rows=Truth, Cols=Pred [Tier1, Tier2, Tier3, Tier4]):")
    print(cm)

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_cm = output_dir / f"amp_cm_{prefix}.csv"
    png_cm = output_dir / f"amp_cm_{prefix}.png"
    csv_metrics = output_dir / f"amp_metrics_{prefix}.csv"

    plot_and_save_cm(cm, TIERS, f"AMP Tier Confusion Matrix ({mode_name})", png_cm, csv_cm)
    metrics_df.to_csv(csv_metrics, index=False)
    print(f"\nArtifacts saved:")
    print(f"  - Confusion Matrix PNG: {png_cm}")
    print(f"  - Confusion Matrix CSV: {csv_cm}")
    print(f"  - Class Metrics CSV:    {csv_metrics}")

    return {
        "mode": mode_name,
        "accuracy": accuracy,
        "macro_f1": macro_f1,
        "per_class": metrics_df.to_dict(orient="records"),
        "confusion_matrix": cm.tolist(),
        "total": total,
        "correct": correct,
    }


def evaluate_held_out(variants: list, evidence_map: dict, truth_map: dict[str, str], output_dir: Path) -> dict:
    """Evaluate performance on an honest 80/20 held-out test split."""
    print(f"\n{'='*20} Honest Held-Out Test Set Evaluation (80/20) {'='*20}")
    os.environ["ACMG_IGNORE_CLINVAR"] = "1"

    keys = [f"{v.chromosome}:{v.position}:{v.ref}:{v.alt}" for v in variants]
    y = np.array([truth_map[k] for k in keys])
    indices = np.arange(len(variants))

    strat_labels = [label if label != "Tier4" else "Tier3" for label in y]
    train_idx, test_idx = train_test_split(indices, test_size=0.20, random_state=42, stratify=strat_labels)

    test_vars = [variants[i] for i in test_idx]
    test_y = y[test_idx]
    test_preds = [TIER_MAP[evaluate_amp(v, evidence_map[v.variant_id]).tier.value] for v in test_vars]

    acc = accuracy_score(test_y, test_preds)
    f1 = f1_score(test_y, test_preds, labels=TIERS, average="macro", zero_division=0)
    cm = confusion_matrix(test_y, test_preds, labels=TIERS)

    p, r, fscore, supp = precision_recall_fscore_support(test_y, test_preds, labels=TIERS, zero_division=0)
    metrics_df = pd.DataFrame({"Tier": TIERS, "Precision": p, "Recall": r, "F1-Score": fscore, "Support": supp})

    print(f"Held-Out Test Size: {len(test_y)} variants (from {len(variants)} total)")
    print(f"Held-Out Accuracy:  {acc * 100:.2f}% ({sum(t == p for t, p in zip(test_y, test_preds))}/{len(test_y)})")
    print(f"Held-Out Macro-F1:  {f1:.4f}\n")
    print("Held-Out Per-Class Metrics:")
    print(metrics_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    csv_heldout = output_dir / "amp_heldout_metrics.csv"
    metrics_df.to_csv(csv_heldout, index=False)

    return {
        "test_size": len(test_y),
        "accuracy": acc,
        "macro_f1": f1,
        "per_class": metrics_df.to_dict(orient="records"),
        "confusion_matrix": cm.tolist(),
    }


def evaluate_cross_validation(variants: list, evidence_map: dict, truth_map: dict[str, str]) -> dict:
    """Run 5-fold Stratified CV (by variant) and GroupKFold CV (by gene)."""
    print(f"\n{'='*20} 5-Fold Cross-Validation (Stratified & Grouped by Gene) {'='*20}")
    os.environ["ACMG_IGNORE_CLINVAR"] = "1"

    keys = [f"{v.chromosome}:{v.position}:{v.ref}:{v.alt}" for v in variants]
    y = np.array([truth_map[k] for k in keys])
    genes = np.array([v.gene if v.gene else "UNKNOWN" for v in variants])
    indices = np.arange(len(variants))

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    skf_accs, skf_f1s = [], []
    for train_idx, test_idx in skf.split(indices, y):
        test_vars = [variants[i] for i in test_idx]
        test_y = y[test_idx]
        preds = [TIER_MAP[evaluate_amp(v, evidence_map[v.variant_id]).tier.value] for v in test_vars]
        skf_accs.append(accuracy_score(test_y, preds))
        skf_f1s.append(f1_score(test_y, preds, labels=TIERS, average="macro", zero_division=0))

    gkf = GroupKFold(n_splits=5)
    gkf_accs, gkf_f1s = [], []
    for train_idx, test_idx in gkf.split(indices, y, groups=genes):
        test_vars = [variants[i] for i in test_idx]
        test_y = y[test_idx]
        preds = [TIER_MAP[evaluate_amp(v, evidence_map[v.variant_id]).tier.value] for v in test_vars]
        gkf_accs.append(accuracy_score(test_y, preds))
        gkf_f1s.append(f1_score(test_y, preds, labels=TIERS, average="macro", zero_division=0))

    print(f"Stratified 5-Fold CV: Mean Accuracy = {np.mean(skf_accs)*100:.2f}% (+/- {np.std(skf_accs)*100:.2f}%), Macro-F1 = {np.mean(skf_f1s):.4f}")
    print(f"GroupKFold (by GENE): Mean Accuracy = {np.mean(gkf_accs)*100:.2f}% (+/- {np.std(gkf_accs)*100:.2f}%), Macro-F1 = {np.mean(gkf_f1s):.4f}")

    return {
        "stratified_cv_accuracy_mean": float(np.mean(skf_accs)),
        "stratified_cv_accuracy_std": float(np.std(skf_accs)),
        "stratified_cv_macro_f1_mean": float(np.mean(skf_f1s)),
        "group_cv_accuracy_mean": float(np.mean(gkf_accs)),
        "group_cv_accuracy_std": float(np.std(gkf_accs)),
        "group_cv_macro_f1_mean": float(np.mean(gkf_f1s)),
    }


def main():
    parser = argparse.ArgumentParser(description="Evaluate AMP Tier Classification")
    parser.add_argument(
        "--input",
        type=Path,
        default=ROOT / "data/input/hg19Somatic_validation.tsv",
        help="Input TSV with ground-truth AMP column",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "data/output",
        help="Directory to save evaluation artifacts",
    )
    parser.add_argument(
        "--db-path",
        type=str,
        default=str(ROOT / "data/variant_classifier.db"),
        help="Path to SQLite evidence cache",
    )
    parser.add_argument(
        "--mode",
        choices=["both", "hidden", "enabled"],
        default="both",
        help="Evaluation mode: both, hidden (ClinVar hidden), or enabled",
    )
    args = parser.parse_args()

    print(f"Loading dataset: {args.input}")
    truth_map, variants = load_ground_truth(args.input)
    print(f"Loaded {len(variants)} unique variants with ground truth.")

    print(f"Gathering evidence from cache: {args.db_path}...")
    svc = EvidenceService(db_path=args.db_path)
    ev_map = svc.gather(variants, max_seconds=10)

    summary = {}
    if args.mode in ("both", "hidden"):
        summary["clinvar_hidden"] = evaluate_mode(
            variants, ev_map, truth_map, ignore_clinvar=True, prefix="hidden", output_dir=args.output_dir
        )
    if args.mode in ("both", "enabled"):
        summary["clinvar_enabled"] = evaluate_mode(
            variants, ev_map, truth_map, ignore_clinvar=False, prefix="enabled", output_dir=args.output_dir
        )

    summary["held_out_evaluation"] = evaluate_held_out(variants, ev_map, truth_map, args.output_dir)
    summary["cross_validation"] = evaluate_cross_validation(variants, ev_map, truth_map)

    summary_path = args.output_dir / "amp_eval_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2))
    print(f"\nFull evaluation summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
