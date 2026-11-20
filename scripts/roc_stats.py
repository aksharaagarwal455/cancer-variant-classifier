"""
Statistical validation: sensitivity, specificity, confusion matrix, ROC, AUC.
Ground truth: ClinVar CLNSIG (binarized P/LP = 1, else = 0).
"""
import sys
import pandas as pd
from sklearn.metrics import (
    roc_curve, roc_auc_score, confusion_matrix, ConfusionMatrixDisplay
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "data/output/S11_after.tsv"
df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

def truth_binary(clnsig):
    c = clnsig.lower()
    if "conflicting" in c or c == "":
        return None
    if "pathogenic" in c:
        return 1
    return 0

df["y_true"] = df["CLNSIG"].apply(truth_binary)
df["acmg_score"] = pd.to_numeric(df["acmg_score"], errors="coerce")
df["y_pred_binary"] = (df["acmg_classification"].isin(["PATHOGENIC", "LIKELY_PATHOGENIC"])).astype(int)

valid = df[df["y_true"].notna() & df["acmg_score"].notna()].copy()
valid["y_true"] = valid["y_true"].astype(int)

n = len(valid)
n_pos = int(valid["y_true"].sum())
n_neg = n - n_pos

print(f"File: {path}")
print(f"Usable rows: {n}  (positives/P-LP: {n_pos}, negatives/VUS+benign: {n_neg})")
print()

y_true = valid["y_true"].values
y_pred = valid["y_pred_binary"].values
tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

sensitivity = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
specificity = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
accuracy = (tp + tn) / n

print("=== CONFUSION MATRIX ===")
print("                 Predicted P/LP   Predicted not-P/LP")
print(f"Actual P/LP           TP={tp}              FN={fn}")
print(f"Actual not-P/LP       FP={fp}              TN={tn}")
print()
print(f"Sensitivity (Recall) : {sensitivity:.3f}  ({tp}/{tp+fn})")
print(f"Specificity           : {specificity:.3f}  ({tn}/{tn+fp})")
if (tp + fp) > 0:
    print(f"Precision (PPV)       : {precision:.3f}")
else:
    print("Precision (PPV)       : undefined (0 predicted positives)")
print(f"Overall accuracy       : {accuracy:.3f}")
print()

if n_pos > 0 and n_neg > 0:
    fpr, tpr, thresholds = roc_curve(y_true, valid["acmg_score"].values)
    auc = roc_auc_score(y_true, valid["acmg_score"].values)
    print("=== ROC / AUC ===")
    print(f"AUC = {auc:.3f}")

    plt.figure(figsize=(5, 5))
    plt.plot(fpr, tpr, label=f"ACMG score (AUC = {auc:.3f})", linewidth=2)
    plt.plot([0, 1], [0, 1], "--", color="gray", label="Random guess (AUC = 0.5)")
    plt.xlabel("False Positive Rate (1 - Specificity)")
    plt.ylabel("True Positive Rate (Sensitivity)")
    plt.title("ROC Curve: ACMG Score vs ClinVar P/LP")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig("data/output/roc_curve.png", dpi=150)
    print("Saved: data/output/roc_curve.png")

    disp = ConfusionMatrixDisplay(confusion_matrix(y_true, y_pred, labels=[0, 1]),
                                   display_labels=["Not P/LP", "P/LP"])
    disp.plot(cmap="Blues")
    plt.title("Confusion Matrix")
    plt.tight_layout()
    plt.savefig("data/output/confusion_matrix.png", dpi=150)
    print("Saved: data/output/confusion_matrix.png")
else:
    print("Cannot compute ROC/AUC: need both positive and negative cases.")

print()
print("=" * 60)
print(f"CAVEAT: n={n}, only {n_pos} positive cases. This is a proof-of-concept")
print("validation, not a statistically powered clinical accuracy study.")
print("Treat as directional evidence, not a performance guarantee.")
print("=" * 60)
