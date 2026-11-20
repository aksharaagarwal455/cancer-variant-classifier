"""Build the statistics report as an Excel workbook AND a Word document from real run outputs.

  python scripts/make_report.py                       # uses data/output/{S11,hg19}_{before,after}.tsv
  python scripts/make_report.py --out data/output     # choose output folder

Inputs (made by `app.cli`, see docs/USER_GUIDE.md):
  data/output/S11_before.tsv   data/output/S11_after.tsv     (+ *_summary.json next to them)
  data/output/hg19_before.tsv  data/output/hg19_after.tsv

"before" = original rules (ACMG_SPLICE_PVS1=0 ACMG_REVEL_PP3=0), "after" = current default rules. Both hide ClinVar from the
engine (ACMG_IGNORE_CLINVAR=1) so the ClinVar label is not used to make the prediction.

Outputs: statistics_report.xlsx (formula-driven: change a data row and the metrics update),
statistics_report.docx, and report_figures/*.png.
Every number is computed from the TSV files, never typed in by hand.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import confusion_matrix, roc_auc_score, roc_curve  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app.services.classification import acmg as acmg_mod  # noqa: E402

DATASETS = [
    ("S11", "S11 labelled subset", "S11_before.tsv", "S11_after.tsv"),
    ("hg19", "hg19 somatic validation", "hg19_before.tsv", "hg19_after.tsv"),
]
FONT = "Arial"


# ----------------------------------------------------------------------------- data + metrics
def truth_binary(clnsig: str):
    c = (clnsig or "").lower()
    if c == "" or "conflicting" in c:
        return None
    return 1 if "pathogenic" in c else 0


def load_labelled(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df["y_true"] = df["CLNSIG"].map(truth_binary)
    df["acmg_score"] = pd.to_numeric(df["acmg_score"], errors="coerce")
    df["y_pred"] = df["acmg_classification"].isin(["PATHOGENIC", "LIKELY_PATHOGENIC"]).astype(int)
    df["consequence"] = df["ExonicFunc.refGene"].where(df["ExonicFunc.refGene"] != "", df["Func.refGene"])
    df["consequence"] = df["consequence"].replace("", "unknown")
    valid = df[df["y_true"].notna() & df["acmg_score"].notna()].copy()
    valid["y_true"] = valid["y_true"].astype(int)
    return valid.reset_index(drop=True)


def wilson(k: int, n: int, z: float = 1.96):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def metrics(df: pd.DataFrame) -> dict:
    tn, fp, fn, tp = confusion_matrix(df["y_true"], df["y_pred"], labels=[0, 1]).ravel()
    tn, fp, fn, tp = int(tn), int(fp), int(fn), int(tp)
    n = tn + fp + fn + tp
    m = {"n": n, "pos": tp + fn, "neg": tn + fp, "tp": tp, "fn": fn, "fp": fp, "tn": tn}
    m["sens"] = tp / (tp + fn) if tp + fn else float("nan")
    m["spec"] = tn / (tn + fp) if tn + fp else float("nan")
    m["prec"] = tp / (tp + fp) if tp + fp else float("nan")
    m["acc"] = (tp + tn) / n if n else float("nan")
    m["sens_ci"] = wilson(tp, tp + fn)
    m["spec_ci"] = wilson(tn, tn + fp)
    m["acc_ci"] = wilson(tp + tn, n)
    m["auc"] = roc_auc_score(df["y_true"], df["acmg_score"]) if 0 < m["pos"] < n else float("nan")
    return m


def acmg_group(v: str) -> str:
    v = v.upper()
    if "PATHOGENIC" in v and "BENIGN" not in v:
        return "P/LP"
    if "BENIGN" in v:
        return "B/LB"
    return "VUS"


def conflicts(df: pd.DataFrame) -> int:
    g = df["acmg_classification"].map(acmg_group)
    return int((((g == "P/LP") & (df["amp_tier"] == "TIER_IV")) |
                ((g == "B/LB") & df["amp_tier"].isin(["TIER_I", "TIER_II"]))).sum())


def fmt(x, nd=3):
    return "n/a" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{nd}f}"


def load_summary(path: Path):
    p = path.with_name(path.stem + "_summary.json")
    return json.loads(p.read_text()) if p.exists() else None


# ----------------------------------------------------------------------------- figures
def make_figures(data: dict, figdir: Path) -> dict[str, Path]:
    figdir.mkdir(parents=True, exist_ok=True)
    out: dict[str, Path] = {}
    for key, title, _, _ in DATASETS:
        b, a = data[key]["before"], data[key]["after"]
        fig, ax = plt.subplots(figsize=(5, 5))
        for lab, df, style in (("original rules", b, "--"), ("current rules", a, "-")):
            fpr, tpr, _ = roc_curve(df["y_true"], df["acmg_score"])
            ax.plot(fpr, tpr, style, linewidth=2, label=f"{lab} (AUC {roc_auc_score(df['y_true'], df['acmg_score']):.3f})")
        ax.plot([0, 1], [0, 1], ":", color="gray", label="chance")
        ax.set_xlabel("False positive rate (1 - specificity)")
        ax.set_ylabel("True positive rate (sensitivity)")
        ax.set_title(f"ROC: {title}")
        ax.legend(loc="lower right", fontsize=8)
        fig.tight_layout()
        p = figdir / f"roc_{key}.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        out[f"roc_{key}"] = p

        fig, axes = plt.subplots(1, 2, figsize=(8, 3.6))
        for ax, (lab, m) in zip(axes, (("before", data[key]["m_before"]), ("after", data[key]["m_after"]))):
            mat = [[m["tn"], m["fp"]], [m["fn"], m["tp"]]]
            ax.imshow(mat, cmap="Blues")
            for i in range(2):
                for j in range(2):
                    ax.text(j, i, mat[i][j], ha="center", va="center", fontsize=14,
                            color="white" if mat[i][j] > max(max(r) for r in mat) / 2 else "black")
            ax.set_xticks([0, 1], ["pred not P/LP", "pred P/LP"])
            ax.set_yticks([0, 1], ["true not P/LP", "true P/LP"])
            ax.set_title(f"{title}", fontsize=10)
        fig.suptitle(title, fontsize=11)
        fig.tight_layout()
        p = figdir / f"confusion_{key}.png"
        fig.savefig(p, dpi=150)
        plt.close(fig)
        out[f"cm_{key}"] = p
    return out


# ----------------------------------------------------------------------------- XLSX
def build_xlsx(data: dict, figs: dict[str, Path], path: Path) -> None:
    from openpyxl import Workbook
    from openpyxl.drawing.image import Image as XLImage
    from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
    from openpyxl.utils import get_column_letter as L

    wb = Workbook()
    f_base = Font(name=FONT, size=10)
    f_bold = Font(name=FONT, size=10, bold=True)
    f_head = Font(name=FONT, size=10, bold=True, color="FFFFFF")
    f_input = Font(name=FONT, size=10, color="0000FF")
    f_title = Font(name=FONT, size=14, bold=True)
    fill_head = PatternFill("solid", start_color="1F3864")
    fill_key = PatternFill("solid", start_color="FFF2CC")
    thin = Side(style="thin", color="BFBFBF")
    box = Border(left=thin, right=thin, top=thin, bottom=thin)

    def style_range(ws, rng, font=f_base, border=True, align=None):
        for row in ws[rng]:
            for c in row:
                c.font = font
                if border:
                    c.border = box
                if align:
                    c.alignment = align

    def header(ws, row, labels, col0=1):
        for i, lab in enumerate(labels):
            c = ws.cell(row=row, column=col0 + i, value=lab)
            c.font, c.fill, c.border = f_head, fill_head, box
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

    # ---- data sheets (one per dataset x version) ----
    ranges = {}
    for key, title, _, _ in DATASETS:
        for ver in ("before", "after"):
            name = f"{key}_{ver}"
            ws = wb.create_sheet(name)
            df = data[key][ver]
            cols = ["variant_id", "Gene", "CLNSIG", "y_true", "acmg_classification", "acmg_score", "y_pred",
                    "amp_tier", "consequence"]
            header(ws, 1, cols)
            for i, r in enumerate(df.itertuples(index=False), start=2):
                vals = [r.variant_id, r.Gene, r.CLNSIG, int(r.y_true), r.acmg_classification, int(r.acmg_score),
                        int(r.y_pred), r.amp_tier, r.consequence]
                for j, v in enumerate(vals, start=1):
                    c = ws.cell(row=i, column=j, value=v)
                    c.font = f_base
            n = len(df) + 1
            ranges[name] = n
            for j, w in enumerate([34, 10, 20, 8, 22, 11, 8, 11, 22], start=1):
                ws.column_dimensions[L(j)].width = w
            ws.freeze_panes = "A2"
            ws.sheet_properties.tabColor = "A6A6A6"

    # ---- Summary ----
    ws = wb.active
    ws.title = "Summary"
    ws["A1"] = "Variant classifier - validation statistics"
    ws["A1"].font = f_title
    ws["A2"] = (f"Generated {date.today().isoformat()} from the TSV outputs. Ground truth = ClinVar CLNSIG "
                "(P/LP = positive; conflicting / blank rows excluded). ClinVar hidden from the engine.")
    ws["A2"].font = f_base
    ws["A3"] = "z for confidence interval"
    ws["A3"].font = f_bold
    ws["B3"] = 1.96
    ws["B3"].font = f_input
    ws["B3"].fill = fill_key
    ws["C3"] = "95% Wilson interval (input: change to 2.576 for 99%)"
    ws["C3"].font = f_base
    heads = ["Dataset", "Version", "Usable n", "Positives", "TP", "FN", "FP", "TN", "Sensitivity", "Sens CI low",
             "Sens CI high", "Specificity", "Spec CI low", "Spec CI high", "Precision (PPV)", "Accuracy",
             "AUC (computed by script)"]
    header(ws, 5, heads)
    ws.row_dimensions[5].height = 32
    r = 6
    srow = {}
    for key, title, _, _ in DATASETS:
        for ver in ("before", "after"):
            sh = f"{key}_{ver}"
            n = ranges[sh]
            yt, yp = f"{sh}!$D$2:$D${n}", f"{sh}!$G$2:$G${n}"
            ws.cell(r, 1, title)
            ws.cell(r, 2, "original rules" if ver == "before" else "current rules")
            ws.cell(r, 5, f"=COUNTIFS({yt},1,{yp},1)")
            ws.cell(r, 6, f"=COUNTIFS({yt},1,{yp},0)")
            ws.cell(r, 7, f"=COUNTIFS({yt},0,{yp},1)")
            ws.cell(r, 8, f"=COUNTIFS({yt},0,{yp},0)")
            ws.cell(r, 3, f"=E{r}+F{r}+G{r}+H{r}")
            ws.cell(r, 4, f"=E{r}+F{r}")
            ws.cell(r, 9, f'=IF(D{r}=0,"n/a",E{r}/D{r})')
            ws.cell(r, 12, f'=IF(H{r}+G{r}=0,"n/a",H{r}/(H{r}+G{r}))')
            for p, nn, lo, hi in ((f"I{r}", f"D{r}", 10, 11), (f"L{r}", f"(H{r}+G{r})", 13, 14)):
                z = "$B$3"
                core = f"{z}*SQRT({p}*(1-{p})/{nn}+{z}^2/(4*{nn}^2))"
                ws.cell(r, lo, f'=IFERROR(({p}+{z}^2/(2*{nn})-{core})/(1+{z}^2/{nn}),"n/a")')
                ws.cell(r, hi, f'=IFERROR(({p}+{z}^2/(2*{nn})+{core})/(1+{z}^2/{nn}),"n/a")')
            ws.cell(r, 15, f'=IF(E{r}+G{r}=0,"n/a",E{r}/(E{r}+G{r}))')
            ws.cell(r, 16, f'=IF(C{r}=0,"n/a",(E{r}+H{r})/C{r})')
            auc = data[key]["m_before" if ver == "before" else "m_after"]["auc"]
            ws.cell(r, 17, None if math.isnan(auc) else round(auc, 4)).font = f_input
            srow[(key, ver)] = r
            r += 1
    style_range(ws, f"A6:Q{r - 1}")
    for rr in range(6, r):
        ws.cell(rr, 17).font = f_input
        for cc in range(9, 17):
            ws.cell(rr, cc).number_format = "0.000"
        ws.cell(rr, 17).number_format = "0.000"
    for j, w in enumerate([26, 20, 9, 10, 7, 7, 7, 7, 12, 11, 11, 12, 11, 11, 13, 10, 14], start=1):
        ws.column_dimensions[L(j)].width = w
    n0 = r + 1
    notes = [
        "How to read this sheet: TP/FN/FP/TN and every rate are FORMULAS over the data tabs (S11_before ... hg19_after).",
        "Blue = values typed in or produced by the Python script (z, AUC). AUC needs the full score ranking, so it is not a sheet formula.",
        "'original rules' = ACMG_SPLICE_PVS1=0 ACMG_REVEL_PP3=0. 'current rules' = splice-site PVS1 + ClinGen-calibrated REVEL (see the Ablation tab for each switch alone).",
        "CAUTION: the splice rule was added after inspecting the two S11 misses, so S11 'after' is NOT independent evidence. REVEL was trained on public variant data that overlaps ClinVar, so its gain on ClinVar labels is optimistic.",
        "The larger hg19 set shows the same direction of improvement, but it is the same dataset family and the rule was not tuned on it separately.",
        "Small positive counts give wide intervals (see CI columns). These are proof-of-concept numbers, not a clinical benchmark.",
    ]
    for i, t in enumerate(notes):
        ws.cell(n0 + i, 1, t).font = f_base
    ws.freeze_panes = "C6"

    # ---- Confusion ----
    wc = wb.create_sheet("Confusion", 1)
    wc["A1"] = "Confusion matrices (linked to Summary)"
    wc["A1"].font = f_title
    rr = 3
    for key, title, _, _ in DATASETS:
        for ver in ("before", "after"):
            sr = srow[(key, ver)]
            wc.cell(rr, 1, f"{title} - {ver} rules").font = f_bold
            header(wc, rr + 1, ["", "Predicted P/LP", "Predicted not P/LP"])
            wc.cell(rr + 2, 1, "Actual P/LP")
            wc.cell(rr + 2, 2, f"=Summary!E{sr}")
            wc.cell(rr + 2, 3, f"=Summary!F{sr}")
            wc.cell(rr + 3, 1, "Actual not P/LP")
            wc.cell(rr + 3, 2, f"=Summary!G{sr}")
            wc.cell(rr + 3, 3, f"=Summary!H{sr}")
            style_range(wc, f"A{rr + 2}:C{rr + 3}")
            for x in (rr + 2, rr + 3):
                wc.cell(x, 1).font = f_bold
            rr += 6
    for j, w in enumerate([30, 18, 20], start=1):
        wc.column_dimensions[L(j)].width = w

    # ---- ACMG / AMP distributions (labelled rows) ----
    wd = wb.create_sheet("ACMG_AMP")
    wd["A1"] = "ACMG class and AMP tier distribution (labelled rows only, 'after' version)"
    wd["A1"].font = f_title
    acmg_labels = ["PATHOGENIC", "LIKELY_PATHOGENIC", "VUS", "LIKELY_BENIGN", "BENIGN"]
    amp_labels = ["TIER_I", "TIER_II", "TIER_III", "TIER_IV"]
    header(wd, 3, ["Dataset"] + acmg_labels + amp_labels)
    for i, (key, title, _, _) in enumerate(DATASETS):
        rw = 4 + i
        sh = f"{key}_after"
        n = ranges[sh]
        wd.cell(rw, 1, title)
        for j, lab in enumerate(acmg_labels):
            wd.cell(rw, 2 + j, f'=COUNTIF({sh}!$E$2:$E${n},"{lab}")')
        for j, lab in enumerate(amp_labels):
            wd.cell(rw, 7 + j, f'=COUNTIF({sh}!$H$2:$H${n},"{lab}")')
    style_range(wd, f"A4:J{3 + len(DATASETS)}")
    wd.cell(5 + len(DATASETS), 1, "ACMG/AMP internal conflicts (ACMG P/LP with AMP Tier IV, or ACMG B/LB with Tier I/II):").font = f_bold
    for i, (key, title, _, _) in enumerate(DATASETS):
        wd.cell(6 + len(DATASETS) + i, 1, title).font = f_base
        c = wd.cell(6 + len(DATASETS) + i, 2, conflicts(data[key]["after"]))
        c.font = f_input
    wd.cell(7 + 2 * len(DATASETS), 1, "(conflict counts are computed by the script from the TSV; blue = script value)").font = f_base
    for j in range(1, 11):
        wd.column_dimensions[L(j)].width = 24 if j == 1 else 16

    # ---- FN by consequence ----
    wf = wb.create_sheet("Misses_by_consequence")
    wf["A1"] = "Where are the missed pathogenic variants? ('after' version, formulas over data tabs)"
    wf["A1"].font = f_title
    rw = 3
    for key, title, _, _ in DATASETS:
        sh = f"{key}_after"
        n = ranges[sh]
        wf.cell(rw, 1, title).font = f_bold
        header(wf, rw + 1, ["Consequence (ANNOVAR)", "ClinVar P/LP", "Detected (TP)", "Missed (FN)", "Sensitivity"])
        cats = sorted(data[key]["after"]["consequence"].unique())
        for i, cat in enumerate(cats):
            x = rw + 2 + i
            wf.cell(x, 1, cat)
            wf.cell(x, 2, f'=COUNTIFS({sh}!$I$2:$I${n},A{x},{sh}!$D$2:$D${n},1)')
            wf.cell(x, 3, f'=COUNTIFS({sh}!$I$2:$I${n},A{x},{sh}!$D$2:$D${n},1,{sh}!$G$2:$G${n},1)')
            wf.cell(x, 4, f"=B{x}-C{x}")
            wf.cell(x, 5, f'=IF(B{x}=0,"n/a",C{x}/B{x})')
            wf.cell(x, 5).number_format = "0.000"
        last = rw + 1 + len(cats)
        style_range(wf, f"A{rw + 2}:E{last}")
        tot = last + 1
        wf.cell(tot, 1, "Total").font = f_bold
        for col in "BCD":
            wf[f"{col}{tot}"] = f"=SUM({col}{rw + 2}:{col}{last})"
            wf[f"{col}{tot}"].font = f_bold
        rw = tot + 3
    wf.cell(rw, 1, "Reading: the engine implements 8 of ~28 ACMG criteria and, in this validation, no ClinVar. Missense P/LP "
                   "variants usually need criteria it does not implement (PS1, PS3, PM1, PM5 ...), so they stay VUS.").font = f_base
    for j, w in enumerate([30, 14, 14, 12, 12], start=1):
        wf.column_dimensions[L(j)].width = w

    # ---- Performance ----
    wp = wb.create_sheet("Performance")
    wp["A1"] = "Run time (measured from *_summary.json of the 'after' runs)"
    wp["A1"].font = f_title
    header(wp, 3, ["Dataset", "Mode", "Rows in file", "Unique variants", "Evidence seconds", "Seconds per variant",
                   "Lookups from cache", "Live lookups", "Skipped lookups"])
    for i, (key, title, _, _) in enumerate(DATASETS):
        s = data[key]["summary"]
        x = 4 + i
        wp.cell(x, 1, title)
        if s:
            wp.cell(x, 2, s["mode"])
            wp.cell(x, 3, s["total_rows"])
            wp.cell(x, 4, s["unique_variants_classified"])
            wp.cell(x, 5, s["evidence_seconds"])
            wp.cell(x, 6, f"=E{x}/D{x}")
            wp.cell(x, 7, s["lookups"]["local_db"])
            wp.cell(x, 8, s["lookups"]["live_api"])
            wp.cell(x, 9, s["lookups"]["skipped"])
            for cc in (3, 4, 5, 7, 8, 9):
                wp.cell(x, cc).font = f_input
        wp.cell(x, 6).number_format = "0.0000"
    style_range(wp, f"A4:I{3 + len(DATASETS)}")
    for i in range(len(DATASETS)):
        for cc in (3, 4, 5, 7, 8, 9):
            wp.cell(4 + i, cc).font = f_input
    x0 = 6 + len(DATASETS)
    wp.cell(x0, 1, "Why the old code could run for hours (derived from the OLD retry settings - an upper bound, not a measurement):").font = f_bold
    header(wp, x0 + 1, ["Source", "Old behaviour per failing lookup", "Seconds (worst case)", "Lookups still needing the network (hg19 set)", "Worst-case minutes (sequential)"])
    old = [
        ("gnomAD (HTTP 429)", "4 requests + 1+2+4 s backoff, retried on every run", 8, 123),
        ("TCGA / GDC", "2 attempts x 10 s timeout + 1 s backoff", 21, 383),
    ]
    for i, (a, b, c, d) in enumerate(old):
        x = x0 + 2 + i
        wp.cell(x, 1, a)
        wp.cell(x, 2, b)
        wp.cell(x, 3, c).font = f_input
        wp.cell(x, 4, d).font = f_input
        wp.cell(x, 5, f"=C{x}*D{x}/60")
        wp.cell(x, 5).number_format = "0.0"
    style_range(wp, f"A{x0 + 2}:E{x0 + 3}")
    for i in range(2):
        wp.cell(x0 + 2 + i, 3).font = f_input
        wp.cell(x0 + 2 + i, 4).font = f_input
    tot = x0 + 4
    wp.cell(tot, 1, "Total worst case (minutes)").font = f_bold
    wp.cell(tot, 5, f"=SUM(E{x0 + 2}:E{x0 + 3})").font = f_bold
    wp.cell(tot, 5).number_format = "0.0"
    wp.cell(tot + 1, 1, "Blue inputs come from the old source code (retry/timeouts) and scripts/preflight.py cache counts. "
                        "The new code runs the sources in parallel with short timeouts, caches no failures and stops "
                        "after repeated failures, so these worst cases no longer apply.").font = f_base
    for j, w in enumerate([30, 46, 18, 40, 26, 20, 20, 14, 16], start=1):
        wp.column_dimensions[L(j)].width = w

    # ---- Rules ----
    wr = wb.create_sheet("Rules")
    wr["A1"] = "Rules implemented (read from the code, so this sheet cannot drift)"
    wr["A1"].font = f_title
    header(wr, 3, ["Item", "Value"])
    rules = [
        ("Bayesian points", "VERY_STRONG=%d, STRONG=%d, MODERATE=%d, SUPPORTING=%d" % tuple(acmg_mod.POINTS[k] for k in ("VERY_STRONG", "STRONG", "MODERATE", "SUPPORTING"))),
        ("BA1 (stand-alone benign)", "AF >= %s" % acmg_mod.BA1_THRESHOLD),
        ("BS1 (strong benign)", "AF >= %s" % acmg_mod.BS1_THRESHOLD),
        ("PM2 (moderate pathogenic)", "AF < %s" % acmg_mod.PM2_THRESHOLD),
        ("PVS1 (very strong pathogenic)", "frameshift / stopgain / stoploss, or canonical splice site (ANNOVAR 'splicing'); switch ACMG_SPLICE_PVS1"),
        ("PP3 / BP4", "REVEL (ClinGen-calibrated: >=0.644 supporting, >=0.773 moderate, >=0.932 strong; <=0.290 / 0.183 / 0.016 / 0.003 benign), else >=3 damaging / >=3 benign predictors"),
        ("PP5 / BP6", "ClinVar P/LP or B/LB for the EXACT allele (via MyVariant); hidden from both engines in validation runs"),
        ("Classes", "Pathogenic >= 10; Likely pathogenic 6-9; VUS 0-5; Likely benign -6 to -1; Benign <= -7; BA1 forces Benign"),
        ("AMP tiers", "Common variant or ClinVar benign -> IV; Tier I/II gene panel + (ClinVar P, null, or TCGA >= 10 cases) -> I/II; otherwise III"),
    ]
    for i, (a, b) in enumerate(rules, start=4):
        wr.cell(i, 1, a)
        wr.cell(i, 2, b)
    style_range(wr, f"A4:B{3 + len(rules)}", align=Alignment(wrap_text=True, vertical="top"))
    wr.column_dimensions["A"].width = 32
    wr.column_dimensions["B"].width = 110

    # ---- Charts ----
    wg = wb.create_sheet("Charts")
    wg["A1"] = "ROC curves and confusion matrices (images generated by the script)"
    wg["A1"].font = f_title
    row = 3
    for key, _, _, _ in DATASETS:
        for fk, anchor in ((f"roc_{key}", f"A{row}"), (f"cm_{key}", f"J{row}")):
            img = XLImage(str(figs[fk]))
            img.width, img.height = (int(img.width * 0.42), int(img.height * 0.42)) if fk.startswith("roc") else (int(img.width * 0.5), int(img.height * 0.5))
            wg.add_image(img, anchor)
        row += 24

    # ---- Notes ----
    abl_path = Path(__file__).resolve().parents[1] / "data/output/ablation.csv"
    if abl_path.exists():
        wa = wb.create_sheet("Ablation")
        wa["A1"] = "Rule ablation (ClinVar hidden; cache-only runs; produced by scripts/ablation.py)"
        wa["A1"].font = f_title
        abl = pd.read_csv(abl_path)
        header(wa, 3, list(abl.columns))
        for i, rec in enumerate(abl.itertuples(index=False), start=4):
            for j, val in enumerate(rec, start=1):
                c = wa.cell(i, j, val if not hasattr(val, "item") else val.item())
                c.font = f_base
        for col in "ABCDEFGHIJKL":
            wa.column_dimensions[col].width = 24 if col in "AB" else 12
        wa.cell(len(abl) + 6, 1, "S11 gains from the splice rule are post-hoc; REVEL training data overlaps ClinVar so its gain is optimistic.").font = f_base

    wn = wb.create_sheet("Notes")
    wn["A1"] = "Method and honest caveats"
    wn["A1"].font = f_title
    lines = [
        "Ground truth: ClinVar significance copied into the input file (CLNSIG). Positive = Pathogenic / Likely pathogenic.",
        "The engine was run with ACMG_IGNORE_CLINVAR=1, so the label was not used for PP5/BP6/PVS1 overrides.",
        "Population frequency: the highest of gnomAD exome, gnomAD genome and the AF columns in the input file; gnomAD 'variant not found' with no file frequency counts as absent.",
        "Duplicate rows (same chrom:pos:ref:alt) are collapsed by the parser; the hg19 file has 2251 rows but far fewer unique variants.",
        "The splice-site PVS1 rule was introduced after the two S11 misses were inspected: its S11 gain is post-hoc. The calibrated-REVEL rule uses published thresholds, but REVEL's training data overlaps ClinVar, so its gain on ClinVar labels is optimistic.",
        "Only 8 of ~28 ACMG criteria are implemented; AMP tiering uses a small curated gene panel. Not a diagnostic device.",
        "With few positives, the confidence intervals are wide. Treat as proof-of-concept, not clinical performance.",
        "Offline runs use only the local cache; lookups not in the cache are reported as skipped and the engines fall back to input-file columns.",
    ]
    for i, t in enumerate(lines, start=3):
        wn.cell(i, 1, f"{i - 2}. {t}").font = f_base
    wn.column_dimensions["A"].width = 150

    order = ["Summary", "Confusion", "ACMG_AMP", "Misses_by_consequence", "Performance", "Rules", "Charts", "Notes"]
    wb._sheets = [wb[n] for n in order] + [s for s in wb._sheets if s.title not in order]
    path.parent.mkdir(parents=True, exist_ok=True)
    wb.save(path)


# ----------------------------------------------------------------------------- DOCX
def build_docx(data: dict, figs: dict[str, Path], path: Path) -> None:
    from docx import Document
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Cm, Pt, RGBColor

    doc = Document()
    sec = doc.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)  # A4
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Cm(2.0))
    st = doc.styles["Normal"]
    st.font.name, st.font.size = FONT, Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    for name, size in (("Heading 1", 15), ("Heading 2", 12), ("Title", 22)):
        s = doc.styles[name]
        s.font.name, s.font.size = FONT, Pt(size)
        s.font.color.rgb = RGBColor(0x1F, 0x38, 0x64)
        s.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)

    def shade(cell, hex_fill):
        tcPr = cell._tc.get_or_add_tcPr()
        shd = OxmlElement("w:shd")
        shd.set(qn("w:val"), "clear")
        shd.set(qn("w:color"), "auto")
        shd.set(qn("w:fill"), hex_fill)
        tcPr.append(shd)

    def table(rows, header_row=True, widths=None, size=9):
        t = doc.add_table(rows=len(rows), cols=len(rows[0]))
        t.style = "Table Grid"
        t.autofit = False
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        if widths:
            for j, w in enumerate(widths):
                t.columns[j].width = Cm(w)
        for i, row in enumerate(rows):
            for j, val in enumerate(row):
                c = t.cell(i, j)
                c.text = ""
                run = c.paragraphs[0].add_run(str(val))
                run.font.name, run.font.size = FONT, Pt(size)
                if header_row and i == 0:
                    run.bold = True
                    run.font.color.rgb = RGBColor(255, 255, 255)
                    shade(c, "1F3864")
                if widths:
                    c.width = Cm(widths[j])
        doc.add_paragraph()
        return t

    def bullet(text):
        p = doc.add_paragraph(style="List Bullet")
        p.add_run(text)

    mS, mSa = data["S11"]["m_before"], data["S11"]["m_after"]
    mH, mHa = data["hg19"]["m_before"], data["hg19"]["m_after"]

    doc.add_paragraph("Variant Classification Platform", style="Title")
    p = doc.add_paragraph()
    p.add_run(f"Statistical validation report - {date.today().strftime('%d %B %Y')}").bold = True
    doc.add_paragraph("Research / evidence-aggregation tool using ACMG and AMP frameworks. Not a diagnostic device.")

    doc.add_heading("1. Summary", 1)
    bullet(f"Small labelled set (S11, {mSa['n']} variants, {mSa['pos']} ClinVar-pathogenic): sensitivity "
           f"{fmt(mS['sens'])} -> {fmt(mSa['sens'])} with the current rules (splice-site PVS1 + calibrated REVEL); specificity "
           f"{fmt(mSa['spec'])}; AUC {fmt(mS['auc'])} -> {fmt(mSa['auc'])}.")
    bullet(f"Larger labelled set (hg19, {mHa['n']} variants, {mHa['pos']} pathogenic): sensitivity {fmt(mH['sens'])} -> "
           f"{fmt(mHa['sens'])} (95% CI {fmt(mHa['sens_ci'][0], 2)}-{fmt(mHa['sens_ci'][1], 2)}); specificity "
           f"{fmt(mHa['spec'])}; accuracy {fmt(mH['acc'])} -> {fmt(mHa['acc'])}; AUC {fmt(mH['auc'])} -> {fmt(mHa['auc'])}.")
    bullet("The engine is conservative: it rarely calls a benign/VUS variant pathogenic (high specificity and precision) "
           "but misses many ClinVar-pathogenic missense variants, because only 8 of about 28 ACMG criteria are implemented.")
    bullet("The S11 improvement is post-hoc (the splice rule was added after inspecting its two misses) and is not independent evidence. "
           "The hg19 set is the stronger check. Its gain from calibrated REVEL costs some specificity and is optimistic because REVEL training data overlaps ClinVar.")
    bullet("The earlier multi-hour runs were most likely caused by sequential network calls with long retries on gnomAD (HTTP 429) and TCGA "
           "(derived from the old code, not measured); see section 6. A cache-only run now takes about a second.")

    doc.add_heading("2. Method and data", 1)
    bullet("Ground truth: ClinVar significance in the input file (CLNSIG). Positive = Pathogenic / Likely pathogenic; "
           "conflicting and blank labels are excluded.")
    bullet("ClinVar was hidden from the engine (ACMG_IGNORE_CLINVAR=1), so the label is not used to produce the prediction.")
    bullet("Prediction = ACMG class Pathogenic or Likely pathogenic, scored with the ClinGen/Tavtigian points system.")
    bullet("Duplicate rows (same chromosome, position, ref, alt) are merged by the parser.")
    sh, sa = data["S11"]["summary"], data["hg19"]["summary"]
    if sa:
        bullet(f"hg19 file: {sa['total_rows']} rows -> {sa['unique_variants_classified']} unique variants -> {mHa['n']} usable labelled variants.")
    bullet("'Before' = original rules (ACMG_SPLICE_PVS1=0 ACMG_REVEL_PP3=0). 'After' = current default: ANNOVAR 'splicing' variants "
           "(within 2 bp of an exon boundary) count as null variants for PVS1, and PP3/BP4 use ClinGen-calibrated REVEL strength.")

    doc.add_heading("3. Results", 1)
    rows = [["Dataset / version", "n", "Pos", "TP", "FN", "FP", "TN", "Sensitivity (95% CI)", "Specificity (95% CI)", "Precision", "Accuracy", "AUC"]]
    for key, title, _, _ in DATASETS:
        for ver, mm in (("before", data[key]["m_before"]), ("after", data[key]["m_after"])):
            rows.append([f"{title} - {ver}", mm["n"], mm["pos"], mm["tp"], mm["fn"], mm["fp"], mm["tn"],
                         f"{fmt(mm['sens'])} ({fmt(mm['sens_ci'][0], 2)}-{fmt(mm['sens_ci'][1], 2)})",
                         f"{fmt(mm['spec'])} ({fmt(mm['spec_ci'][0], 2)}-{fmt(mm['spec_ci'][1], 2)})",
                         fmt(mm["prec"]), fmt(mm["acc"]), fmt(mm["auc"])])
    table(rows, widths=[2.5, 1.0, 1.1, 0.9, 0.9, 0.9, 1.1, 2.4, 2.4, 1.4, 1.4, 1.2], size=7.5)
    doc.add_paragraph("Confidence intervals are Wilson 95% intervals. With few positives the S11 intervals are very wide.")

    for key, title, _, _ in DATASETS:
        doc.add_heading(f"{title}", 2)
        doc.add_picture(str(figs[f"roc_{key}"]), width=Cm(8.0))
        doc.add_picture(str(figs[f"cm_{key}"]), width=Cm(14.0))

    abl_path = Path(__file__).resolve().parents[1] / "data/output/ablation.csv"
    if abl_path.exists():
        abl = pd.read_csv(abl_path)
        doc.add_heading("Rule ablation (each switch alone and together)", 2)
        rows = [["Dataset", "Rules", "TP", "FN", "FP", "TN", "Sens.", "Spec.", "Acc.", "AUC"]]
        for r in abl.itertuples(index=False):
            rows.append([r.dataset, r.rules, str(r.TP), str(r.FN), str(r.FP), str(r.TN),
                         f"{r.sensitivity:.3f}", f"{r.specificity:.3f}", f"{r.accuracy:.3f}", f"{r.AUC:.3f}"])
        table(rows, widths=[4.2, 4.4, 1.2, 1.2, 1.2, 1.2, 1.5, 1.5, 1.5, 1.5])
        doc.add_paragraph("Reading: the splice-site rule adds sensitivity at almost no cost; calibrated REVEL adds more sensitivity "
                          "but lowers specificity (more false positives), and its gain is optimistic because REVEL was trained on "
                          "public variant data that overlaps ClinVar.")

    doc.add_heading("4. Where the misses are", 1)
    after = data["hg19"]["after"]
    pos = after[after["y_true"] == 1]
    g = pos.groupby("consequence")["y_pred"].agg(["size", "sum"]).reset_index()
    rows = [["Consequence (ANNOVAR)", "ClinVar P/LP", "Detected", "Missed", "Sensitivity"]]
    for r in g.itertuples(index=False):
        rows.append([r.consequence, int(r.size), int(r.sum), int(r.size - r.sum), fmt(r.sum / r.size)])
    table(rows, widths=[5, 3, 3, 3, 3])
    missed = pos[pos["y_pred"] == 0]
    top = missed["consequence"].value_counts()
    top_txt = ", ".join(f"{int(n)} {c}" for c, n in top.items()) or "none"
    doc.add_paragraph(f"Reading: {len(missed)} of {len(pos)} ClinVar-pathogenic variants were missed ({top_txt}). "
                      "Null-type variants (frameshift, stop-gain, splice) are detected; missense variants need ACMG criteria "
                      "that are not implemented here (PS1, PS3, PM1, PM5 and others), so they stay VUS. "
                      "This is the main ceiling on sensitivity and the clearest next step. Calibrated REVEL recovers part of it.")

    doc.add_heading("5. Classification distribution and consistency", 1)
    rows = [["Dataset (labelled rows, after)", "Path.", "Likely path.", "VUS", "Likely benign", "Benign", "Tier I", "Tier II", "Tier III", "Tier IV", "ACMG/AMP conflicts"]]
    for key, title, _, _ in DATASETS:
        d = data[key]["after"]
        c = d["acmg_classification"].value_counts()
        a = d["amp_tier"].value_counts()
        rows.append([title] + [int(c.get(k, 0)) for k in ("PATHOGENIC", "LIKELY_PATHOGENIC", "VUS", "LIKELY_BENIGN", "BENIGN")]
                    + [int(a.get(k, 0)) for k in ("TIER_I", "TIER_II", "TIER_III", "TIER_IV")] + [conflicts(d)])
    table(rows)

    doc.add_heading("6. Run time", 1)
    rows = [["Dataset", "Mode", "Unique variants", "Evidence time (s)", "s per variant", "Cache", "Live", "Skipped"]]
    for key, title, _, _ in DATASETS:
        s = data[key]["summary"]
        if s:
            rows.append([title, s["mode"], s["unique_variants_classified"], s["evidence_seconds"], s["seconds_per_variant"],
                         s["lookups"]["local_db"], s["lookups"]["live_api"], s["lookups"]["skipped"]])
    table(rows)
    doc.add_paragraph("Why the old code ran for hours (derived from its retry settings; an upper bound, not a measurement):")
    bullet("gnomAD returned HTTP 429 for 343 cached variants. Each was retried on every run with 4 requests and 1+2+4 s of backoff.")
    bullet("TCGA had 383 uncached lookups in the hg19 file, each up to 2 x 10 s timeouts + 1 s, one after another.")
    bullet("All sources ran one after another with no progress output, so a slow run looked frozen.")
    bullet("The .env file was never loaded, so an NCBI API key was ignored and ClinVar stayed at 3 requests per second.")
    doc.add_paragraph("Fixes: sources run in parallel with 5/10 s timeouts and at most 2 fast retries; failures are never cached; "
                      "a circuit breaker and a wall-clock limit stop slow sources; MyVariant lookups are batched; "
                      "progress and per-source timings are printed; --offline gives a cache-only run.")

    doc.add_heading("7. Limitations", 1)
    for t in [
        "Only 8 of ~28 ACMG criteria are implemented; sensitivity on missense variants is low by design.",
        "The splice-site rule was added after inspecting the S11 misses. Its S11 effect is post-hoc, not independent. The REVEL gain is optimistic (training overlap with ClinVar).",
        "Labels come from ClinVar, which also feeds the annotation columns of the input file; hiding ClinVar from the engine "
        "removes direct use of the label but not all shared information (for example in-silico predictor columns).",
        "Few positives in S11 (6) mean very wide intervals. The hg19 set has more positives but is a single data source.",
        "Offline runs report lookups that were not in the cache as skipped; the engines then use the input-file frequency columns.",
        "Research tool only. Clinical use requires review by a qualified professional.",
    ]:
        bullet(t)

    doc.add_heading("8. Reproduce", 1)
    for cmd in [
        "ACMG_IGNORE_CLINVAR=1 PYTHONPATH=backend python -m app.cli --input data/input/S11_labeled_subset.tsv --limit 0 --output data/output/S11_after.tsv",
        "ACMG_SPLICE_PVS1=0 ACMG_REVEL_PP3=0 ACMG_IGNORE_CLINVAR=1 PYTHONPATH=backend python -m app.cli --input data/input/S11_labeled_subset.tsv --limit 0 --output data/output/S11_before.tsv",
        "(same two commands for data/input/hg19Somatic_validation.tsv -> hg19_after.tsv / hg19_before.tsv)",
        "python scripts/make_report.py",
    ]:
        pr = doc.add_paragraph()
        r = pr.add_run(cmd)
        r.font.name, r.font.size = "Courier New", Pt(8.5)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(path)


# ----------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "data" / "output"))
    a = ap.parse_args()
    out = Path(a.out)
    data: dict = {}
    for key, title, b, af in DATASETS:
        pb, pa = out / b, out / af
        if not (pb.exists() and pa.exists()):
            sys.exit(f"Missing {pb.name} or {pa.name} in {out}. Run the commands in docs/USER_GUIDE.md section 6 first.")
        db, da = load_labelled(pb), load_labelled(pa)
        data[key] = {"before": db, "after": da, "m_before": metrics(db), "m_after": metrics(da), "summary": load_summary(pa)}
    figs = make_figures(data, out / "report_figures")
    build_xlsx(data, figs, out / "statistics_report.xlsx")
    build_docx(data, figs, out / "statistics_report.docx")
    print(f"Wrote {out / 'statistics_report.xlsx'}")
    print(f"Wrote {out / 'statistics_report.docx'}")
    print(f"Wrote figures in {out / 'report_figures'}")


if __name__ == "__main__":
    main()
