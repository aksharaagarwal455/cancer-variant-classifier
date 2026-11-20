"""
Validation + consistency report for classifier output.

Usage:
  python scripts/validate.py data/output/sample_S11_results.tsv

1. Internal conflicts: ACMG says (likely) pathogenic but AMP says Tier IV (or the reverse).
2. Agreement with the CLNSIG column that came with the input file.
   Note: PP5 uses ClinVar, so this is agreement with ClinVar, not independent accuracy.
"""
import sys
import pandas as pd

path = sys.argv[1] if len(sys.argv) > 1 else "data/output/sample_S11_results.tsv"
df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
print(f"File: {path}   rows: {len(df)}\n")


def acmg_group(v):
    v = v.upper()
    if "PATHOGENIC" in v and "BENIGN" not in v:
        return "P/LP"
    if "BENIGN" in v:
        return "B/LB"
    if "VUS" in v or "UNCERTAIN" in v:
        return "VUS"
    return "other"


def truth_group(v):
    s = v.lower().replace("_", " ")
    if not s or "conflicting" in s:
        return None
    if "pathogenic" in s and "benign" not in s:
        return "P/LP"
    if "benign" in s and "pathogenic" not in s:
        return "B/LB"
    if "uncertain" in s or "vus" in s:
        return "VUS"
    return None


df["acmg_group"] = df["acmg_classification"].map(acmg_group)

print("=== 1. INTERNAL CONFLICTS (ACMG vs AMP) ===")
conf = df[((df["acmg_group"] == "P/LP") & (df["amp_tier"] == "TIER_IV")) |
          ((df["acmg_group"] == "B/LB") & df["amp_tier"].isin(["TIER_I", "TIER_II"]))]
print(f"{len(conf)} of {len(df)} variants have conflicting ACMG/AMP calls")
for _, r in conf.head(10).iterrows():
    print(f"- {r['variant_id']}  gene={r.get('Gene','')}  CLNSIG={r.get('CLNSIG','')}")
    print(f"    ACMG: {r['acmg_classification']} | {r['acmg_criteria']}")
    print(f"    AMP : {r['amp_tier']} | {r['amp_explanation']}")

print("\n=== 2. AGREEMENT WITH CLNSIG (input file) ===")
if "CLNSIG" not in df.columns:
    print("No CLNSIG column in this file.")
else:
    df["truth"] = df["CLNSIG"].map(truth_group)
    t = df[df["truth"].notna()]
    print(f"{len(t)} of {len(df)} rows have a usable CLNSIG label (conflicting/empty skipped)")
    if len(t):
        print(pd.crosstab(t["truth"], t["acmg_group"], rownames=["CLNSIG"], colnames=["ACMG"]))
        agree = (t["truth"] == t["acmg_group"]).mean() * 100
        print(f"\nExact agreement (P/LP vs VUS vs B/LB): {agree:.1f}%")
        pl = t[t["truth"] == "P/LP"]
        if len(pl):
            print(f"ClinVar P/LP recovered as P/LP by ACMG: {(pl['acmg_group'] == 'P/LP').mean() * 100:.1f}% ({len(pl)} variants)")
