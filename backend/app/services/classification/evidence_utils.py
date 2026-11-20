"""Evidence helpers shared by the ACMG and AMP engines (so both read evidence the same way)."""

from __future__ import annotations

import os

# -------- ClinVar ---------------------------------------------------------------------------------
def clinvar_ignored() -> bool:
    """ACMG_IGNORE_CLINVAR=1 hides ClinVar from BOTH engines (used for independent validation)."""
    return os.environ.get("ACMG_IGNORE_CLINVAR") == "1"


def clinvar_records(evidence) -> list[dict] | None:
    """Records of the exact-allele ClinVar result, or None when ClinVar is hidden/unavailable."""
    if clinvar_ignored():
        return None
    cv = getattr(evidence, "clinvar", None)
    if not (cv and cv.status.value == "AVAILABLE"):
        return None
    return cv.data.get("records", []) or []


# -------- population frequency --------------------------------------------------------------------
LOCAL_AF_COLUMNS = ("AF", "AF_popmax", "AF_raw", "ExAC_ALL")


def _float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def population_af(evidence) -> tuple[float | None, str]:
    """Highest allele frequency among all available sources (the conservative choice).

    Sources: gnomAD exome AND genome (the old code looked at genome only and silently ignored
    exome-only variants), then the frequency columns already present in the input file.
    A variant that gnomAD explicitly reports as 'Variant not found' and has no file frequency is
    treated as absent (AF = 0), which is what PM2 asks for.
    """
    candidates: list[tuple[float, str]] = []
    gn = getattr(evidence, "gnomad", None)
    gnomad_absent = False
    if gn and gn.status.value == "AVAILABLE":
        v = gn.data.get("variant")
        if v:
            for kind in ("exome", "genome"):
                af = _float((v.get(kind) or {}).get("af"))
                if af is not None:
                    candidates.append((af, f"gnomAD {kind}"))
        else:
            msgs = [str(e.get("message", "")) for e in (gn.data.get("errors") or [])]
            gnomad_absent = bool(msgs) and all("not found" in m.lower() for m in msgs) or (not msgs and gn.data.get("found") is False)
    for col in LOCAL_AF_COLUMNS:
        af = _float(evidence.local_annotation.get(col))
        if af is not None:
            candidates.append((af, f"file {col}"))
    if candidates:
        return max(candidates, key=lambda t: t[0])
    if gnomad_absent:
        return 0.0, "gnomAD (variant not found = absent)"
    return None, "no population frequency data"


# -------- REVEL (ClinGen SVI calibration, Pejaver et al. 2022) ------------------------------------
REVEL_PATHOGENIC = ((0.932, "STRONG"), (0.773, "MODERATE"), (0.644, "SUPPORTING"))
REVEL_BENIGN = ((0.003, "VERY_STRONG"), (0.016, "STRONG"), (0.183, "MODERATE"), (0.290, "SUPPORTING"))


def revel_strength(score) -> tuple[str, str] | None:
    """('pathogenic'|'benign', strength) for a calibrated REVEL score, or None if indeterminate."""
    s = _float(score)
    if s is None:
        return None
    for cut, strength in REVEL_PATHOGENIC:
        if s >= cut:
            return "pathogenic", strength
    for cut, strength in REVEL_BENIGN:
        if s <= cut:
            return "benign", strength
    return None
