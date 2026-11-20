"""Shared definition of a 'null' (loss-of-function) variant for the ACMG and AMP engines.

Two kinds are recognised:
  * exonic null consequences from ANNOVAR's ExonicFunc.refGene (frameshift, stopgain, stoploss)
  * canonical splice-site variants: ANNOVAR Func.refGene == "splicing" (within 2 bp of an
    exon/intron boundary). Rule PVS1 for these is switchable with ACMG_SPLICE_PVS1=0.

IMPORTANT (honesty note): the splice rule was added after a validation run showed that the only
two ClinVar-pathogenic variants missed were canonical splice-site variants. It is clinically
reasonable (ClinGen PVS1 covers canonical splice sites) but it was found by looking at those
misses, so the improved sensitivity on that same set is NOT independent evidence.
"""

from __future__ import annotations

import os

EXONIC_NULL = {"frameshift_insertion", "frameshift_deletion", "frameshift_substitution",
               "stopgain", "stoploss"}


def splice_pvs1_enabled() -> bool:
    return os.environ.get("ACMG_SPLICE_PVS1", "1") != "0"


def null_variant_kind(local_annotation: dict) -> str | None:
    """Returns the null-variant description, or None when the variant is not null."""
    exonic = (local_annotation.get("ExonicFunc.refGene") or "").strip().lower()
    if exonic in EXONIC_NULL:
        return exonic
    func = (local_annotation.get("Func.refGene") or "").strip().lower()
    if splice_pvs1_enabled() and func == "splicing":
        return "canonical splice-site (ANNOVAR splicing, within 2 bp of exon boundary)"
    return None
