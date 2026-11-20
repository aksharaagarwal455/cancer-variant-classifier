"""Measure the real time per variant for each source, bypassing the cache.

  python scripts/benchmark.py --input data/input/sample_S12.tsv --n 10

Sends N live requests per source (one at a time, so the numbers are per-call latency), prints
min / median / max seconds and how many failed. Nothing is written to the cache.
"""
from __future__ import annotations

import argparse
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import httpx

from app import config
from app.services.normalization.coordinates import apply_genome_build
from app.services.parser.tsv_parser import ParseStatus, parse_tsv
from app.services.sources.base import make_timeout
from app.services.sources.clinvar import ClinVarClient
from app.services.sources.gnomad import _BUILD_TO_DATASET, GnomadClient, build_gnomad_variant_id
from app.services.sources.myvariant import _BUILD_TO_ASSEMBLY, MyVariantClient, resolve_hgvs_genomic_id
from app.services.sources.tcga import fetch_tcga_cases, parse_protein_change


def timed(fn):
    t0 = time.monotonic()
    try:
        status, data = fn()
    except Exception as e:
        status, data = "SOURCE_UNAVAILABLE", {"error": type(e).__name__}
    return time.monotonic() - t0, status == "AVAILABLE", data.get("error", "")


def report(name: str, rows: list[tuple[float, bool, str]]) -> None:
    if not rows:
        print(f"{name:<10} no variants applicable")
        return
    ts = sorted(r[0] for r in rows)
    fails = [r for r in rows if not r[1]]
    errs = sorted({r[2] for r in fails})
    print(f"{name:<10} n={len(rows):<3} min={ts[0]:.2f}s  median={statistics.median(ts):.2f}s  max={ts[-1]:.2f}s  "
          f"failed={len(fails)} {errs if errs else ''}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data/input/sample_S12.tsv")
    ap.add_argument("--n", type=int, default=10)
    a = ap.parse_args()
    variants, _ = parse_tsv(a.input)
    ok = [v for v in variants if v.parse_status == ParseStatus.OK]
    ok, _ = apply_genome_build(ok)
    sample = ok[: a.n]
    build = sample[0].genome_build
    print(f"Benchmarking {len(sample)} variants from {a.input} (build {build.value}); live calls, cache bypassed\n")

    with httpx.Client(timeout=make_timeout(15)) as c:
        cv, mv, gn, tc = ClinVarClient(), MyVariantClient(), GnomadClient(), []
        report("clinvar", [timed(lambda v=v: cv._live_fetch(v, c)) for v in sample])
        rows = []
        for v in sample:
            h, _ = resolve_hgvs_genomic_id(v)
            if h and build in _BUILD_TO_ASSEMBLY:
                rows.append(timed(lambda h=h: mv._live_fetch(h, _BUILD_TO_ASSEMBLY[build], c)))
        report("myvariant", rows)
        rows = []
        for v in sample:
            g = build_gnomad_variant_id(v)
            if g and build in _BUILD_TO_DATASET:
                rows.append(timed(lambda v=v, g=g: gn._live_fetch(v, g, _BUILD_TO_DATASET[build], c)))
                time.sleep(1.0)
        report("gnomad", rows)
        rows = []
        for v in sample:
            aa = parse_protein_change(v.hgvsp)
            if v.gene and aa:
                rows.append(timed(lambda v=v, aa=aa: fetch_tcga_cases(v.gene, aa, c)))
        report("tcga", rows)
    print("\nA full uncached run costs roughly (median x variants) per source; the sources run in parallel,\n"
          "so the wall time is about the slowest source. Cached variants cost ~0.")


if __name__ == "__main__":
    main()
