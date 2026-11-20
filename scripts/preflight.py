"""Pre-flight check: run this BEFORE a long run or a demo.

  python scripts/preflight.py --input data/input/hg19Somatic_validation.tsv

Reports: Python version, .env / NCBI key, Redis, how many lookups per source are already
cached vs would need the network, whether each API answers (with latency), and an ESTIMATE of
how long the live lookups would take. Nothing is classified and nothing is written to the cache.
"""
from __future__ import annotations

import argparse
import sqlite3
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

from app import config  # noqa: E402,F401
import os  # noqa: E402
import httpx  # noqa: E402

from app.services.normalization.coordinates import apply_genome_build  # noqa: E402
from app.services.parser.tsv_parser import ParseStatus, parse_tsv  # noqa: E402
from app.services.sources.clinvar import _variant_key  # noqa: E402
from app.services.sources.gnomad import _BUILD_TO_DATASET, build_gnomad_variant_id  # noqa: E402
from app.services.sources.myvariant import _BUILD_TO_ASSEMBLY, resolve_hgvs_genomic_id  # noqa: E402
from app.services.sources.tcga import parse_protein_change  # noqa: E402

PINGS = {
    "myvariant": ("GET", "https://myvariant.info/v1/metadata", None),
    "clinvar": ("GET", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?db=clinvar&retmode=json", None),
    "gnomad": ("POST", "https://gnomad.broadinstitute.org/api", {"query": "{__typename}"}),
    "tcga": ("GET", "https://api.gdc.cancer.gov/status", None),
}


def ping(name: str) -> tuple[bool, float, str]:
    method, url, body = PINGS[name]
    times = []
    last = ""
    for _ in range(3):
        t0 = time.monotonic()
        try:
            r = httpx.request(method, url, json=body, timeout=httpx.Timeout(10.0, connect=5.0))
            times.append(time.monotonic() - t0)
            last = f"HTTP {r.status_code}"
            if r.status_code != 200:
                return False, statistics.median(times), last
        except Exception as e:  # noqa: BLE001
            return False, time.monotonic() - t0, f"{type(e).__name__}"
    return True, statistics.median(times), last


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="data/input/sample_S11.tsv")
    ap.add_argument("--db-path", default="data/variant_classifier.db")
    ap.add_argument("--no-network", action="store_true", help="skip API pings")
    a = ap.parse_args()

    print(f"Python {sys.version.split()[0]}  (needs 3.12)")
    print(f".env present: {(ROOT / '.env').exists()}   NCBI_API_KEY set: {bool(os.environ.get('NCBI_API_KEY'))}"
          f"   (without a key ClinVar is limited to 3 requests/s)")
    try:
        import redis
        redis.Redis(socket_connect_timeout=1).ping()
        print("Redis: running")
    except Exception:
        print("Redis: not running (fine: SQLite cache is used)")

    variants, rep = parse_tsv(a.input)
    ok = [v for v in variants if v.parse_status == ParseStatus.OK]
    ok, _ = apply_genome_build(ok)
    build = ok[0].genome_build
    print(f"\nInput {a.input}: {rep.total_rows} rows, {len(ok)} unique variants, build {build.value}")

    con = sqlite3.connect(a.db_path)

    def have(table: str) -> set[str]:
        try:
            return {r[0] for r in con.execute(f"SELECT variant_key FROM {table} WHERE status='AVAILABLE'")}
        except sqlite3.Error:
            return set()

    cv, gn, mv, tc = have("clinvar_evidence"), have("gnomad_evidence"), have("myvariant_evidence"), have("tcga_evidence")
    need = {"clinvar": 0, "myvariant": 0, "gnomad": 0, "tcga": 0}
    seen_tc: set[str] = set()
    for v in ok:
        if _variant_key(v) and _variant_key(v) not in cv:
            need["clinvar"] += 0  # NCBI position search is cache-only by default; exact-allele ClinVar comes via MyVariant
        if build in _BUILD_TO_ASSEMBLY:
            h, _ = resolve_hgvs_genomic_id(v)
            if h and f"{_BUILD_TO_ASSEMBLY[build]}:{h}" not in mv:
                need["myvariant"] += 1
        if build in _BUILD_TO_DATASET:
            g = build_gnomad_variant_id(v)
            if g and f"{_BUILD_TO_DATASET[build]}:{g}" not in gn:
                need["gnomad"] += 1
        aa = parse_protein_change(v.hgvsp)
        if v.gene and aa and f"{v.gene}:{aa}" not in tc and f"{v.gene}:{aa}" not in seen_tc:
            seen_tc.add(f"{v.gene}:{aa}")
            need["tcga"] += 1

    lat: dict[str, float] = {}
    print("\nsource      live lookups needed   API status")
    for n in need:
        if a.no_network:
            print(f"{n:<10}  {need[n]:>10}")
            continue
        ok_, t, msg = ping(n)
        lat[n] = t
        print(f"{n:<10}  {need[n]:>10}            {'OK ' if ok_ else 'FAIL'} {t:.2f}s  {msg}")

    if lat:
        workers = {"clinvar": 3, "myvariant": 4, "gnomad": 1, "tcga": 4}
        per_call = {"clinvar": 2 * max(lat.get("clinvar", 1), 0.35), "myvariant": lat.get("myvariant", 1) / 50 + 0.05,
                    "gnomad": max(lat.get("gnomad", 1), 1.0), "tcga": lat.get("tcga", 1) * 2}
        est = {n: need[n] * per_call[n] / (workers[n] if n != "clinvar" else 1) for n in need}
        slowest = max(est, key=est.get)
        print(f"\nESTIMATE (rough, from ping latency): slowest source is {slowest}, about {est[slowest] / 60:.1f} min "
              f"of live lookups; sources run in parallel. Cached lookups are instant.")
        print("Run with --offline for a cache-only run, or GNOMAD_ENABLED=0 / TCGA_ENABLED=0 to skip a slow source.")


if __name__ == "__main__":
    main()
