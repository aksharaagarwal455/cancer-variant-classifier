# Architecture

```
TSV input
  -> parser            services/parser/tsv_parser.py        duplicate headers, "." nulls, dedup on chrom:pos:ref:alt
  -> genome build      services/normalization/coordinates.py  landmark-gene coordinates (TP53, BRCA1/2, KRAS, EGFR)
  -> evidence          services/evidence/evidence_service.py  4 sources IN PARALLEL, DB-first
        myvariant.py   MyVariant.info   batched POSTs (50 ids per request)
        clinvar.py     NCBI E-utilities by POSITION (informational only, cache-only unless CLINVAR_NCBI_ENABLED=1)
        gnomad.py      gnomAD GraphQL    2 workers, calls spaced >= 1 s, circuit breaker
        tcga.py        NCI GDC           unique (gene, protein change) pairs, 4 workers
  -> ACMG engine       services/classification/acmg.py       8 criteria, Tavtigian points
  -> AMP engine        services/classification/amp.py        Tier I-IV, independent of ACMG
  -> CLI (TSV/JSON/summary JSON)  or  FastAPI + browser GUI
```

## The cache (DB-first)

Every lookup checks Redis (if running), then SQLite (`data/variant_classifier.db`), and only then the network.
Only **successful** results are stored. Failures (timeouts, HTTP 429, outages) are never cached, so a bad moment
cannot poison later runs, and they are never retried again inside the same run.
Tables: `clinvar_evidence`, `gnomad_evidence`, `myvariant_evidence`, `tcga_evidence`.

## Safeguards that keep a run from hanging (`services/sources/base.py`)

| Safeguard | Behaviour |
|---|---|
| Timeouts | 5 s connect, 8-20 s read (per source) |
| Retries | at most 2 fast retries, sleeps capped at 4 s (a huge `Retry-After` is ignored beyond that) |
| Circuit breaker | after 3-5 consecutive failures the source is skipped for the rest of the run |
| Wall-clock limit | `--max-seconds` (default 900): after that, remaining cache misses are skipped |
| Offline mode | `--offline` or `EVIDENCE_OFFLINE=1`: cache only, zero network calls |
| Isolation | a crashing source returns nothing for that source; the run still completes |
| Honest degradation | skipped lookups are reported as `SOURCE_UNAVAILABLE`; the engines then use the allele-frequency and annotation columns already in the input file. Nothing is invented. |

## Parallelism

The four sources run at the same time, so wall time is about the slowest source, not the sum.
Inside a source: MyVariant sends 50 ids per request (4 workers); gnomAD uses 2 workers behind a shared limiter that spaces
requests at least 1 s apart; TCGA uses 4 workers on unique (gene, protein change) pairs; the NCBI position search is cache-only
unless `CLINVAR_NCBI_ENABLED=1`. Tune gnomAD spacing with `GNOMAD_MIN_INTERVAL` (seconds).
Skip a source entirely with `GNOMAD_ENABLED=0` or `TCGA_ENABLED=0`.

## Configuration

`.env` in the project root is loaded automatically (`backend/app/config.py`). Put `NCBI_API_KEY` there for faster ClinVar.
Earlier versions never read `.env`, so a key placed there was ignored.

## Outputs

`<name>_results.tsv/.json` (one row per unique variant, with ACMG class, score, criteria, AMP tier, source badges)
and `<name>_results_summary.json` (timings per source, lookup counts, class counts, flags used).

## API (FastAPI)

| Route | Purpose |
|---|---|
| `GET /` | browser GUI |
| `GET /health` | liveness |
| `POST /api/classify` (multipart `file`, query `limit`) | classify an uploaded TSV; returns results plus `evidence_seconds` and per-source `evidence_timings` |

The endpoint is a plain `def`, so FastAPI runs it in a worker thread and the server stays responsive during a long run.


## ClinVar evidence: exact allele, not position

The scoring engines read `evidence.clinvar`, which is built from the MyVariant result (MyVariant matches on the full genomic
HGVS id: chromosome, position, ref and alt). The NCBI position search is kept as `evidence.clinvar_position` and only raises a
review flag ("ClinVar has an entry at this position but not for this exact allele"). On the hg19 set the position search
returned pathogenic records for 55 variants that the file's own exact-allele CLNSIG labels VUS or leaves empty, which is why it is no
longer scored. Full explanation: `docs/LOGIC.md`.
