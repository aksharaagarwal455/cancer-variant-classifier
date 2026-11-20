# Performance: why runs took hours, and what changed

## What was measured here

Cache-only runs (`--offline`) in this environment:

| Run | Unique variants | Evidence time | Per variant |
|---|---|---|---|
| S11 labelled subset | 42 | 0.14 s | 0.0033 s |
| hg19 somatic validation | 473 | 0.64 s | 0.0014 s |

Live measurements on the author's MacBook (`scripts/benchmark.py`, 10 variants of `sample_S12.tsv`, cache bypassed):

| Source | Median per call | Max | Notes |
|---|---|---|---|
| MyVariant | 0.34 s | 4.30 s | one call covers up to 50 variants |
| gnomAD | 0.47 s | 1.56 s | 3 of 10 calls returned HTTP 429 (rate limit) |
| TCGA / GDC | 0.57 s average | 7.3 s | 108 new lookups succeeded in the full run |
| NCBI ClinVar (position search) | 1.34 s | 1.94 s | informational only, off by default |

Full live run of `sample_S12.tsv` (611 variants, 122 live lookups, the rest cached, `--max-seconds 300`): **15.7 s total, about
0.026 s per variant**. The build environment itself has no access to the APIs, so those numbers come from the author's machine.
Check readiness with `python scripts/preflight.py`.

## Why the old code could run for hours (derived from its code and your cache; an upper bound)

| Cause | Detail |
|---|---|
| gnomAD HTTP 429 retried every run | 343 cached failures expired after 30 s and were retried each run: 4 requests plus 1+2+4 s backoff each |
| TCGA lookups uncached | the hg19 file needs 383 more live GDC lookups; each could take 2 x 10 s timeouts + 1 s |
| Everything sequential | four sources one after another, one variant at a time, no progress output, so a slow run looked frozen |
| Failures stored | rows with `SOURCE_UNAVAILABLE` were written to the cache (and counted as work) |
| `.env` never loaded | an `NCBI_API_KEY` in `.env` was ignored; ClinVar stayed at 3 requests/s |
| GUI blocking | the classify endpoint was `async` but made blocking calls, freezing the server |

Worst case for the hg19 set from the first two rows: about 16 min (gnomAD) + 134 min (TCGA) = 2.5 h sequentially, consistent
with a run that "never finished" after 3 hours.

## What changed

* Sources run in parallel; MyVariant sends 50 ids per request; TCGA uses 4 workers; gnomAD is spaced at 1 s.
* Timeouts 5 s connect / 8-20 s read, at most 2 fast retries, sleeps capped at 4 s.
* Failures are never cached or retried within a run; stale failure rows are purged on start.
* Circuit breaker (3-5 consecutive failures) and `--max-seconds` wall-clock limit per source.
* `--offline`, `GNOMAD_ENABLED=0`, `TCGA_ENABLED=0` for fast, reproducible runs.
* Progress bars per source, a timing table, and seconds per variant in every run (`*_summary.json`).
* `.env` is loaded.

## Expected live cost after the fix (estimates, not measurements)

* Exact-allele ClinVar comes through MyVariant (50 variants per request), so ClinVar no longer dominates the run time. The NCBI
  position search is cache-only unless `CLINVAR_NCBI_ENABLED=1` (then it is rate-limited: about 0.7 s per variant, 0.25 s with `NCBI_API_KEY`).
* gnomAD: at least 1 s per uncached variant (rate limit), so 100 uncached variants is about 2 minutes; if gnomAD keeps
  returning 429 the breaker stops it and the run continues with input-file frequencies.
* TCGA: 383 lookups with 4 workers; the time depends on GDC latency (run the benchmark). It is capped by `--max-seconds`.
* Cached variants cost almost nothing, and each live result is saved, so the second run of the same file is fast.

## Hang test (measured)

A local server that accepts connections and never answers was used as a dead source (`TcgaClient`, 2 s read timeout,
12 lookups, 4 workers). The run finished in 9.1 s: 7 lookups failed after 2 attempts each (about 4.5 s each), the circuit
breaker opened, and the other 5 were skipped. With the real 10 s timeout the same situation is bounded at roughly 45 s
instead of unbounded. Before this fix the cached run showed no limit at all.
