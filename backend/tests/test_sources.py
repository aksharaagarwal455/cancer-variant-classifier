import json
import time

import httpx

from app.models.variant import GenomeBuild, Variant
from app.services.sources import base
from app.services.sources.base import CircuitBreaker, SourceStats, cached_lookup, request_with_retry, run_parallel
from app.services.sources.myvariant import MyVariantClient
from app.services.sources.tcga import parse_protein_change
from app.database.db import get_local


def client_for(handler):
    return httpx.Client(transport=httpx.MockTransport(handler))


def test_retry_is_bounded_and_fast(monkeypatch):
    sleeps = []
    monkeypatch.setattr(base.time, "sleep", lambda s: sleeps.append(s))
    calls = []

    def handler(req):
        calls.append(1)
        return httpx.Response(429, headers={"Retry-After": "600"})

    with client_for(handler) as c:
        r = request_with_retry(c, "GET", "http://x/", max_retries=2)
    assert r.status_code == 429
    assert len(calls) == 3
    assert max(sleeps) <= 4.0


def test_failures_are_not_cached_but_successes_are(tmp_db):
    stats = SourceStats()
    fail = cached_lookup(source="s", table="clinvar_evidence", key="k1", db_path=tmp_db,
                         fetch=lambda: ("SOURCE_UNAVAILABLE", {"error": "HTTP 429"}), breaker=None, stats=stats)
    assert fail.status.value == "SOURCE_UNAVAILABLE" and get_local("clinvar_evidence", "k1", db_path=tmp_db) is None
    ok = cached_lookup(source="s", table="clinvar_evidence", key="k2", db_path=tmp_db,
                       fetch=lambda: ("AVAILABLE", {"found": False}), breaker=None, stats=stats)
    assert ok.source_type == "live_api"
    again = cached_lookup(source="s", table="clinvar_evidence", key="k2", db_path=tmp_db,
                          fetch=lambda: (_ for _ in ()).throw(AssertionError("must not hit network")), breaker=None, stats=stats)
    assert again.source_type == "local_db"


def test_circuit_breaker_stops_calling_a_failing_source(tmp_db):
    breaker = CircuitBreaker(threshold=3)
    stats = SourceStats()
    calls = []

    def fetch():
        calls.append(1)
        return "SOURCE_UNAVAILABLE", {"error": "HTTP 429"}

    results = [cached_lookup(source="g", table="gnomad_evidence", key=f"k{i}", db_path=tmp_db, fetch=fetch,
                             breaker=breaker, stats=stats) for i in range(10)]
    assert len(calls) == 3
    assert [r.source_type for r in results].count("skipped") == 7


def test_offline_mode_never_calls_network(tmp_db, monkeypatch):
    monkeypatch.setenv("EVIDENCE_OFFLINE", "1")
    r = cached_lookup(source="g", table="gnomad_evidence", key="k", db_path=tmp_db,
                      fetch=lambda: (_ for _ in ()).throw(AssertionError("network used")), breaker=None, stats=SourceStats())
    assert r.source_type == "skipped"


def test_deadline_skips_remaining_lookups(tmp_db):
    breaker = CircuitBreaker(threshold=99, deadline_seconds=0.05)
    time.sleep(0.1)
    r = cached_lookup(source="g", table="gnomad_evidence", key="k", db_path=tmp_db,
                      fetch=lambda: (_ for _ in ()).throw(AssertionError("network used")), breaker=breaker, stats=SourceStats())
    assert r.source_type == "skipped"


def test_run_parallel_is_faster_than_sequential_and_keeps_order():
    t0 = time.monotonic()
    out = run_parallel(list(range(8)), lambda i: (time.sleep(0.1), i)[1], workers=4)
    assert out == list(range(8))
    assert time.monotonic() - t0 < 0.6


def variants(n):
    return [Variant(variant_id=f"chr1:{100 + i}:A:T", original_row_id=i, chromosome="1", position=100 + i, ref="A",
                    alt="T", genome_build=GenomeBuild.GRCH37, gene="G") for i in range(n)]


def test_myvariant_uses_one_batch_request_for_many_variants(tmp_db):
    seen = []

    def handler(req):
        body = req.content.decode()
        seen.append(body)
        ids = [x for x in httpx.QueryParams(body)["ids"].split(",")]
        return httpx.Response(200, json=[{"query": i, "cadd": {"phred": 1}} for i in ids])

    mv = MyVariantClient(db_path=tmp_db)
    orig = httpx.Client

    def fake_client(*a, **k):
        k.pop("timeout", None)
        return orig(transport=httpx.MockTransport(handler))

    import app.services.sources.myvariant as mod
    mod.httpx.Client = fake_client
    try:
        res = mv.annotate_batch(variants(20))
    finally:
        mod.httpx.Client = orig
    assert len(seen) == 1
    assert all(r.status.value == "AVAILABLE" for r in res.values())
    res2 = mv.annotate_batch(variants(20))
    assert all(r.source_type == "local_db" for r in res2.values())


def test_tcga_protein_change_parsing():
    assert parse_protein_change("p.V600E") == "V600E"
    assert parse_protein_change("p.Val600Glu") == "V600E"
    assert parse_protein_change("p.E210fs") is None
    assert parse_protein_change(None) is None
