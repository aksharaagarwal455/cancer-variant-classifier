"""End-to-end checks that need no internet: the GUI API, and the whole LIVE code path against a fake server."""
from __future__ import annotations

from pathlib import Path

import httpx
import pytest
from fastapi.testclient import TestClient

from app.services.classification.acmg import evaluate_acmg
from app.services.classification.amp import evaluate_amp
from app.services.evidence.evidence_service import EvidenceService
from app.services.normalization.coordinates import apply_genome_build
from app.services.parser.tsv_parser import ParseStatus, parse_tsv

ROOT = Path(__file__).resolve().parents[2]
SAMPLE = ROOT / "data/input/S11_labeled_subset.tsv"


def small_file(tmp_path, n=6):
    lines = SAMPLE.read_text(encoding="utf-8").splitlines()
    p = tmp_path / "small.tsv"
    p.write_text("\n".join(lines[: n + 1]) + "\n", encoding="utf-8")
    return p


# ---------------------------------------------------------------- GUI API
def test_api_classify_offline_returns_tiers_flags_and_timings(tmp_db, tmp_path, monkeypatch):
    from app.api import main
    monkeypatch.setattr(main, "DB_PATH", tmp_db)
    client = TestClient(main.app)
    with open(small_file(tmp_path), "rb") as fh:
        r = client.post("/api/classify?limit=5&offline=true", files={"file": ("small.tsv", fh)})
    assert r.status_code == 200
    body = r.json()
    assert body["total_variants"] == 5 and set(body["source_timings"]) == {"myvariant", "clinvar", "gnomad", "tcga"}
    first = body["results"][0]
    for key in ("acmg_classification", "amp_tier", "amp_tier_label", "evidence_completeness", "review_flags", "population_af"):
        assert key in first
    assert first["amp_tier"].startswith("TIER_")
    import os
    assert os.environ.get("EVIDENCE_OFFLINE") != "1"          # the offline switch must not leak after the request


def test_api_rejects_a_file_without_variants(tmp_db, tmp_path, monkeypatch):
    from app.api import main
    monkeypatch.setattr(main, "DB_PATH", tmp_db)
    bad = tmp_path / "bad.tsv"
    bad.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with open(bad, "rb") as fh:
        r = TestClient(main.app).post("/api/classify?limit=5", files={"file": ("bad.tsv", fh)})
    assert r.status_code == 422


def test_health_and_index(tmp_db):
    from app.api import main
    c = TestClient(main.app)
    assert c.get("/health").json() == {"status": "ok"}
    assert c.get("/").status_code == 200


def test_api_export_tsv_and_xlsx(tmp_db):
    from app.api import main
    c = TestClient(main.app)
    payload = {
        "results": [
            {
                "gene": "TP53",
                "variant_id": "chr17:7577120:C:T",
                "chromosome": "17",
                "position": 7577120,
                "ref": "C",
                "alt": "T",
                "genome_build": "GRCh37",
                "acmg_classification": "PATHOGENIC",
                "acmg_score": 11,
                "acmg_criteria": "PVS1=SUPPORTED",
                "acmg_explanation": "Score 11 >= 10.",
                "amp_tier": "TIER_I",
                "amp_tier_label": "Tier I - Strong Clinical Significance",
                "amp_explanation": "TP53 established biomarker with supporting evidence.",
                "evidence_completeness": "HIGH",
                "review_flags": [],
                "population_af": 0.00001,
                "population_af_source": "gnomad",
                "clinvar_exact_allele": ["Pathogenic"],
                "tcga_cases": 45,
            }
        ],
        "format": "tsv",
    }
    # TSV export
    r_tsv = c.post("/api/export", json=payload)
    assert r_tsv.status_code == 200
    assert "TP53" in r_tsv.text
    assert "chr17:7577120:C:T" in r_tsv.text

    # XLSX export
    payload["format"] = "xlsx"
    r_xlsx = c.post("/api/export", json=payload)
    assert r_xlsx.status_code == 200
    assert len(r_xlsx.content) > 1000
    assert r_xlsx.headers["content-type"] == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"


# ---------------------------------------------------------------- live path against a fake server
def test_live_path_end_to_end_with_fake_servers(tmp_db, tmp_path, monkeypatch):
    calls = {"myvariant": 0, "gnomad": 0, "gdc": 0, "eutils": 0}

    def handler(request: httpx.Request) -> httpx.Response:
        host = request.url.host
        if "myvariant" in host:
            calls["myvariant"] += 1
            ids = httpx.QueryParams(request.content.decode())["ids"].split(",")
            out = []
            for i, vid in enumerate(ids):
                if i == 0:   # first variant: exact-allele ClinVar pathogenic
                    out.append({"query": vid, "clinvar": {"variant_id": 1, "rcv": [
                        {"accession": "RCV1", "clinical_significance": "Pathogenic", "conditions": {"name": "X"}}]}})
                else:
                    out.append({"query": vid, "notfound": True})
            return httpx.Response(200, json=out)
        if "gnomad" in host:
            calls["gnomad"] += 1
            return httpx.Response(200, json={"errors": [{"message": "Variant not found"}], "data": {"variant": None}})
        if "gdc" in host:
            calls["gdc"] += 1
            return httpx.Response(200, json={"data": {"pagination": {"total": 3}}})
        calls["eutils"] += 1
        return httpx.Response(500)

    real = httpx.Client
    monkeypatch.setattr(httpx, "Client", lambda *a, **k: real(transport=httpx.MockTransport(handler)))
    monkeypatch.setenv("GNOMAD_MIN_INTERVAL", "0")

    variants, _ = parse_tsv(str(small_file(tmp_path, 6)))
    ok = [v for v in variants if v.parse_status == ParseStatus.OK]
    ok, _ = apply_genome_build(ok, forced_build="hg19")
    svc = EvidenceService(db_path=tmp_db)
    ev = svc.gather(ok)

    assert calls["eutils"] == 0, "the NCBI position search must stay off by default"
    assert calls["myvariant"] >= 1 and calls["gnomad"] >= 1
    assert not svc.errors
    first = ev[ok[0].variant_id]
    assert first.clinvar.data["match"] == "exact_allele"
    assert first.clinvar.data["records"][0]["clinical_significance"] == "Pathogenic"
    for v in ok:   # both engines run on live-path evidence without errors
        assert evaluate_acmg(v, ev[v.variant_id]).classification
        assert evaluate_amp(v, ev[v.variant_id]).tier

    before = dict(calls)
    svc2 = EvidenceService(db_path=tmp_db)
    svc2.gather(ok)                                   # second run: everything from the cache, zero network calls
    assert calls == before
    assert all(t["live_calls"] == 0 for t in svc2.timings.values())


def test_api_warns_when_genome_build_cannot_be_detected(tmp_db, tmp_path, monkeypatch):
    from app.api import main
    monkeypatch.setattr(main, "DB_PATH", tmp_db)
    client = TestClient(main.app)
    with open(small_file(tmp_path, 3), "rb") as fh:          # too few landmark genes to detect the build
        body = client.post("/api/classify?limit=3&offline=true", files={"file": ("s.tsv", fh)}).json()
    assert "UNKNOWN" in body["genome_build"] and any("Genome build" in w for w in body["warnings"])
    with open(small_file(tmp_path, 3), "rb") as fh:          # the user supplies it
        body = client.post("/api/classify?limit=3&offline=true&genome_build=hg19", files={"file": ("s.tsv", fh)}).json()
    assert body["genome_build"] == "GRCh37" and not body["warnings"]


def test_api_rejects_unknown_genome_build(tmp_db, tmp_path, monkeypatch):
    from app.api import main
    monkeypatch.setattr(main, "DB_PATH", tmp_db)
    with open(small_file(tmp_path, 3), "rb") as fh:
        r = TestClient(main.app).post("/api/classify?genome_build=hg99", files={"file": ("s.tsv", fh)})
    assert r.status_code == 422
