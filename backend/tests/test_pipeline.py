import time
from pathlib import Path

from app.models.variant import GenomeBuild
from app.services.evidence.evidence_service import EvidenceService
from app.services.normalization.coordinates import apply_genome_build
from app.services.parser.tsv_parser import ParseStatus, parse_tsv
from app.services.sources.base import EvidenceResult, SourceStatus

ROOT = Path(__file__).resolve().parents[2]


def load(name):
    vs, rep = parse_tsv(str(ROOT / "data" / "input" / name))
    ok = [v for v in vs if v.parse_status == ParseStatus.OK]
    ok, _ = apply_genome_build(ok)
    return ok, rep


def test_labelled_subset_parses_and_is_grch37():
    ok, rep = load("S11_labeled_subset.tsv")
    assert len(ok) == 42 and ok[0].genome_build == GenomeBuild.GRCH37


def test_duplicate_rows_are_collapsed():
    ok, rep = load("hg19Somatic_validation.tsv")
    assert rep.total_rows > len(ok) and rep.duplicates == rep.total_rows - len(ok)


def test_offline_gather_never_crashes_and_reports_skips(tmp_db, monkeypatch):
    monkeypatch.setenv("EVIDENCE_OFFLINE", "1")
    ok, _ = load("S11_labeled_subset.tsv")
    svc = EvidenceService(db_path=tmp_db)
    ev = svc.gather(ok[:5])
    assert len(ev) == 5
    assert all(e.clinvar.source_type == "skipped" for e in ev.values())
    assert set(svc.timings) == {"myvariant", "clinvar", "gnomad", "tcga"}


def test_sources_run_concurrently_not_one_after_another(tmp_db):
    ok, _ = load("S11_labeled_subset.tsv")
    svc = EvidenceService(db_path=tmp_db)

    class Slow:
        def __init__(self, name):
            self.stats = type("S", (), {"summary": lambda s: {}})()
            self.name = name

        def annotate_batch(self, variants, progress=None):
            time.sleep(0.3)
            return {v.variant_id: EvidenceResult(self.name, SourceStatus.NOT_ASSESSABLE) for v in variants}

    svc.clients = {n: Slow(n) for n in ("myvariant", "clinvar", "gnomad", "tcga")}
    t0 = time.monotonic()
    svc.gather(ok[:3])
    assert time.monotonic() - t0 < 0.9   # 4 x 0.3 s sequentially would be 1.2 s


def test_a_crashing_source_does_not_break_the_run(tmp_db):
    ok, _ = load("S11_labeled_subset.tsv")
    svc = EvidenceService(db_path=tmp_db)

    class Boom:
        stats = type("S", (), {"summary": lambda s: {}})()

        def annotate_batch(self, variants, progress=None):
            raise RuntimeError("boom")

    svc.clients["tcga"] = Boom()
    svc.clients["myvariant"] = Boom()
    svc.clients["clinvar"] = Boom()
    svc.clients["gnomad"] = Boom()
    ev = svc.gather(ok[:2])
    assert len(ev) == 2 and all(e.tcga is None for e in ev.values())


def test_rules_folder_matches_the_code():
    import json
    from export_rules import current_rules   # scripts/ is on the pytest path

    for name, content in current_rules().items():
        on_disk = json.loads((ROOT / "rules" / name).read_text())
        assert on_disk == content, f"rules/{name} is out of date: run python scripts/export_rules.py"


def test_wilson_interval_known_value():
    from make_report import wilson

    lo, hi = wilson(6, 6)
    assert abs(lo - 0.6097) < 0.001 and hi == 1.0
