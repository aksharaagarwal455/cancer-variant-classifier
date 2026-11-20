import pytest

from app.database import db


@pytest.fixture()
def tmp_db(tmp_path, monkeypatch):
    """Fresh SQLite cache per test; Redis disabled; no env flags leaking between tests."""
    db._engine = None
    db._redis_client = False
    for k in ("EVIDENCE_OFFLINE", "ACMG_IGNORE_CLINVAR", "ACMG_SPLICE_PVS1", "TCGA_ENABLED", "GNOMAD_ENABLED"):
        monkeypatch.delenv(k, raising=False)
    path = str(tmp_path / "test.db")
    db.get_engine(path)
    yield path
    db._engine = None
