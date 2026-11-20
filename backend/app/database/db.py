"""
SQLite local-first evidence database.

Design: every evidence lookup checks this DB FIRST. Only on a genuine
cache miss does the caller fall through to a live API -- and when it
does, the live result is written back here before being returned, so
the next lookup for the same variant is served locally.

This is intentionally a real, inspectable SQLite file (not an opaque
JSON blob) so it can be queried directly for auditing:
  sqlite3 data/variant_classifier.db "SELECT * FROM clinvar_evidence LIMIT 5;"
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from sqlalchemy import create_engine, event, text
from sqlalchemy.engine import Engine
import redis as redis_lib

_redis_client = None

def get_redis():
    global _redis_client
    if _redis_client is None:
        try:
            _redis_client = redis_lib.Redis(host="localhost", port=6379, db=0, decode_responses=True,
                                            socket_connect_timeout=1, socket_timeout=1)
            _redis_client.ping()
        except Exception:
            _redis_client = False
    return _redis_client

DEFAULT_DB_PATH = "data/variant_classifier.db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS clinvar_evidence (
    variant_key TEXT PRIMARY KEY,
    status TEXT NOT NULL,
    data_json TEXT NOT NULL,
    retrieved_at TEXT NOT NULL,
    source_type TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS gnomad_evidence (
    variant_key TEXT PRIMARY KEY,
    status TEXT NOT NULL,
    data_json TEXT NOT NULL,
    retrieved_at TEXT NOT NULL,
    source_type TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS myvariant_evidence (
    variant_key TEXT PRIMARY KEY,
    status TEXT NOT NULL,
    data_json TEXT NOT NULL,
    retrieved_at TEXT NOT NULL,
    source_type TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS tcga_evidence (
    variant_key TEXT PRIMARY KEY,
    status TEXT NOT NULL,
    data_json TEXT NOT NULL,
    retrieved_at TEXT NOT NULL,
    source_type TEXT NOT NULL
);
"""

_engine: Engine | None = None


def get_engine(db_path: str = DEFAULT_DB_PATH) -> Engine:
    global _engine
    if _engine is None:
        Path(db_path).parent.mkdir(parents=True, exist_ok=True)
        _engine = create_engine(f"sqlite:///{db_path}", future=True, connect_args={"timeout": 30})

        @event.listens_for(_engine, "connect")
        def _pragmas(dbapi_conn, _):
            cur = dbapi_conn.cursor()
            cur.execute("PRAGMA journal_mode=WAL")
            cur.execute("PRAGMA synchronous=NORMAL")
            cur.close()

        with _engine.begin() as conn:
            for stmt in SCHEMA.strip().split(";"):
                if stmt.strip():
                    conn.execute(text(stmt))
    return _engine


def _stale_failure(status, retrieved_at, ttl_seconds: int = 30) -> bool:
    """A SOURCE_UNAVAILABLE row older than ttl is treated as a cache miss so it gets retried."""
    if status != "SOURCE_UNAVAILABLE":
        return False
    try:
        age = (datetime.now(timezone.utc) - datetime.fromisoformat(retrieved_at)).total_seconds()
    except Exception:
        return True
    return age > ttl_seconds


def get_local(table: str, variant_key: str, db_path: str = DEFAULT_DB_PATH) -> dict | None:
    """
    Three-tier lookup: Redis (fastest) -> SQLite (persistent) -> None (caller falls to live API).
    A SQLite hit repopulates Redis so the next lookup is fast too.
    """
    r = get_redis()
    if r:
        rkey = f"{table}:{variant_key}"
        cached = r.get(rkey)
        if cached:
            entry = json.loads(cached)
            if entry.get("status") == "AVAILABLE":
                entry["source_type"] = "redis"
                return entry

    engine = get_engine(db_path)
    with engine.connect() as conn:
        row = conn.execute(
            text(f"SELECT status, data_json, retrieved_at, source_type FROM {table} WHERE variant_key = :k"),
            {"k": variant_key},
        ).fetchone()
    if row is None or row[0] != "AVAILABLE":
        return None
    entry = {"status": row[0], "data": json.loads(row[1]), "retrieved_at": row[2], "source_type": row[3]}
    if r:
        r.set(f"{table}:{variant_key}", json.dumps(entry), ex=86400)
    return entry


def save_local(
    table: str,
    variant_key: str,
    status: str,
    data: dict,
    source_type: str = "live_api",
    db_path: str = DEFAULT_DB_PATH,
) -> None:
    if status != "AVAILABLE":
        return
    engine = get_engine(db_path)
    now = datetime.now(timezone.utc).isoformat()
    r = get_redis()
    if r:
        r.set(f"{table}:{variant_key}", json.dumps({"status": status, "data": data, "retrieved_at": now, "source_type": source_type}), ex=86400)
    with engine.begin() as conn:
        conn.execute(
            text(
                f"INSERT OR REPLACE INTO {table} (variant_key, status, data_json, retrieved_at, source_type) "
                f"VALUES (:k, :s, :d, :t, :src)"
            ),
            {"k": variant_key, "s": status, "d": json.dumps(data), "t": now, "src": source_type},
        )


def count_rows(table: str, db_path: str = DEFAULT_DB_PATH) -> int:
    engine = get_engine(db_path)
    with engine.connect() as conn:
        return conn.execute(text(f"SELECT COUNT(*) FROM {table}")).scalar_one()


EVIDENCE_TABLES = ("clinvar_evidence", "gnomad_evidence", "myvariant_evidence", "tcga_evidence")


def purge_failures(db_path: str = DEFAULT_DB_PATH) -> int:
    """Delete failure rows left behind by older versions. Returns how many were removed."""
    engine = get_engine(db_path)
    removed = 0
    with engine.begin() as conn:
        for t in EVIDENCE_TABLES:
            removed += conn.execute(text(f"DELETE FROM {t} WHERE status != 'AVAILABLE'")).rowcount or 0
    return removed
