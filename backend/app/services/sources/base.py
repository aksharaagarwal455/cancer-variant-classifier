"""Shared building blocks for every evidence source.

Design goals (these are what stop a run from hanging for hours):
  * every HTTP call has a hard timeout (5 s connect, 8-20 s read) and at most 2 fast retries;
  * no sleep is ever longer than 4 s, even if a server sends a huge Retry-After;
  * failures are never cached, so a bad moment cannot poison later runs;
  * a circuit breaker + wall-clock deadline stops a failing / slow source quickly;
  * lookups run in a thread pool and the four sources run concurrently;
  * EVIDENCE_OFFLINE=1 means "cache only, zero network calls".
"""

from __future__ import annotations

import os
import statistics
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable

import httpx

from app.database.db import get_local, save_local

MAX_SLEEP_SECONDS = 4.0


class SourceStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    UNAVAILABLE = "SOURCE_UNAVAILABLE"
    NOT_ASSESSABLE = "NOT_ASSESSABLE"


@dataclass
class EvidenceResult:
    source: str
    status: SourceStatus
    data: dict[str, Any] = field(default_factory=dict)
    error: str | None = None
    retrieved_at: str | None = None
    source_type: str = "live_api"  # local_db | redis | live_api | skipped | n/a


def not_assessable(source: str, reason: str) -> EvidenceResult:
    return EvidenceResult(source=source, status=SourceStatus.NOT_ASSESSABLE, error=reason, source_type="n/a")


def skipped(source: str, reason: str) -> EvidenceResult:
    return EvidenceResult(source=source, status=SourceStatus.UNAVAILABLE, data={"error": reason},
                          error=reason, source_type="skipped")


def make_timeout(read: float = 15.0) -> httpx.Timeout:
    return httpx.Timeout(read, connect=5.0)


def is_offline() -> bool:
    return os.environ.get("EVIDENCE_OFFLINE", "0") == "1"


def source_enabled(name: str) -> bool:
    return os.environ.get(f"{name.upper()}_ENABLED", "1") != "0"


# ----------------------------------------------------------------------------- HTTP
def _sleep(seconds: float) -> None:
    time.sleep(min(max(seconds, 0.0), MAX_SLEEP_SECONDS))


def request_with_retry(client: httpx.Client, method: str, url: str, max_retries: int = 2,
                       base_backoff_seconds: float = 0.5, **kwargs: Any) -> httpx.Response | None:
    """Bounded retry: network errors, HTTP 429 and 5xx are retried at most `max_retries` times."""
    for attempt in range(max_retries + 1):
        last = attempt == max_retries
        wait = base_backoff_seconds * (2 ** attempt)
        try:
            resp = client.request(method, url, **kwargs)
        except httpx.HTTPError:
            if last:
                return None
            _sleep(wait)
            continue
        if resp.status_code == 429 or resp.status_code >= 500:
            if last:
                return resp
            if resp.status_code == 429:
                retry_after = resp.headers.get("Retry-After")
                try:
                    wait = float(retry_after) if retry_after else wait
                except ValueError:
                    pass
            _sleep(wait)
            continue
        return resp
    return None


class RateLimiter:
    """Thread-safe minimum spacing between requests to one host."""

    def __init__(self, min_interval: float):
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._next = 0.0

    def wait(self) -> None:
        with self._lock:
            now = time.monotonic()
            start = max(now, self._next)
            self._next = start + self.min_interval
        delay = start - now
        if delay > 0:
            time.sleep(delay)


# ----------------------------------------------------------------------------- bookkeeping
class SourceStats:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self.cache_hits = 0
        self.live_ok = 0
        self.live_failed = 0
        self.skipped = 0
        self.times: list[float] = []

    def hit(self, n: int = 1) -> None:
        with self._lock:
            self.cache_hits += n

    def skip(self, n: int = 1) -> None:
        with self._lock:
            self.skipped += n

    def live(self, ok: bool, seconds: float, n: int = 1) -> None:
        with self._lock:
            if ok:
                self.live_ok += n
            else:
                self.live_failed += n
            self.times.append(seconds)

    def summary(self) -> dict[str, Any]:
        with self._lock:
            t = sorted(self.times)
            return {
                "cache_hits": self.cache_hits, "live_ok": self.live_ok, "live_failed": self.live_failed,
                "skipped": self.skipped, "live_calls": len(t),
                "live_avg_s": round(sum(t) / len(t), 3) if t else None,
                "live_median_s": round(statistics.median(t), 3) if t else None,
                "live_max_s": round(t[-1], 3) if t else None,
            }


class CircuitBreaker:
    """Stops calling a source after `threshold` consecutive failures or after `deadline_seconds`."""

    def __init__(self, threshold: int = 5, deadline_seconds: float | None = None):
        self.threshold = threshold
        self.deadline_seconds = deadline_seconds
        self._t0 = time.monotonic()
        self._fails = 0
        self._lock = threading.Lock()

    def record(self, ok: bool) -> None:
        with self._lock:
            self._fails = 0 if ok else self._fails + 1

    def reason(self) -> str | None:
        with self._lock:
            if self._fails >= self.threshold:
                return f"circuit open after {self._fails} consecutive failures"
        if self.deadline_seconds is not None and time.monotonic() - self._t0 > self.deadline_seconds:
            return f"time limit of {self.deadline_seconds:.0f}s reached"
        return None

    def allow(self) -> bool:
        return self.reason() is None


def run_parallel(items: list, fn: Callable, workers: int = 4, progress: Callable[[int], None] | None = None) -> list:
    """Apply fn to every item in a thread pool; results keep the input order."""
    if not items:
        return []
    with ThreadPoolExecutor(max_workers=max(1, min(workers, len(items)))) as ex:
        futures = [ex.submit(fn, i) for i in items]
        if progress:
            for _ in as_completed(futures):
                progress(1)
        return [f.result() for f in futures]


# ----------------------------------------------------------------------------- cache-first lookup
def cached_lookup(*, source: str, table: str, key: str, db_path: str,
                  fetch: Callable[[], tuple[str, dict]], breaker: CircuitBreaker | None,
                  stats: SourceStats, cache_only_reason: str | None = None) -> EvidenceResult:
    """cache -> (offline? skip) -> (breaker open? skip) -> live fetch. Only successes are cached."""
    local = get_local(table, key, db_path=db_path)
    if local is not None:
        stats.hit()
        return EvidenceResult(source=source, status=SourceStatus(local["status"]), data=local["data"],
                              retrieved_at=local["retrieved_at"], source_type="local_db")
    if cache_only_reason:                      # source intentionally never calls the network
        return not_assessable(source, cache_only_reason)
    if is_offline():
        stats.skip()
        return skipped(source, "offline mode: not in cache")
    if breaker is not None and not breaker.allow():
        stats.skip()
        return skipped(source, breaker.reason() or "skipped")

    t0 = time.monotonic()
    try:
        status_str, data = fetch()
    except Exception as e:  # noqa: BLE001  (a source must never crash the run)
        status_str, data = SourceStatus.UNAVAILABLE.value, {"error": type(e).__name__}
    ok = status_str == SourceStatus.AVAILABLE.value
    stats.live(ok, time.monotonic() - t0)
    if breaker is not None:
        breaker.record(ok)
    if ok:
        save_local(table, key, status_str, data, source_type="live_api", db_path=db_path)
    from datetime import datetime, timezone
    return EvidenceResult(source=source, status=SourceStatus(status_str), data=data,
                          error=None if ok else data.get("error"),
                          retrieved_at=datetime.now(timezone.utc).isoformat(), source_type="live_api")


def run_lookups(*, source: str, table: str, groups: dict[str, list], fetch_for: Callable[[Any], tuple[str, dict]],
                db_path: str, breaker: CircuitBreaker, stats: SourceStats, workers: int,
                progress: Callable[[int], None] | None,
                cache_only_reason: str | None = None) -> dict[str, EvidenceResult]:
    """One cached lookup per unique key (groups: key -> variants sharing it). Returns variant_id -> result."""
    keys = list(groups)

    def one(key: str) -> EvidenceResult:
        v0 = groups[key][0]
        res = cached_lookup(source=source, table=table, key=key, db_path=db_path,
                            fetch=lambda: fetch_for(v0), breaker=breaker, stats=stats,
                            cache_only_reason=cache_only_reason)
        if progress:
            progress(len(groups[key]))
        return res

    outs = run_parallel(keys, one, workers)
    results: dict[str, EvidenceResult] = {}
    for key, res in zip(keys, outs):
        for v in groups[key]:
            results[v.variant_id] = res
    return results
