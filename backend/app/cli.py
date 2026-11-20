"""
CLI for the DB-first variant classifier.

Interactive:
  PYTHONPATH=backend python -m app.cli

Non-interactive:
  PYTHONPATH=backend python -m app.cli --input data/input/sample_S11.tsv --limit 30
"""

from __future__ import annotations

import csv
import json
import os
import time
from collections import Counter
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.progress import track
from rich.table import Table

from app import config

from app.services.classification.acmg import evaluate_acmg
from app.services.classification.amp import evaluate_amp
from app.services.classification.evidence_utils import population_af
from app.services.classification.quality import TIER_LABELS, assess_quality
from app.services.evidence.evidence_service import EvidenceService
from app.services.normalization.coordinates import apply_genome_build
from app.services.parser.tsv_parser import ParseStatus, parse_tsv

app = typer.Typer()
console = Console()

DISCLAIMER = (
    "Research/evidence-aggregation tool, not a clinical diagnostic device. "
    "ACMG uses a documented subset of criteria (PVS1, PM2, PP3, PP5, BA1, BS1, BP4, BP6). "
    "AMP tiering is evaluated independently, not copied from ACMG. "
    "Requires review by a qualified clinical genomics professional before clinical use."
)


def _clean_path(raw: str) -> str:
    p = raw.strip().strip("'").strip('"')
    p = p.replace("\\ ", " ")
    return str(Path(p).expanduser())


def _count_rows(path: Path) -> int:
    try:
        with open(path, "rb") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except OSError:
        return 0


def _choose_input_file() -> str:
    input_dir = Path("data/input") if Path("data/input").exists() else Path("../data/input")
    available = sorted(input_dir.glob("*.tsv"))
    tbl = Table(title="Choose an input file", show_lines=False)
    tbl.add_column("#", justify="right")
    tbl.add_column("file")
    tbl.add_column("variants", justify="right")
    for i, f in enumerate(available, 1):
        tbl.add_row(str(i), f.name, str(_count_rows(f)))
    console.print(tbl)
    console.print("  Type the number, a file name, or drag a file here / paste a full path.")

    while True:
        answer = typer.prompt("Which file do you want to classify").strip()
        if answer.isdigit() and 1 <= int(answer) <= len(available):
            return str(available[int(answer) - 1])
        candidate = _clean_path(answer)
        if Path(candidate).is_file():
            return candidate
        if (input_dir / answer).is_file():
            return str(input_dir / answer)
        console.print(f"[red]Could not find '{answer}'. Try again.[/red]")


REQUIRED_COLUMNS = {"chromosome": ["CHROM", "chrom", "chr", "chromosome"]}
HELPFUL_COLUMNS = ["Gene", "ExonicFunc.refGene", "Func.refGene", "AF", "REVEL_score", "CADD_phred",
                   "SIFT_pred", "Polyphen2_HVAR_pred", "MetaSVM_pred", "CLNSIG"]


def _input_readiness(path: str) -> None:
    """Tell the user up front which annotation columns exist, because missing ones weaken the result."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    have = set(header)
    missing = [c for c in HELPFUL_COLUMNS if c not in have]
    if missing:
        console.print(f"  [yellow]Input check: {len(HELPFUL_COLUMNS) - len(missing)}/{len(HELPFUL_COLUMNS)} expected "
                      f"annotation columns present. Missing: {', '.join(missing)}. "
                      f"Results will rely more on live sources and be flagged LOW/MEDIUM completeness.[/yellow]")
    else:
        console.print(f"  [green]Input check: all {len(HELPFUL_COLUMNS)} expected annotation columns present.[/green]")


@app.command()
def classify(
    input: Optional[str] = typer.Option(None, "--input", help="Path to input TSV. If omitted, you will be asked to choose."),
    output: Optional[str] = typer.Option(None, "--output", help="Path to output TSV. If omitted, derived from the input filename."),
    limit: Optional[int] = typer.Option(None, "--limit", help="Max variants to classify (0 = all). If omitted, you will be asked."),
    db_path: str = typer.Option("data/variant_classifier.db", "--db-path"),
    genome_build: Optional[str] = typer.Option(None, "--genome-build", help="hg19 or hg38. Overrides auto-detection."),
    offline: bool = typer.Option(False, "--offline", help="Cache only: zero network calls (fast, reproducible)."),
    max_seconds: int = typer.Option(900, "--max-seconds", help="Wall-clock limit for live lookups per source; after it remaining lookups are skipped."),
    show: Optional[int] = typer.Option(None, "--show", help="Print the first N variants with their tiers (default 15 when interactive, 0 otherwise; -1 = all)."),
):
    if offline:
        os.environ["EVIDENCE_OFFLINE"] = "1"
    is_offline = os.environ.get("EVIDENCE_OFFLINE", "0") == "1"
    interactive = input is None
    if interactive:
        input = _choose_input_file()
    else:
        input = _clean_path(input)
        if not Path(input).is_file():
            console.print(f"[red]File not found: {input}[/red]")
            raise typer.Exit(code=1)

    if limit is None:
        if interactive:
            limit = typer.prompt("How many variants to classify? (0 = all)", default=0, type=int)
        else:
            limit = 0
    if interactive and not is_offline:
        if typer.confirm("Cache only (offline, instant, no network)?", default=False):
            os.environ["EVIDENCE_OFFLINE"] = "1"
            is_offline = True
    if show is None:
        show = 15 if interactive else 0

    if output is None:
        stem = Path(input).stem
        out_dir = Path("data/output") if Path("data/output").exists() else Path("../data/output")
        output = str(out_dir / f"{stem}_results.tsv")

    console.print(f"[bold]Parsing[/bold] {input} ...")
    _input_readiness(input)
    variants, parse_report = parse_tsv(input)
    console.print(f"  total_rows={parse_report.total_rows} parsed_ok={parse_report.parsed_ok} "
                  f"invalid={parse_report.invalid} duplicates={parse_report.duplicates}")

    ok_variants = [v for v in variants if v.parse_status == ParseStatus.OK]
    if not ok_variants:
        console.print("[red]No valid variants found in this file. Check the column names (CHROM, POS, REF, ALT, Gene).[/red]")
        raise typer.Exit(code=1)

    console.print("[bold]Detecting genome build[/bold] (whole file, before --limit) ...")
    ok_variants, build_evidence = apply_genome_build(ok_variants, forced_build=genome_build)
    console.print(f"  build: {ok_variants[0].genome_build.value}  evidence: {build_evidence}")
    if ok_variants[0].genome_build.name == "UNKNOWN" and interactive:
        answer = typer.prompt("Genome build could not be detected. Type hg19 or hg38 (Enter = skip MyVariant/gnomAD)", default="").strip()
        if answer:
            try:
                ok_variants, build_evidence = apply_genome_build(ok_variants, forced_build=answer)
            except ValueError as e:
                console.print(f"[red]{e}[/red]")
                raise typer.Exit(code=1)
            console.print(f"  build set by you: {ok_variants[0].genome_build.value}")
    if ok_variants[0].genome_build.name == "UNKNOWN":
        console.print("[red]  WARNING: genome build not determined; MyVariant lookups are disabled. "
                      "Re-run with --genome-build hg19 or hg38.[/red]")

    if limit > 0:
        ok_variants = ok_variants[:limit]
        console.print(f"  [yellow]--limit applied: {limit} variants[/yellow]")
    elif len(ok_variants) > 100:
        console.print(f"  [yellow]Classifying all {len(ok_variants)} variants. Variants not already in the cache "
                      f"call live APIs, so a first run can take several minutes.[/yellow]")

    mode = "offline (cache only, no network)" if is_offline else f"live (cache first; each source stops after {max_seconds}s)"
    console.print(f"[bold]Gathering evidence[/bold] -- mode: {mode}")
    svc = EvidenceService(db_path=db_path)
    if svc.purged:
        console.print(f"  [dim]removed {svc.purged} stale failure rows left by an older version[/dim]")
    t_ev = time.monotonic()
    evidence_map = svc.gather(ok_variants, show_progress=True, max_seconds=max_seconds)
    evidence_seconds = time.monotonic() - t_ev
    for name, err in svc.errors.items():
        console.print(f"  [red]source '{name}' crashed and was skipped: {err}[/red]")

    def _count(pred):
        return sum(1 for e in evidence_map.values()
                   for r in (e.clinvar_position, e.gnomad, e.myvariant, e.tcga) if r and pred(r))

    local_hits = _count(lambda r: r.status.value != "NOT_ASSESSABLE" and r.source_type == "local_db")
    live_hits = _count(lambda r: r.status.value != "NOT_ASSESSABLE" and r.source_type == "live_api")
    skipped_n = _count(lambda r: r.source_type == "skipped")
    not_assessable = _count(lambda r: r.status.value == "NOT_ASSESSABLE")
    console.print(f"  evidence lookups: {local_hits} from local DB, {live_hits} from live API, "
                  f"{skipped_n} skipped, {not_assessable} not assessable (no lookup attempted)")
    per_variant = evidence_seconds / max(1, len(ok_variants))
    console.print(f"  evidence time: {evidence_seconds:.2f}s total = {per_variant:.4f}s per variant")
    tbl = Table(title="Per-source timing", show_lines=False)
    for col in ("source", "cache", "live ok", "live failed", "skipped", "avg live s", "max live s", "wall s"):
        tbl.add_column(col, justify="right" if col != "source" else "left")
    for name, t in svc.timings.items():
        tbl.add_row(name, str(t.get("cache_hits", 0)), str(t.get("live_ok", 0)), str(t.get("live_failed", 0)),
                    str(t.get("skipped", 0)), str(t.get("live_avg_s") or "-"), str(t.get("live_max_s") or "-"),
                    str(t.get("wall_s", "-")))
    console.print(tbl)

    console.print("[bold]Running ACMG + AMP engines[/bold] ...")
    rows = []
    acmg_counts: Counter[str] = Counter()
    amp_counts: Counter[str] = Counter()
    completeness_counts: Counter[str] = Counter()

    for v in track(ok_variants, description="Classifying"):
        ev = evidence_map[v.variant_id]
        acmg_result = evaluate_acmg(v, ev)
        amp_result = evaluate_amp(v, ev)
        acmg_counts[acmg_result.classification.value] += 1
        amp_counts[amp_result.tier.value] += 1

        row = dict(v.original_input)
        row["variant_id"] = v.variant_id
        row["genome_build"] = v.genome_build.value
        row["acmg_classification"] = acmg_result.classification.value
        row["acmg_score"] = acmg_result.score
        row["acmg_criteria"] = "; ".join(f"{c.code}={c.status}" for c in acmg_result.criteria if c.status == "SUPPORTED")
        row["acmg_explanation"] = acmg_result.explanation
        row["amp_tier"] = amp_result.tier.value
        row["amp_tier_label"] = TIER_LABELS[amp_result.tier.value]
        row["amp_explanation"] = amp_result.explanation
        pvs1 = any(cr.code == "PVS1" and cr.status == "SUPPORTED" for cr in acmg_result.criteria)
        completeness, flags = assess_quality(ev, acmg_result.classification.value, amp_result.tier.value,
                                             pvs1=pvs1, gene=v.gene)
        af, af_src = population_af(ev)
        cv_recs = (ev.clinvar.data.get("records", []) if ev.clinvar and ev.clinvar.status.value == "AVAILABLE" else [])
        row["evidence_completeness"] = completeness
        row["review_flags"] = " | ".join(flags)
        row["population_af"] = "" if af is None else f"{af:.6g}"
        row["population_af_source"] = af_src
        row["clinvar_exact_allele"] = "; ".join(sorted({r.get("clinical_significance") for r in cv_recs if r.get("clinical_significance")}))
        row["myvariant_source"] = ev.myvariant.source_type if ev.myvariant else "N/A"
        row["clinvar_source"] = ev.clinvar.source_type if ev.clinvar else "N/A"
        row["gnomad_source"] = ev.gnomad.source_type if ev.gnomad else "N/A"
        row["tcga_source"] = ev.tcga.source_type if ev.tcga else "N/A"
        row["tcga_cases"] = (ev.tcga.data or {}).get("cases", "") if ev.tcga and ev.tcga.status.value == "AVAILABLE" else ""
        rows.append(row)
        completeness_counts[completeness] += 1

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    json_path = output_path.with_suffix(".json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2, default=str)

    summary = {
        "input": input, "output": str(output_path), "genome_build": ok_variants[0].genome_build.value,
        "total_rows": parse_report.total_rows, "unique_variants_classified": len(rows),
        "mode": "offline" if is_offline else "live", "max_seconds": max_seconds,
        "env": {k: os.environ.get(k, "") for k in ("ACMG_IGNORE_CLINVAR", "ACMG_SPLICE_PVS1", "GNOMAD_ENABLED",
                                                   "TCGA_ENABLED", "EVIDENCE_OFFLINE")},
        "evidence_seconds": round(evidence_seconds, 2), "seconds_per_variant": round(per_variant, 4),
        "lookups": {"local_db": local_hits, "live_api": live_hits, "skipped": skipped_n, "not_assessable": not_assessable},
        "source_timings": svc.timings, "acmg_counts": dict(acmg_counts), "amp_counts": dict(amp_counts),
        "evidence_completeness": dict(completeness_counts),
    }
    summary_path = output_path.with_name(output_path.stem + "_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    console.print()
    console.print("=" * 60)
    console.print("[bold]VARIANT CLASSIFICATION COMPLETE[/bold]")
    console.print("=" * 60)
    console.print(f"Input: {input}")
    console.print(f"Total parsed: {parse_report.total_rows}  OK: {parse_report.parsed_ok}  "
                  f"Invalid: {parse_report.invalid}  Duplicates: {parse_report.duplicates}")
    console.print(f"Classified: {len(rows)}")
    console.print(f"Evidence source: {local_hits} local_db, {live_hits} live_api")
    console.print()
    console.print("[bold]ACMG:[/bold]")
    for k, cnt in acmg_counts.items():
        console.print(f"  {k}: {cnt}")
    console.print()
    console.print("[bold]AMP:[/bold]")
    for k, cnt in amp_counts.items():
        console.print(f"  {k}: {cnt}")
    console.print()
    tier_tbl = Table(title="Variant tiers (AMP/ASCO/CAP-style) and ACMG class", show_lines=False)
    tier_tbl.add_column("AMP tier")
    tier_tbl.add_column("meaning")
    tier_tbl.add_column("variants", justify="right")
    tier_tbl.add_column("%", justify="right")
    for t in ("TIER_I", "TIER_II", "TIER_III", "TIER_IV"):
        n = amp_counts.get(t, 0)
        tier_tbl.add_row(t.replace("_", " "), TIER_LABELS[t].split(" - ")[1].split(" (")[0], str(n),
                         f"{100 * n / max(1, len(rows)):.1f}")
    console.print(tier_tbl)
    console.print("Evidence completeness: " + ", ".join(f"{k}={completeness_counts.get(k, 0)}" for k in ("HIGH", "MEDIUM", "LOW")))
    n_show = len(rows) if show is not None and show < 0 else min(show or 0, len(rows))
    if n_show:
        vt = Table(title=f"First {n_show} classified variants", show_lines=False)
        for col in ("Gene", "variant", "ACMG", "AMP tier", "completeness"):
            vt.add_column(col)
        for r in rows[:n_show]:
            vt.add_row(str(r.get("Gene", "")), r["variant_id"][:42], r["acmg_classification"],
                       r["amp_tier"].replace("_", " "), r["evidence_completeness"])
        console.print(vt)
    console.print(f"Outputs: {output_path}, {json_path}, {summary_path}")
    console.print("=" * 60)
    console.print(f"[dim]{DISCLAIMER}[/dim]")


if __name__ == "__main__":
    app()
