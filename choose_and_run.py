#!/usr/bin/env python3
"""Interactive launcher for the variant classifier.

Run from the project root (the folder that contains `backend/`, `data/`, `scripts/`):

    python3 choose_and_run.py

It lists the TSV files it finds, lets you pick by number, and runs the same
commands you use in the demo:
  1. classify   -> python -m app.cli --input <file> [--limit N]
  2. validate   -> python scripts/validate.py <labelled output file>
  3. roc stats  -> python scripts/roc_stats.py <labelled output file>
  4. start GUI  -> uvicorn app.api.main:app --app-dir backend --port 8000
"""
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path.cwd()
ENV = dict(os.environ, PYTHONPATH="backend")


def find_files(folder, pattern="*.tsv"):
    d = ROOT / folder
    return sorted(d.glob(pattern)) if d.is_dir() else []


def pick(prompt, files):
    """Show a numbered list and return the chosen Path (or None to go back)."""
    if not files:
        print(f"  No files found for: {prompt}")
        manual = input("  Type a path manually (or press Enter to cancel): ").strip()
        return Path(manual) if manual else None
    print(f"\n{prompt}")
    for i, f in enumerate(files, 1):
        try:
            label = f.relative_to(ROOT)
        except ValueError:
            label = f
        print(f"  {i}) {label}")
    print("  m) type a path manually   q) cancel")
    while True:
        ans = input("Choose: ").strip().lower()
        if ans == "q":
            return None
        if ans == "m":
            manual = input("Path: ").strip()
            return Path(manual) if manual else None
        if ans.isdigit() and 1 <= int(ans) <= len(files):
            return files[int(ans) - 1]
        print("  Not a valid choice, try again.")


def run(cmd):
    print("\n$", " ".join(str(c) for c in cmd), "\n")
    try:
        subprocess.run(cmd, env=ENV, check=False)
    except FileNotFoundError as e:
        print("Could not run command:", e)
    except KeyboardInterrupt:
        print("\nStopped.")


def do_classify():
    f = pick("Pick an input TSV to classify:", find_files("data/input"))
    if not f:
        return
    n = input("How many variants? (number, or Enter for all): ").strip()
    cmd = [sys.executable, "-m", "app.cli", "--input", str(f)]
    if n.isdigit() and int(n) > 0:
        cmd += ["--limit", n]
    run(cmd)


def pick_labelled():
    files = find_files("data/output", "*labeled*.tsv") or find_files("data/output")
    return pick("Pick a labelled output TSV:", files)


def do_validate():
    f = pick_labelled()
    if f:
        run([sys.executable, "scripts/validate.py", str(f)])


def do_roc():
    f = pick_labelled()
    if f:
        run([sys.executable, "scripts/roc_stats.py", str(f)])


def do_gui():
    print("GUI will be at http://127.0.0.1:8000  (Ctrl+C to stop)")
    run([sys.executable, "-m", "uvicorn", "app.api.main:app",
         "--app-dir", "backend", "--port", "8000"])


MENU = [
    ("Classify a TSV", do_classify),
    ("Check ACMG/AMP consistency (validate.py)", do_validate),
    ("Compute accuracy stats (roc_stats.py)", do_roc),
    ("Start the browser GUI", do_gui),
]


def main():
    if not (ROOT / "backend").is_dir():
        print("Run this from the project root (the folder containing backend/).")
        sys.exit(1)
    while True:
        print("\n=== Cancer Variant Classifier ===")
        for i, (name, _) in enumerate(MENU, 1):
            print(f"  {i}) {name}")
        print("  q) quit")
        ans = input("Choose: ").strip().lower()
        if ans in ("q", "quit", "exit"):
            break
        if ans.isdigit() and 1 <= int(ans) <= len(MENU):
            MENU[int(ans) - 1][1]()
        else:
            print("  Not a valid choice.")


if __name__ == "__main__":
    main()
