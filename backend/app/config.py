"""Loads <project root>/.env into the environment (NCBI_API_KEY etc.). Import this first.

Before this module existed, .env was never read, so an NCBI_API_KEY placed in
.env was silently ignored and ClinVar stayed at the slow 3 requests/second limit.
"""
from __future__ import annotations

from pathlib import Path

from dotenv import load_dotenv

PROJECT_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(PROJECT_ROOT / ".env", override=False)
