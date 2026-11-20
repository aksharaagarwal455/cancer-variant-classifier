"""
Canonical Variant model.

Strongly-typed, normalized representation used by every downstream
stage. Built by the parser from whatever raw schema the input file
uses; always keeps the untouched original row so nothing is silently
discarded.
"""

from __future__ import annotations

from enum import Enum
from typing import Any, Optional

from pydantic import BaseModel, Field


class GenomeBuild(str, Enum):
    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"
    UNKNOWN = "GENOME_BUILD_UNKNOWN"


class ParseStatus(str, Enum):
    OK = "OK"
    INVALID = "INVALID"
    DUPLICATE = "DUPLICATE"


class Variant(BaseModel):
    variant_id: str = Field(..., description="Deterministic ID: chrom:pos:ref:alt[:transcript]")
    original_row_id: int

    chromosome: Optional[str] = None
    position: Optional[int] = None
    ref: Optional[str] = None
    alt: Optional[str] = None
    genome_build: GenomeBuild = GenomeBuild.UNKNOWN

    gene: Optional[str] = None
    transcript: Optional[str] = None
    hgvsg: Optional[str] = None
    hgvsc: Optional[str] = None
    hgvsp: Optional[str] = None
    rsid: Optional[str] = None

    variant_type: Optional[str] = None
    consequence: Optional[str] = None

    parse_status: ParseStatus = ParseStatus.OK
    parse_errors: list[str] = Field(default_factory=list)

    original_input: dict[str, Any] = Field(default_factory=dict)

    model_config = {"use_enum_values": False}
