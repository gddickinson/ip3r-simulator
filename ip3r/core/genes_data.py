"""Read-only access to the ``ip3r_genes`` committed results tables.

The findings checks and exhibits read the publication project's tables
**live** — the point of a check is to recompute a published number from the
table it rests on, so reading a copy would check the copy. Every loader
raises :class:`GenesDataMissing` (not ``FileNotFoundError``) when the project
or a table is absent, which the check runner turns into ``not_run`` rather
than a failure: a missing sibling checkout is not a scientific result.

Parsing is deliberately stdlib-only (``csv`` + ``json``), so the loaders have
no dependency the checks could inherit an error from.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from ..config import genes_results

__all__ = ["GenesDataMissing", "results_dir", "path", "read_tsv",
           "read_json", "read_text", "available", "sha256", "as_float",
           "as_bool"]


class GenesDataMissing(RuntimeError):
    """The ip3r_genes project, or one of its tables, is not available."""


def results_dir() -> Path:
    d = genes_results()
    if not d.is_dir():
        raise GenesDataMissing(
            f"ip3r_genes results not found at {d} — clone ip3r_genes beside "
            f"this project or set IP3R_GENES_DIR")
    return d


def path(rel: str) -> Path:
    p = results_dir() / rel
    if not p.exists():
        raise GenesDataMissing(f"ip3r_genes table missing: results/{rel}")
    return p


def available(rel: str) -> bool:
    try:
        path(rel)
        return True
    except GenesDataMissing:
        return False


def read_tsv(rel: str) -> list[dict]:
    with open(path(rel), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(rel: str):
    return json.loads(path(rel).read_text(encoding="utf-8"))


def read_text(rel: str) -> str:
    return path(rel).read_text(encoding="utf-8")


def sha256(rel: str) -> str:
    return hashlib.sha256(path(rel).read_bytes()).hexdigest()


def as_float(value: str) -> float:
    """Parse a table cell; blank, NA and nan all become NaN."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def as_bool(value: str) -> bool:
    return str(value).strip() in {"True", "true", "1", "yes"}
