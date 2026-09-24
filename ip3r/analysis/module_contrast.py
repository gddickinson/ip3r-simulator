"""Paper 6's paired module contrast, re-derived from the deep alignments.

For each orthologue in a paralog's deep alignment (S17's
``aln_<gene>.fasta``), one identity to the human reference in the ligand core
and one in the pore module; the test is on the paired difference, so a tip
that is divergent everywhere cancels out.

What is independent of S22 here:

* the **residue → column map** is built by walking the reference row of the
  alignment and checking its ungapped sequence against the UniProt sequence
  this project holds — not read from S17's ``deep_col`` column;
* the **module spans** come from :mod:`ip3r.core.modules`;
* the **statistics** are :func:`stats.sign_test` and
  :func:`stats.signed_rank_test`, written here.

The identity rule is S22's, stated in its prose: over the module's reference
columns that the tip also covers, the fraction with the same residue; a tip
must cover ``ligand.module_min_coverage`` of each module to enter the pair.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from ..config import PARALOG_ACC
from ..core import genes_data as G
from ..core.annotations import reference_sequence
from ..core.modules import module
from ..parameters import PARAMETERS as _P
from .stats import sign_test, signed_rank_test

__all__ = ["alignment_path", "read_alignment", "reference_columns",
           "tip_identities", "PairedContrast", "paired_contrast",
           "AlignmentMismatch", "clear_caches"]


class AlignmentMismatch(RuntimeError):
    """The alignment's reference row is not the reference sequence."""


def alignment_path(paralog: str) -> str:
    return f"constraint/aln_{paralog}.fasta"


@lru_cache(maxsize=None)
def read_alignment(paralog: str) -> tuple[tuple[str, str], ...]:
    """``((label, row), ...)`` in file order."""
    out, label, chunks = [], None, []
    for line in G.read_text(alignment_path(paralog)).splitlines():
        if line.startswith(">"):
            if label is not None:
                out.append((label, "".join(chunks)))
            label, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if label is not None:
        out.append((label, "".join(chunks)))
    return tuple(out)


def clear_caches() -> None:
    """Forget cached alignments (the calibration test swaps the project)."""
    read_alignment.cache_clear()
    reference_columns.cache_clear()


def _reference_row(paralog: str) -> tuple[str, str]:
    acc = PARALOG_ACC[paralog]
    refs = [(lab, row) for lab, row in read_alignment(paralog)
            if lab.startswith("REF|") and acc in lab]
    if len(refs) != 1:
        raise AlignmentMismatch(f"{paralog}: {len(refs)} reference rows for {acc}")
    return refs[0]


@lru_cache(maxsize=None)
def reference_columns(paralog: str) -> dict[int, int]:
    """``residue number -> alignment column`` for the reference row.

    Refuses if the ungapped row is not the canonical sequence: a map built
    on the wrong sequence would put every module in the wrong place.
    """
    _, row = _reference_row(paralog)
    cols = [i for i, ch in enumerate(row) if ch != "-"]
    ungapped = "".join(row[i] for i in cols).upper()
    if ungapped != reference_sequence(paralog).upper():
        raise AlignmentMismatch(
            f"{paralog}: reference row ({len(ungapped)} aa) differs from the "
            f"UniProt sequence ({len(reference_sequence(paralog))} aa)")
    return {k + 1: c for k, c in enumerate(cols)}


def tip_identities(paralog: str, definition: str) -> dict[str, tuple[float, float]]:
    """``tip -> (identity, coverage)`` for one module, every non-reference tip.

    Identity is NaN for a tip covering none of the module.
    """
    ref_label, ref = _reference_row(paralog)
    colmap = reference_columns(paralog)
    cols = np.array([colmap[r] for r in module(paralog, definition).residues])
    ref_aa = np.frombuffer(ref.upper().encode(), "S1")[cols]
    out = {}
    for label, row in read_alignment(paralog):
        if label == ref_label:
            continue
        aa = np.frombuffer(row.upper().encode(), "S1")[cols]
        covered = aa != b"-"
        n = int(covered.sum())
        ident = float((aa[covered] == ref_aa[covered]).sum() / n) if n else float("nan")
        out[label] = (ident, n / len(cols))
    return out


@dataclass
class PairedContrast:
    paralog: str
    core_definition: str
    pore_definition: str
    tips: list = field(default_factory=list)
    core: np.ndarray = field(default_factory=lambda: np.zeros(0))
    pore: np.ndarray = field(default_factory=lambda: np.zeros(0))
    dropped: list = field(default_factory=list)   # (tip, lost_on)

    @property
    def diff(self) -> np.ndarray:
        return self.core - self.pore

    @property
    def mean_difference(self) -> float:
        return float(self.diff.mean())

    @property
    def direction(self) -> str:
        return "core > pore" if self.mean_difference > 0 else "pore > core"

    def stats(self) -> dict:
        s, w = sign_test(self.diff), signed_rank_test(self.diff)
        return {"n_tips": len(self.tips), "n_dropped": len(self.dropped),
                "mean_core_identity": float(self.core.mean()),
                "mean_pore_identity": float(self.pore.mean()),
                "mean_difference": self.mean_difference,
                "n_core_more_conserved": s["n_pos"],
                "n_pore_more_conserved": s["n_neg"], "n_ties": s["n_ties"],
                "p_sign": s["p"], "p_wilcoxon": w["p"],
                "direction": self.direction}


def paired_contrast(paralog: str, core_definition: str = "contact_span",
                    pore_definition: str = "channel_minus_luminal") -> PairedContrast:
    """The paired test's inputs: tips resolving both modules, and the rest."""
    floor = _P.value("ligand.module_min_coverage")
    a = tip_identities(paralog, core_definition)
    b = tip_identities(paralog, pore_definition)
    out = PairedContrast(paralog, core_definition, pore_definition)
    core, pore = [], []
    for tip, (ia, ca) in a.items():
        ib, cb = b[tip]
        if not np.isfinite(ia) or ca < floor:
            out.dropped.append((tip, core_definition))
        elif not np.isfinite(ib) or cb < floor:
            out.dropped.append((tip, pore_definition))
        else:
            out.tips.append(tip)
            core.append(ia)
            pore.append(ib)
    out.core, out.pore = np.array(core), np.array(pore)
    return out
