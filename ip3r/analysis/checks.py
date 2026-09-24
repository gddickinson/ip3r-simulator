"""The findings checks: every published result this application re-derives.

A **check** pairs a statement ``ip3r_genes`` publishes with an *independent*
route to the same number, and reports whether the two agree. The routes come
in three kinds, and the kind is part of the result because they are not
equally strong:

* ``recomputed`` — from primary data with this project's own code: the
  deposited coordinates for the structural results. The strongest kind; the
  two projects share no code on this path.
* ``rederived`` — from the publication's committed *input* tables (per-residue
  conservation, per-variant scores, the tree), with this project's own
  arithmetic, compared with its committed *output* table. It checks the step
  from input to output, not the inputs.
* ``read`` — the published value is read and tested against the statement
  made about it (e.g. "only one AU hypothesis survives"). It checks that the
  prose matches the table, and nothing more.

Outcomes: ``confirmed``, ``discrepancy``, ``not_run`` (data or structure
unavailable — never a scientific result), ``error`` (the check itself broke).

A check that cannot say "no" is not a check, so every check function is
exercised in ``tests/test_checks_calibration.py`` on a planted discrepancy.
And checks refuse to report ``confirmed`` against a modified parameter
registry: the published numbers were produced at the defaults.
"""

from __future__ import annotations

import time
import traceback
from dataclasses import dataclass, field
from typing import Callable

from ..core.genes_data import GenesDataMissing
from ..parameters import PARAMETERS

__all__ = ["Outcome", "Check", "CheckResult", "REGISTRY", "register",
           "run_check", "run_checks", "PAPERS", "confirmed", "discrepancy",
           "not_run", "agree", "all_checks"]

#: The six papers of the ip3r_genes series, plus the structural baseline.
PAPERS = {
    "structure": "S0 baseline — the channel measured on 6DQN",
    "range": "Paper 1 — ancestrally eukaryotic, lost repeatedly outside animals",
    "origin": "Paper 2 — two separate early-vertebrate duplications",
    "retention": "Paper 3 — no vertebrate has lost a paralogue",
    "archive": "Paper 4 — four in five genes unreachable from protein databases",
    "constraint": "Paper 5 — selection strongest at gate and filter",
    "ligand": "Paper 6 — binding core not more constrained than the pore",
    "ledger": "The publication's own claims ledgers",
}


@dataclass
class Outcome:
    status: str                      # confirmed | discrepancy | not_run | error
    published: str = ""
    found: str = ""
    detail: str = ""
    data: dict = field(default_factory=dict)


def confirmed(published, found, detail="", **data) -> Outcome:
    return Outcome("confirmed", str(published), str(found), detail, data)


def discrepancy(published, found, detail="", **data) -> Outcome:
    return Outcome("discrepancy", str(published), str(found), detail, data)


def not_run(detail: str) -> Outcome:
    return Outcome("not_run", detail=detail)


def agree(ok: bool, published, found, detail="", **data) -> Outcome:
    return (confirmed if ok else discrepancy)(published, found, detail, **data)


@dataclass(frozen=True)
class Check:
    id: str
    paper: str
    claim: str
    method: str
    kind: str                        # recomputed | rederived | read
    sources: tuple = ()
    structures: tuple = ()
    func: Callable[[], Outcome] | None = None


@dataclass
class CheckResult:
    check: Check
    outcome: Outcome
    seconds: float

    @property
    def status(self) -> str:
        return self.outcome.status


REGISTRY: dict[str, Check] = {}


def register(id: str, paper: str, claim: str, method: str, kind: str,
             sources=(), structures=()):
    """Decorator: add a check function to the registry."""
    if paper not in PAPERS:
        raise ValueError(f"unknown paper {paper!r}")
    if kind not in ("recomputed", "rederived", "read"):
        raise ValueError(f"unknown kind {kind!r}")

    def wrap(fn: Callable[[], Outcome]):
        if id in REGISTRY:
            raise ValueError(f"duplicate check id {id}")
        REGISTRY[id] = Check(id, paper, claim, method, kind, tuple(sources),
                             tuple(structures), fn)
        return fn
    return wrap


def _load_all() -> None:
    # Importing the modules registers their checks.
    from . import (checks_constraint, checks_evolution, checks_modules,  # noqa: F401
                   checks_shells, checks_structure)


def run_check(check: Check, allow_modified: bool = False) -> CheckResult:
    t0 = time.perf_counter()
    if PARAMETERS.modified and not allow_modified:
        out = not_run("parameters differ from their documented defaults ("
                      + PARAMETERS.override_summary() + "); the published "
                      "numbers were produced at the defaults")
        return CheckResult(check, out, 0.0)
    try:
        out = check.func()
    except GenesDataMissing as exc:
        out = not_run(str(exc))
    except FileNotFoundError as exc:
        out = not_run(f"missing file: {exc}")
    except Exception as exc:                          # the check itself broke
        out = Outcome("error", detail=f"{type(exc).__name__}: {exc}",
                      data={"traceback": traceback.format_exc()})
    return CheckResult(check, out, time.perf_counter() - t0)


def run_checks(ids=None, paper: str | None = None, allow_modified: bool = False,
               progress=None) -> list[CheckResult]:
    """Run the selected checks (all by default) in registry order."""
    _load_all()
    todo = [c for c in REGISTRY.values()
            if (ids is None or c.id in ids) and (paper is None or c.paper == paper)]
    results = []
    for i, c in enumerate(todo):
        if progress:
            progress(i, len(todo), c)
        results.append(run_check(c, allow_modified))
    return results


def all_checks() -> list[Check]:
    _load_all()
    return list(REGISTRY.values())
