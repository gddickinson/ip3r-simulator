"""Paper 6 §8: the ligand pocket and the two modules asked of the
substitution rate.

S17's FEL run (``constraint/fel_sites.tsv``: per codon site α, β, ω, q and
FEL's verdict, in each paralog's UniProt numbering) is the committed *input*.
Everything joined to it is this project's: the pocket is the one
:mod:`ip3r.analysis.shell_constraint` measured on six ITPR3 deposits and
carried by its own alignment, the modules are :mod:`ip3r.core.modules`, and
the rank tests and q-values are :mod:`ip3r.analysis.stats`.

As S22 does, β is compared by rank (it is zero at most sites), a site is
purifying at FEL's own q ≤ ``check.fel_q``, and ω is summarised only over
sites whose α is not at FEL's bound.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from ..core import genes_data as G
from ..core.modules import CORE_DEFINITIONS, module
from ..parameters import PARAMETERS as _P
from ..structure.shells import SHELLS
from .shell_constraint import pocket
from .stats import benjamini_hochberg, mann_whitney_less, mann_whitney_two_sided

__all__ = ["FEL_TSV", "Site", "fel_sites", "shell_rates", "MODULE_PAIRS",
           "module_rates", "clear_caches"]

FEL_TSV = "constraint/fel_sites.tsv"

#: S22's module comparisons: (core, pore, primary), in its table's order.
MODULE_PAIRS = (("contact_span", "channel_minus_luminal", True),
                ("contact_span", "channel_all", False),
                ("ibc_literature", "channel_minus_luminal", False))
assert {c for c, _, _ in MODULE_PAIRS} == set(CORE_DEFINITIONS)


@dataclass(frozen=True)
class Site:
    resi: int
    alpha: float
    beta: float
    omega: float
    q: float
    verdict: str
    alpha_at_bound: bool


@lru_cache(maxsize=None)
def fel_sites(gene: str) -> dict[int, Site]:
    """S17's FEL sites of one paralog, keyed by UniProt residue."""
    out = {}
    for r in G.read_tsv(FEL_TSV):
        if r["paralog"] != gene or r["resi"] in ("", "None"):
            continue
        out[int(r["resi"])] = Site(int(r["resi"]), G.as_float(r["alpha"]),
                                   G.as_float(r["beta"]), G.as_float(r["omega"]),
                                   G.as_float(r["q_value"]), r["verdict"],
                                   r["alpha_at_bound"] == "True")
    return out


def _summary(sites: list[Site]) -> dict:
    beta = np.array([s.beta for s in sites], float)
    usable = [s.omega for s in sites if not s.alpha_at_bound
              and np.isfinite(s.omega) and np.isfinite(s.alpha)]
    q = _P.value("check.fel_q")
    pur = sum(1 for s in sites if np.isfinite(s.q) and s.q <= q and s.verdict == "purifying")
    fin = beta[np.isfinite(beta)]
    return {"n_sites": len(sites),
            "median_beta": float(np.median(fin)) if len(fin) else float("nan"),
            "mean_beta": float(fin.mean()) if len(fin) else float("nan"),
            "n_beta_zero": int((fin == 0.0).sum()),
            "n_usable_alpha": len(usable),
            "n_alpha_at_bound": sum(s.alpha_at_bound for s in sites),
            "median_omega_usable": float(np.median(usable)) if usable else float("nan"),
            "frac_purifying_q05": pur / len(sites) if sites else float("nan")}


def _betas(sites) -> np.ndarray:
    return np.array([s.beta for s in sites], float)


def shell_rates(genes) -> list[dict]:
    """One row per (paralog, shell): S22's ``omega_by_shell.tsv`` fields.
    Each shell's β against the whole protein's (one-sided, lower), with
    q-values over every row asked for together, as S22 corrects them."""
    rows = []
    for gene in genes:
        fel, pk = fel_sites(gene), pocket(gene)
        protein = _betas(fel.values())
        for name in SHELLS:
            sub = [fel[int(r)] for r in sorted(pk.resi[pk.in_shell(name)]) if int(r) in fel]
            row = {"paralog": gene, "shell": name, **_summary(sub)}
            row["whole_protein_median_beta"] = float(np.nanmedian(protein))
            row["p_less_than_protein"] = mann_whitney_less(_betas(sub), protein)["p"]
            rows.append(row)
    for r, q in zip(rows, benjamini_hochberg([r["p_less_than_protein"] for r in rows])):
        r["q_less_than_protein"] = float(q)
    return rows


def module_rates(genes) -> list[dict]:
    """Core against pore on β for every pair in :data:`MODULE_PAIRS` and
    paralog (two-sided rank test; q over all rows): ``omega_module_test.tsv``."""
    rows = []
    for core_def, pore_def, primary in MODULE_PAIRS:
        for gene in genes:
            fel = fel_sites(gene)
            a = _betas(fel[r] for r in module(gene, core_def).residues if r in fel)
            b = _betas(fel[r] for r in module(gene, pore_def).residues if r in fel)
            a, b = a[np.isfinite(a)], b[np.isfinite(b)]
            mw = mann_whitney_two_sided(a, b)
            rows.append({"paralog": gene, "core_definition": core_def,
                         "pore_definition": pore_def, "is_primary": primary,
                         "n_core_sites": len(a), "n_pore_sites": len(b),
                         "mean_core_beta": float(a.mean()),
                         "mean_pore_beta": float(b.mean()),
                         "cles_core_gt_pore": mw["cles"], "p_mannwhitney": mw["p"],
                         "direction": ("core evolves faster" if a.mean() > b.mean()
                                       else "pore evolves faster")})
    for r, q in zip(rows, benjamini_hochberg([r["p_mannwhitney"] for r in rows])):
        r["q_mannwhitney"] = float(q)
    return rows


def clear_caches() -> None:
    fel_sites.cache_clear()
