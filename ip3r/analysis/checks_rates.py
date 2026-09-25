"""Paper 6 §8 (S22): the pocket and the modules asked of FEL's β.

Both checks join S17's per-site FEL rates (the committed input) to structure
this project built: the pocket recomputed from six ITPR3 deposits and carried
by its own alignment, and the modules rebuilt from the imported annotation.
S22's ``omega_by_shell.tsv`` and ``omega_module_test.tsv`` are read only to
compare with (:mod:`ip3r.analysis.fel_rates`).
"""

from __future__ import annotations

import math

from ..config import PARALOGS
from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .checks import agree, register
from .checks_shells import PER_RESIDUE
from .fel_rates import FEL_TSV, module_rates, shell_rates
from .shell_constraint import DEPOSITS

BY_SHELL = "ligand_site/omega_by_shell.tsv"
MODULE_TEST = "ligand_site/omega_module_test.tsv"

_COUNTS = ("n_sites", "n_beta_zero", "n_usable_alpha", "n_alpha_at_bound")


def _p_close(a: float, b: float) -> bool:
    return abs(math.log10(a) - math.log10(b)) <= _P.value("check.log_p_tol")


def _differs(mine: dict, theirs: dict, counts, values, pvalues) -> list[str]:
    tol = _P.value("check.stat_tol")
    bad = [k for k in counts if int(mine[k]) != int(theirs[k])]
    bad += [k for k in values if abs(mine[k] - float(theirs[k])) > tol]
    bad += [k for k in pvalues if not _p_close(mine[k], float(theirs[k]))]
    return bad


@register("P6.shell_rates", "ligand",
          "The shell result replicates on the substitution rate: FEL calls "
          "every contact-shell site purifying (q ≤ 0.05) in all three "
          "paralogues, and the share steps down past the second shell to a "
          "floor. It is not monotone: the fourth shell sits a little above "
          "the third in all three (omega_by_shell.tsv, S22 §8).",
          "S17's per-site FEL rates joined to the pocket recomputed here "
          "(six ITPR3 deposits, carried by this project's alignment); per "
          "shell the purifying share, β summaries and a one-sided rank test "
          "against the whole protein with BH q over the twelve rows (own "
          "code), compared field by field; the pattern tested on ours.",
          "rederived", (FEL_TSV, BY_SHELL) + PER_RESIDUE, DEPOSITS)
def shell_rates_check():
    pub = {(r["paralog"], r["shell"]): r for r in G.read_tsv(BY_SHELL)}
    rows = shell_rates(PARALOGS)
    ok, lines, data = True, [], {}
    for gene in PARALOGS:
        mine = {r["shell"]: r for r in rows if r["paralog"] == gene}
        bad = [f"{s} {k}" for s, r in mine.items()
               for k in _differs(r, pub[(gene, s)], _COUNTS,
                                 ("mean_beta", "median_beta", "frac_purifying_q05"),
                                 ("p_less_than_protein", "q_less_than_protein"))]
        f = {s: r["frac_purifying_q05"] for s, r in mine.items()}
        shape = (f["contact"] == 1.0 and f["second"] > f["third"]
                 and f["fourth"] > f["third"])
        ok &= shape and not bad
        lines.append(f"{gene}: " + " / ".join(f"{v:.3f}" for v in f.values())
                     + f", drop {f['contact'] - min(f.values()):.3f}"
                     + ("" if shape else " — NOT the stated shape")
                     + (f" — differs in {', '.join(bad)}" if bad else ""))
        data[gene] = f
    return agree(ok, "contact 1.000 everywhere; a step past the second shell, "
                 "fourth above third", "; ".join(lines),
                 "Shares are contact / second / third / fourth, at FEL's q ≤ "
                 f"{_P.value('check.fel_q'):g}. β is zero at most sites, so "
                 "the rank tests are on ties.", fractions=data)


@register("P6.module_rates", "ligand",
          "The module comparison does not replicate on β: the pore evolves "
          "faster than the ligand core in ITPR1 (q 6.17e-08), and neither "
          "direction is significant in ITPR2 (pore faster, q 0.317) or ITPR3 "
          "(core faster, q 0.058) (omega_module_test.tsv, S22 §8).",
          "S17's per-site FEL β over both modules rebuilt here (the primary "
          "pair and S22's two sensitivity definitions, the literature core "
          "carried by this project's alignment); two-sided rank test, "
          "common-language effect and BH q over the nine rows (own code), "
          "compared field by field; the direction and significance per "
          "paralogue tested on ours.",
          "rederived", (FEL_TSV, MODULE_TEST))
def module_rates_check():
    key = ("paralog", "core_definition", "pore_definition")
    pub = {tuple(r[k] for k in key): r for r in G.read_tsv(MODULE_TEST)}
    alpha = _P.value("check.alpha")
    ok, lines, data = True, [], {}
    for r in module_rates(PARALOGS):
        p = pub[tuple(r[k] for k in key)]
        bad = _differs(r, p, ("n_core_sites", "n_pore_sites"),
                       ("mean_core_beta", "mean_pore_beta", "cles_core_gt_pore"),
                       ("p_mannwhitney", "q_mannwhitney"))
        if r["direction"] != p["direction"]:
            bad.append("direction")
        ok &= not bad
        if bad:
            lines.append(f"{'/'.join(r[k] for k in key)} differs in {', '.join(bad)}")
        if r["is_primary"]:
            data[r["paralog"]] = r
    # the prose's pattern, on the re-derived primary rows
    pattern = (data["ITPR1"]["direction"] == "pore evolves faster"
               and data["ITPR1"]["q_mannwhitney"] < alpha
               and all(data[g]["q_mannwhitney"] >= alpha for g in ("ITPR2", "ITPR3"))
               and data["ITPR2"]["direction"] == "pore evolves faster"
               and data["ITPR3"]["direction"] == "core evolves faster")
    if not pattern:
        lines.append("the stated directions and significance do not hold")
    found = "; ".join(
        f"{g}: core β {d['mean_core_beta']:.4f} vs pore {d['mean_pore_beta']:.4f}, "
        f"{d['direction']}, q {d['q_mannwhitney']:.3g}" for g, d in data.items())
    return agree(ok and pattern, "ITPR1 pore faster (q 6.2e-08); ITPR2, ITPR3 not "
                 "significant", "; ".join([found] + lines),
                 "Primary pair: the span of the ten contacts against PF00520 less "
                 "the luminal loop. All nine rows (three definitions) are "
                 "compared, since their q-values are corrected together.",
                 modules={g: {k: d[k] for k in ("mean_core_beta", "mean_pore_beta",
                                                 "q_mannwhitney")}
                          for g, d in data.items()})
