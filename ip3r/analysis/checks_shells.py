"""Paper 6's ligand shells: the pocket measured, then its conservation.

``P6.shell_distances`` recomputes S22's consensus pocket from the six
IP3-bound depositions (all atoms, own subunit, median over deposits). The
other three carry that *recomputed* pocket to each paralog with this
project's own alignment, join it to S17's per-residue conservation and
re-derive the per-shell means, the trend and the absence of a step at the
contact radius (:mod:`ip3r.analysis.shell_constraint`). The published
``shell_constraint.tsv`` and ``shell_trend.tsv`` are only read to compare.
"""

from __future__ import annotations

import math

from ..config import PARALOG_ACC, PARALOGS
from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from ..structure.shells import SHELLS
from .checks import agree, register
from .shell_constraint import DEPOSITS, contact_step, measured_shells, pocket, shell_rows, trend

SHELLS_TSV = "ligand_site/ligand_shells.tsv"
CONSTRAINT = "ligand_site/shell_constraint.tsv"
TREND = "ligand_site/shell_trend.tsv"
PER_RESIDUE = tuple(f"constraint/constraint_{g}_{PARALOG_ACC[g]}.tsv" for g in PARALOGS)


def _p_close(a: float, b: float) -> bool:
    return abs(math.log10(a) - math.log10(b)) <= _P.value("check.log_p_tol")


@register("P6.shell_distances", "ligand",
          "IP3 binding defines a pocket of 125 residues within 15 Å, in four "
          "shells (12 contact, 14 second, 40 third, 59 fourth) by the median "
          "all-atom distance over six ITPR3 depositions (ligand_shells.tsv).",
          "Every residue's all-atom distance to its own subunit's IP3 "
          "re-measured in all six depositions with this project's reader; "
          "median, shell, and the count of depositions placing it in contact "
          "compared residue by residue.",
          "recomputed", (SHELLS_TSV,), DEPOSITS)
def shell_distances():
    pub = {int(r["resi"]): r for r in G.read_tsv(SHELLS_TSV)}
    mine = {s.resi: s for s in measured_shells()}
    tol = _P.value("check.length_tol")
    only_pub, only_mine = sorted(set(pub) - set(mine)), sorted(set(mine) - set(pub))
    bad = [r for r in set(pub) & set(mine)
           if mine[r].shell != pub[r]["shell"]
           or abs(mine[r].median - float(pub[r]["median_distance_A"])) > tol
           or mine[r].n_structures != int(pub[r]["n_structures"])
           or mine[r].n_contact != int(pub[r]["n_structures_contact"])]
    counts = {s: sum(v.shell == s for v in mine.values()) for s in SHELLS}
    worst = max((abs(mine[r].median - float(pub[r]["median_distance_A"]))
                 for r in set(pub) & set(mine)), default=float("nan"))
    found = (f"{len(mine)} residues ({', '.join(f'{n} {s}' for s, n in counts.items())}); "
             f"largest median difference {worst:.1e} Å")
    if only_pub or only_mine or bad:
        found += (f"; only published {only_pub[:8]}, only here {only_mine[:8]}, "
                  f"differing {sorted(bad)[:8]}")
    return agree(not (only_pub or only_mine or bad), f"{len(pub)} residues", found,
                 "Residues are ITPR3 numbers; a deposit's value is the best "
                 "over the subunits whose own site is occupied.", counts=counts)


@register("P6.shell_constraint", "ligand",
          "Every shell out to 15 Å sits above the whole-protein mean "
          "conservation, in all three paralogues (shell_constraint.tsv).",
          "The recomputed pocket carried to ITPR1/ITPR2 by this project's "
          "Gotoh/BLOSUM62 alignment (not S17's MAFFT transfer), joined to the "
          "per-residue deep JSD; per-shell means and a one-sided Mann-Whitney "
          "against every scored residue (own code) compared field by field.",
          "rederived", (CONSTRAINT,) + PER_RESIDUE, DEPOSITS)
def shell_constraint_check():
    pub = {(r["paralog"], r["shell"]): r for r in G.read_tsv(CONSTRAINT)}
    tol, alpha = _P.value("check.stat_tol"), _P.value("check.alpha")
    ok, lines, data = True, [], {}
    for gene in PARALOGS:
        rows, bad, above = shell_rows(gene), [], True
        for r in rows:
            p = pub[(gene, r["shell"])]
            if r["n_residues"] != int(p["n_residues"]):
                bad.append(f"{r['shell']} n")
            bad += [f"{r['shell']} {k}" for k in ("mean_jsd", "median_jsd",
                                                  "mean_frac_modal", "whole_protein_mean_jsd")
                    if abs(r[k] - float(p[k])) > tol]
            if not _p_close(r["p_greater_than_protein"], float(p["p_greater_than_protein"])):
                bad.append(f"{r['shell']} p")
            above &= (r["mean_jsd"] > r["whole_protein_mean_jsd"]
                      and r["p_greater_than_protein"] < alpha)
        ok &= not bad and above
        lines.append(f"{gene}: " + " / ".join(f"{r['mean_jsd']:.3f}" for r in rows)
                     + f" vs protein {rows[0]['whole_protein_mean_jsd']:.4f}"
                     + (f", all p < {alpha:g}" if above else ", NOT all above")
                     + (f" — differs in {', '.join(bad)}" if bad else "")
                     + (f"; {pocket(gene).n_unaligned} unaligned" if pocket(gene).n_unaligned else ""))
        data[gene] = {r["shell"]: r["mean_jsd"] for r in rows}
        data[gene]["protein"] = rows[0]["whole_protein_mean_jsd"]
    return agree(ok, "all 12 shell means above the protein, p < 0.05",
                 "; ".join(lines), "Means are contact / second / third / fourth.",
                 shells=data)


@register("P6.shell_trend", "ligand",
          "Conservation falls with distance from IP3 in all three paralogues, "
          "but shallowly: the Spearman correlation is significant only in "
          "ITPR2 (ρ −0.436; ITPR1 −0.175, p 0.0505; ITPR3 −0.168, p 0.061).",
          "Spearman's ρ of deep JSD against the recomputed median distance "
          "over the carried pocket (own ranks; t-distribution p), compared "
          "with shell_trend.tsv.",
          "rederived", (TREND,) + PER_RESIDUE, DEPOSITS)
def shell_trend():
    pub = {r["paralog"]: r for r in G.read_tsv(TREND)}
    tol, alpha = _P.value("check.stat_tol"), _P.value("check.alpha")
    ok, lines, data, sig = True, [], {}, {}
    for gene in PARALOGS:
        t, p = trend(gene), pub[gene]
        bad = [k for k in ("rho_distance_vs_jsd", "rho_distance_vs_frac_modal")
               if abs(t[k] - float(p[k])) > tol]
        bad += [k for k in ("p_distance_vs_jsd", "p_distance_vs_frac_modal")
                if not _p_close(t[k], float(p[k]))]
        if t["n_residues"] != int(p["n_residues"]):
            bad.append("n")
        # the pattern read off the published table as well as ours
        sig[gene] = (t["p_distance_vs_jsd"] < alpha, float(p["p_distance_vs_jsd"]) < alpha)
        ok &= not bad and t["rho_distance_vs_jsd"] < 0
        lines.append(f"{gene}: ρ {t['rho_distance_vs_jsd']:+.4f}, p {t['p_distance_vs_jsd']:.3g}"
                     + (f" — differs in {', '.join(bad)}" if bad else ""))
        pk = pocket(gene)
        data[gene] = {"distance": pk.distance.tolist(), "jsd": pk.jsd.tolist(),
                      "protein": float(pk.whole_protein.mean()),
                      "rho": t["rho_distance_vs_jsd"]}
    pattern = all(s == (g == "ITPR2") for g, pair in sig.items() for s in pair)
    if not pattern:
        lines.append("significant only in ITPR2 does not hold")
    return agree(ok and pattern, "ρ < 0 in all three; significant in ITPR2 only",
                 "; ".join(lines), trend=data,
                 edges=[_P.value("ligand.contact_cutoff"), _P.value("ligand.shell_second_edge"),
                        _P.value("ligand.shell_third_edge"), _P.value("ligand.shell_radius")])


@register("P6.no_contact_step", "ligand",
          "The step a contact-driven model predicts at 4.5 Å is not there: "
          "the contact shell is no more conserved than the shell just "
          "outside it.",
          "Per paralogue, a one-sided Mann-Whitney of the contact shell's "
          "deep JSD against the second shell's, and the drop in mean JSD "
          "across each shell boundary; the claim holds if the contact shell "
          "does not beat the second and the 4.5 Å boundary does not carry "
          "the largest drop. Not a table in the paper: the prose is tested "
          "on the per-residue input.",
          "rederived", PER_RESIDUE, DEPOSITS)
def no_contact_step():
    alpha = _P.value("check.alpha")
    ok, lines, data = True, [], {}
    for gene in PARALOGS:
        s = contact_step(gene)
        first = f"{SHELLS[0]}→{SHELLS[1]}"
        holds = s["p_contact_gt_second"] >= alpha and s["largest_drop"] != first
        ok &= holds
        lines.append(f"{gene}: drop at 4.5 Å {s['drops'][first]:+.4f} "
                     f"(p {s['p_contact_gt_second']:.2f}); largest at "
                     f"{s['largest_drop']} ({s['drops'][s['largest_drop']]:+.4f})")
        data[gene] = s
    return agree(ok, "no step at 4.5 Å", "; ".join(lines),
                 "12 contact against 14 second-shell residues is a small test, "
                 "so 'not significant' alone would be weak; the drop between "
                 "shells is therefore also located.", steps=data)
