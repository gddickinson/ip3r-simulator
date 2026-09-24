"""Paper 2's tree claims, re-derived from ``rooted.nwk`` with :mod:`.tree`.

Each check reads the published numbers from the table that states them and
recomputes them from the tree with this project's reader and its own
definitions (a paralogue's clade is the MRCA of every record naming it;
cyclostomes are recognised by genus; support is counted per bipartition).
The Newick text travels in the outcome so the exhibit draws the tree the
verdict was reached on.
"""

from __future__ import annotations

from ..config import PARALOGS
from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .checks import agree, register
from .tree import (TREE, bipartitions, cyclostome_clades, is_cyclostome, paralog_clades,
                   support, vertebrate_tips)
from .newick import parse

CLADES = "phylogeny/paralog_clades.tsv"
CYCLO = "phylogeny/cyclostome_placement.tsv"
SUPPORT = "phylogeny/support_summary.tsv"


def _fmt(s) -> str:
    return "unlabelled" if s is None else f"{s[0]:g}/{s[1]:g}"


@register("P2.paralog_clades", "origin",
          "ITPR1, ITPR2 and ITPR3 are each a clade at 100/100 holding 19, 13 "
          "and 19 proteins; seven unnamed shark, chimaera and coelacanth genes "
          "fall inside them (4/2/1), and the only vertebrate genes outside "
          "every clade are the six from hagfish and lamprey.",
          "For each paralogue, the MRCA in rooted.nwk of every tip whose "
          "record names it; its size, SH-aLRT/UFBoot label, unnamed "
          "vertebrate tips and any tip that does not belong (another "
          "paralogue, a cyclostome, an invertebrate). Compared with "
          "paralog_clades.tsv; the vertebrates left outside all three clades "
          "are listed.",
          "rederived", (TREE, CLADES))
def clades():
    text = G.read_text(TREE)
    root = parse(text)
    pub = {r["paralog"]: r for r in G.read_tsv(CLADES)}
    cl = paralog_clades(root)
    ok, pub_s, found_s = True, [], []
    for g in PARALOGS:
        c, p = cl[g], pub[g]
        want = (int(p["n_extended"]), (float(p["alrt"]), float(p["ufboot"])),
                sorted(x for x in p["added"].split(";") if x))
        got = (len(c.tips), c.support, c.unnamed)
        ok &= got == want and not c.foreign
        pub_s.append(f"{g} {want[0]} at {_fmt(want[1])} (+{len(want[2])} unnamed)")
        found_s.append(f"{g} {got[0]} at {_fmt(got[1])} (+{len(got[2])} unnamed"
                       + (f", {len(c.foreign)} foreign" if c.foreign else "") + ")")
    inside = set().union(*(c.tips for c in cl.values()))
    outside = sorted(vertebrate_tips(root) - inside)
    ok &= bool(outside) and all(is_cyclostome(t) for t in outside) and len(outside) == 6
    found_s.append(f"{len(outside)} vertebrate tips outside, "
                   f"{sum(map(is_cyclostome, outside))} of them cyclostome")
    return agree(ok, "; ".join(pub_s) + "; 6 cyclostome tips outside",
                 "; ".join(found_s), newick=text)


@register("P2.cyclostome_lineages", "origin",
          "The six hagfish and lamprey loci form two cyclostome-only clades, "
          "each holding both species: four branch first among the vertebrates "
          "(99.5/100), and two join the ITPR2 + ITPR3 clade (83.1/77).",
          "Maximal clades of cyclostome tips (recognised by genus) in "
          "rooted.nwk; for each, the genera it holds, its support, and what "
          "its sister subtree contains. The published supports are read "
          "from cyclostome_placement.tsv by the clade size it reports.",
          "rederived", (TREE, CYCLO))
def cyclostomes():
    text = G.read_text(TREE)
    root = parse(text)
    rows = G.read_tsv(CYCLO)
    pub = {int(r["clade_size"]): (float(r["alrt"]), float(r["ufboot"])) for r in rows}
    cc = cyclostome_clades(root)
    both = [c for c in cc if len(c.genera) >= 2]
    first = [c for c in cc if c.first_among_vertebrates]
    joins = [c for c in cc if set(c.sister) <= {"ITPR2", "ITPR3", "vertebrate_basal"}
             and c.sister["ITPR2"] and c.sister["ITPR3"]]
    found = [f"{len(cc)} cyclostome-only clade(s) of sizes "
             f"{'/'.join(str(len(c.tips)) for c in cc)}, "
             f"{len(both)} holding both genera"]
    ok = len(rows) == 6 and len(cc) == 2 and len(both) == 2
    if first:
        c = first[0]
        found.append(f"{len(c.tips)} branch first among the vertebrates at {_fmt(c.support)}")
        ok &= pub.get(len(c.tips)) == c.support
    if joins:
        c = joins[0]
        n = len(c.tips) + sum(c.sister.values())
        found.append(f"{len(c.tips)} join ITPR2 + ITPR3 (a {n}-tip clade) at "
                     f"{_fmt(c.parent_support)}")
        ok &= pub.get(n) == c.parent_support
    ok &= len(first) == 1 and len(joins) == 1 and first[0] is not joins[0]
    published = (f"{len(rows)} loci; 2 clades, both species in each; "
                 + "; ".join(f"{k}-tip clade at {_fmt(v)}" for k, v in sorted(pub.items())))
    return agree(ok, published, "; ".join(found), newick=text)


@register("P2.support_bar", "origin",
          "Of the 131 internal nodes, 91 (69.5 %) clear both SH-aLRT ≥ 80 and "
          "UFBoot ≥ 95.",
          "Every labelled bipartition of rooted.nwk (the root's two edges "
          "are one bipartition and counted once), tested against the "
          "registered bars tree.alrt_min and tree.ufboot_min; compared with "
          "support_summary.tsv.",
          "rederived", (TREE, SUPPORT))
def support_bar():
    text = G.read_text(TREE)
    parts = [support(n) for n in bipartitions(parse(text))]
    a_min, u_min = _P.value("tree.alrt_min"), _P.value("tree.ufboot_min")
    got = {"internal nodes": len(parts),
           "SH-aLRT >= 80": sum(a >= a_min for a, _ in parts),
           "UFBoot >= 95": sum(u >= u_min for _, u in parts),
           "both thresholds": sum(a >= a_min and u >= u_min for a, u in parts)}
    pub = {r["statistic"]: r["value"] for r in G.read_tsv(SUPPORT)}
    want = {k: int(float(pub[k])) for k in got}
    ok = got == want
    return agree(ok, ", ".join(f"{k} {v}" for k, v in want.items()),
                 ", ".join(f"{k} {v}" for k, v in got.items()),
                 f"{got['both thresholds'] / got['internal nodes']:.1%} clear both",
                 newick=text)
