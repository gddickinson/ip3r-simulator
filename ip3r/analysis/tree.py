"""Paper 2's tree, read with this project's Newick reader.

What the origin paper says about its tree is a set of clade questions: which
tips form each paralogue's clade and at what support, where the jawless-fish
(cyclostome) loci sit, and how many nodes clear the joint support bar. They
are answered here from ``rooted.nwk`` alone, with tip groups read from the
census prefix every label carries (``ITPR1_…``, ``vertebrate_basal_…``,
``RYR_…``) and cyclostomes recognised by genus, so no ip3r_genes clade table
enters the answer. The checks in ``checks_tree`` compare these answers with
the published tables; ``tree_figure`` draws them.

A rooted tree's root joins two edges that are one bipartition of the
unrooted tree (IQ-TREE labels both), so support is counted per bipartition,
not per node — the first run counted 92 of 132 where the paper has 91 of 131.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from ..config import PARALOGS
from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .newick import Node, leaves, mrca, parse

__all__ = ["TREE", "GROUPS", "CYCLOSTOME_GENERA", "load_tree", "group_of",
           "is_cyclostome", "genus_of", "support", "is_supported", "bipartitions",
           "Clade", "paralog_clades", "outgroup_clade", "CyclostomeClade",
           "cyclostome_clades", "vertebrate_tips"]

TREE = "phylogeny/rooted.nwk"
#: Census groups, in the order their prefixes are tried (longest first).
GROUPS = ("vertebrate_basal", "invert_metazoa", "ITPR1", "ITPR2", "ITPR3",
          "protist", "plant", "fungi", "RYR")
#: Hagfish and lamprey genera; a tip is a cyclostome if one is in its label.
CYCLOSTOME_GENERA = ("Myxine", "Eptatretus", "Petromyzon", "Lethenteron",
                     "Entosphenus", "Geotria", "Mordacia")
_VERTEBRATE = frozenset(("ITPR1", "ITPR2", "ITPR3", "vertebrate_basal"))


def load_tree() -> Node:
    return parse(G.read_text(TREE))


def group_of(label: str) -> str:
    for g in GROUPS:
        if label.startswith(g + "_"):
            return g
    return "other"


def genus_of(label: str) -> str | None:
    parts = label.split("_")
    return next((p for p in parts if p in CYCLOSTOME_GENERA), None)


def is_cyclostome(label: str) -> bool:
    return genus_of(label) is not None


def support(node: Node) -> tuple[float, float] | None:
    """(SH-aLRT, UFBoot) from an IQ-TREE ``a/b`` node label, else None."""
    parts = node.label.split("/")
    if len(parts) != 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


def is_supported(node: Node) -> bool:
    s = support(node)
    return s is not None and s[0] >= _P.value("tree.alrt_min") \
        and s[1] >= _P.value("tree.ufboot_min")


def bipartitions(root: Node) -> list[Node]:
    """One labelled internal node per bipartition of the unrooted tree."""
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        stack.extend(n.children)
        if n is not root and n.children and support(n) is not None:
            out.append(n)
    if len(root.children) == 2:                     # the rooted edge, twice
        twin = root.children[1]
        out = [n for n in out if n is not twin] if support(root.children[0]) else out
    return out


def vertebrate_tips(root: Node) -> set[str]:
    return {lf.label for lf in leaves(root) if group_of(lf.label) in _VERTEBRATE}


@dataclass
class Clade:
    name: str
    node: Node = field(repr=False)
    tips: list[str]
    named: list[str]            # records naming this paralogue
    unnamed: list[str]          # vertebrate records naming none
    foreign: list[str]          # anything else (another paralogue, a cyclostome, …)

    @property
    def support(self) -> tuple[float, float] | None:
        return support(self.node)


def _clade(name: str, node: Node) -> Clade:
    tips = sorted(lf.label for lf in leaves(node))
    named = [t for t in tips if group_of(t) == name]
    unnamed = [t for t in tips if group_of(t) == "vertebrate_basal" and not is_cyclostome(t)]
    foreign = [t for t in tips if t not in named and t not in unnamed]
    return Clade(name, node, tips, named, unnamed, foreign)


def paralog_clades(root: Node) -> dict[str, Clade]:
    """Each paralogue's whole clade: the MRCA of every tip whose record names it."""
    out = {}
    for g in PARALOGS:
        tips = {lf.label for lf in leaves(root) if group_of(lf.label) == g}
        out[g] = _clade(g, mrca(root, tips))
    return out


def outgroup_clade(root: Node) -> Clade:
    tips = {lf.label for lf in leaves(root) if group_of(lf.label) == "RYR"}
    return _clade("RYR", mrca(root, tips))


@dataclass
class CyclostomeClade:
    node: Node = field(repr=False)
    tips: list[str]
    genera: list[str]
    support: tuple[float, float] | None
    parent_support: tuple[float, float] | None
    sister: Counter           # census groups of the sister subtree
    first_among_vertebrates: bool


def cyclostome_clades(root: Node) -> list[CyclostomeClade]:
    """The maximal clades holding only cyclostome tips (a lone tip is its own)."""
    cyc = {lf.label for lf in leaves(root) if is_cyclostome(lf.label)}
    verts = vertebrate_tips(root)
    found: list[Node] = []
    for lf in leaves(root):
        if lf.label not in cyc:
            continue
        n = lf
        while n.parent is not None and all(x.label in cyc for x in leaves(n.parent)):
            n = n.parent
        if not any(n is f for f in found):
            found.append(n)
    out = []
    for n in found:
        tips = sorted(lf.label for lf in leaves(n))
        p = n.parent
        sister = Counter(group_of(lf.label) for s in (p.children if p else [])
                         if s is not n for lf in leaves(s))
        under_parent = {lf.label for lf in leaves(p)} if p else set()
        out.append(CyclostomeClade(
            n, tips, sorted({genus_of(t) for t in tips}), support(n),
            support(p) if p else None, sister,
            bool(p) and verts <= under_parent and under_parent <= verts))
    return sorted(out, key=lambda c: -len(c.tips))
