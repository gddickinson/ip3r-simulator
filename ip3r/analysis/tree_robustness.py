"""Paper 2's clade claims re-asked of the ``--bnni`` tree (S7 §5.7).

``--bnni`` gives every UFBoot replicate an extra round of NNI, stripping the
support that came from the model rather than the data. The question is
whether the claims the reported tree carries survive it.

Nothing here reads S7's claim tables to decide what a claim *is*. Each
claim's tip set is rebuilt from the census prefixes by this project's own
rules (:mod:`.tree`): a paralogue's core is its largest pure clade (every
tip's record names it; S7 §5.1 shows five shark and coelacanth records
falling outside it), its clade is the MRCA of every record naming it, the pairs and the triple are unions of the
clades, the vertebrates are every vertebrate tip, and the family is every
tip but the RyR outgroup. Both trees are rooted on that outgroup with
:func:`.newick.reroot`, and a claim holds where a node's tips are exactly
its set. ``claim_members.tsv`` is read only to compare the sets.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..config import PARALOGS
from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from .newick import Node, clade_node, leaves, parse, reroot
from .tree import TREE, group_of, paralog_clades, support, vertebrate_tips

__all__ = ["BNNI_TREE", "load_pair", "largest_pure_clade", "claim_sets", "Claim", "compare", "strong",
           "VERDICTS"]

BNNI_TREE = "phylogeny/itpr_ml_bnni.treefile"

VERDICTS = ("held", "weakened", "lost", "strengthened", "unsupported in both",
            "not a clade in the main tree")


def _outgroup(root: Node) -> set[str]:
    return {lf.label for lf in leaves(root) if group_of(lf.label) == "RYR"}


def load_pair() -> tuple[Node, Node]:
    """(reported tree, ``--bnni`` tree), both rooted on the RyR outgroup."""
    main = parse(G.read_text(TREE))
    alt = parse(G.read_text(BNNI_TREE))
    main, alt = reroot(main, _outgroup(main)), reroot(alt, _outgroup(alt))
    if {lf.label for lf in leaves(main)} != {lf.label for lf in leaves(alt)}:
        raise ValueError("the two trees do not carry the same tips")
    return main, alt


def largest_pure_clade(root: Node, group: str) -> set[str]:
    """The biggest clade whose every tip's census group is ``group``."""
    best: set[str] = set()
    stack = [root]
    while stack:
        n = stack.pop()
        tips = {lf.label for lf in leaves(n)}
        if all(group_of(t) == group for t in tips):
            best = max(best, tips, key=len)       # pure: nothing below is bigger
        else:
            stack.extend(n.children)
    return best


def claim_sets(root: Node) -> dict[str, set[str]]:
    """Every clade claim of the origin paper, as a tip set, by this project's
    rules, in the order S7's tables list them."""
    tips = {lf.label for lf in leaves(root)}
    core = {g: largest_pure_clade(root, g) for g in PARALOGS}
    clade = {g: {lf.label for lf in leaves(c.node)} for g, c in paralog_clades(root).items()}
    out = {f"{g} core clade": core[g] for g in PARALOGS}
    out |= {f"{g} clade incl. unlabelled tips": clade[g] for g in PARALOGS}
    for a, b in (("ITPR1", "ITPR2"), ("ITPR1", "ITPR3"), ("ITPR2", "ITPR3")):
        out[f"{a} + {b}"] = clade[a] | clade[b]
    out["all three paralogs"] = set().union(*clade.values())
    out["vertebrate ITPRs (incl. unassigned)"] = vertebrate_tips(root)
    out["ITPR family (all non-RyR)"] = tips - _outgroup(root)
    out["RyR outgroup"] = _outgroup(root)
    return out


def strong(s) -> bool:
    return s is not None and s[0] >= _P.value("tree.alrt_min") \
        and s[1] >= _P.value("tree.ufboot_min")


@dataclass
class Claim:
    name: str
    tips: set
    main: tuple | None          # (SH-aLRT, UFBoot) if a clade, else None
    alt: tuple | None
    main_clade: bool
    alt_clade: bool

    @property
    def verdict(self) -> str:
        if not self.main_clade:
            return "not a clade in the main tree"
        if not self.alt_clade:
            return "lost"
        m, a = strong(self.main), strong(self.alt)
        return ("held" if m and a else "weakened" if m else
                "strengthened" if a else "unsupported in both")


def _read(root: Node, tips: set) -> tuple[bool, tuple | None]:
    n = clade_node(root, tips)
    if n is None:
        return False, None
    if n.parent is root:                # the root edge: IQ-TREE labels both sides
        return True, next((support(c) for c in root.children if support(c)), None)
    return True, support(n)


def compare(main: Node, alt: Node) -> list[Claim]:
    """Every claim read from both trees (the sets are rebuilt on the main
    tree, so both are asked the same question)."""
    out = []
    for name, tips in claim_sets(main).items():
        (mc, ms), (ac, as_) = _read(main, tips), _read(alt, tips)
        out.append(Claim(name, tips, ms, as_, mc, ac))
    return out
