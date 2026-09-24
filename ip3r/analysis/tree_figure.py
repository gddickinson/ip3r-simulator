"""Paper 2's tree drawn as a rectangular phylogram on a matplotlib ``Axes``.

Everything annotated is computed by :mod:`.tree` on the tree being drawn:
each paralogue's whole clade is boxed and labelled with its size and
support, the ryanodine-receptor outgroup is boxed, tips are marked by what
their record says (a named paralogue, an unnamed vertebrate, a cyclostome,
anything else), and a black dot marks every node clearing the registered
joint support bar. Used by the Tree tab and by the P2 exhibits.
"""

from __future__ import annotations

from dataclasses import dataclass

from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle

from .newick import Node, leaves
from .tree import (bipartitions, cyclostome_clades, group_of, is_cyclostome, is_supported,
                   outgroup_clade, paralog_clades, support)

__all__ = ["CLADE_COLORS", "Layout", "layout", "draw_tree"]

#: One hue per boxed clade (readable on the dark background; not the
#: conservation ramp, so a tree beside a painted structure is not confused).
CLADE_COLORS = {"ITPR1": "#5b9ef2", "ITPR2": "#f28c4d", "ITPR3": "#72cc80",
                "RYR": "#cc80e6"}
_CYCLO = "#f2cc4d"
_GREY = "#8a8f99"
_FG = "#d7dbe3"


@dataclass
class Layout:
    x: dict      # id(node) -> distance from the root
    y: dict      # id(node) -> row (tips 0..n-1)
    tips: list   # leaves in drawing order


def _ladderize(n: Node) -> list[Node]:
    return sorted(n.children, key=lambda c: len(leaves(c)))


def layout(root: Node) -> Layout:
    """Tip rows in ladderized order; x is the summed branch length."""
    x, y, tips = {}, {}, []

    def walk(n: Node, depth: float) -> float:
        x[id(n)] = depth
        if n.is_leaf:
            y[id(n)] = float(len(tips))
            tips.append(n)
        else:
            ys = [walk(c, depth + (c.length or 0.0)) for c in _ladderize(n)]
            y[id(n)] = (ys[0] + ys[-1]) / 2
        return y[id(n)]

    walk(root, 0.0)
    return Layout(x, y, tips)


def _edges(n: Node, lay: Layout, out: list) -> None:
    for c in n.children:
        out.append([(lay.x[id(n)], lay.y[id(c)]), (lay.x[id(c)], lay.y[id(c)])])
        _edges(c, lay, out)
    if n.children:
        ys = [lay.y[id(c)] for c in n.children]
        out.append([(lay.x[id(n)], min(ys)), (lay.x[id(n)], max(ys))])


def _box(ax, lay: Layout, node: Node, colour: str, text: str, right: float) -> None:
    ys = [lay.y[id(t)] for t in leaves(node)]
    x0 = lay.x[id(node)]
    ax.add_patch(Rectangle((x0, min(ys) - 0.5), right - x0, max(ys) - min(ys) + 1,
                           facecolor=colour, alpha=0.14, edgecolor=colour, lw=0.8))
    ax.text(right, (min(ys) + max(ys)) / 2, " " + text, color=colour, fontsize=7,
            va="center", ha="left")


def draw_tree(ax, root: Node, labels: bool = False, supports: bool = True) -> dict:
    """Draw ``root`` on ``ax``; returns the annotation counts it drew."""
    lay = layout(root)
    segs: list = []
    _edges(root, lay, segs)
    ax.add_collection(LineCollection(segs, colors=_GREY, linewidths=0.6))
    right = max(lay.x.values()) * 1.02

    clades = paralog_clades(root)
    in_clade = {t: g for g, c in clades.items() for t in c.tips}
    for g, c in clades.items():
        _box(ax, lay, c.node, CLADE_COLORS[g], f"{g}  {len(c.tips)}  {_sup(c.node)}", right)
    og = outgroup_clade(root)
    _box(ax, lay, og.node, CLADE_COLORS["RYR"], f"RyR  {len(og.tips)}", right)

    groups: dict[str, list] = {"named": [], "unnamed": [], "cyclostome": [], "other": []}
    for t in lay.tips:
        g = group_of(t.label)
        kind = ("cyclostome" if is_cyclostome(t.label) else
                "named" if g in CLADE_COLORS else
                "unnamed" if g == "vertebrate_basal" else "other")
        groups[kind].append(t)
    pts = lambda ts: ([lay.x[id(t)] for t in ts], [lay.y[id(t)] for t in ts])  # noqa: E731
    named = groups["named"]
    ax.scatter(*pts(named), s=9, zorder=3, linewidths=0,
               c=[CLADE_COLORS[group_of(t.label)] for t in named])
    un = groups["unnamed"]
    ax.scatter(*pts(un), s=11, zorder=3, facecolors="none", linewidths=0.8,
               edgecolors=[CLADE_COLORS.get(in_clade.get(t.label), _GREY) for t in un])
    ax.scatter(*pts(groups["cyclostome"]), s=22, zorder=4, marker="D",
               c=_CYCLO, edgecolors="#12151c", linewidths=0.5, label="hagfish / lamprey")
    ax.scatter(*pts(groups["other"]), s=4, zorder=3, c=_GREY, linewidths=0)

    cc = cyclostome_clades(root)
    for c in cc:                   # its own support, and where it joins the rest
        for n in (c.node, c.node.parent):
            if n is not None and _sup(n):
                ax.text(lay.x[id(n)], lay.y[id(n)], f"{_sup(n)} ", color=_CYCLO,
                        fontsize=6, ha="right", va="bottom")
    n_sup = 0
    if supports:
        ok = [n for n in _internal(root) if is_supported(n)]
        n_sup = sum(map(is_supported, bipartitions(root)))
        ax.scatter([lay.x[id(n)] for n in ok], [lay.y[id(n)] for n in ok], s=5,
                   c=_FG, zorder=4, linewidths=0)
    if labels:
        for t in lay.tips:
            ax.text(lay.x[id(t)], lay.y[id(t)], " " + t.label, fontsize=3.2, va="center",
                    clip_on=True,
                    color=CLADE_COLORS.get(in_clade.get(t.label), _GREY))
    ax.set_xlim(-0.02 * right, right * (1.45 if labels else 1.22))
    ax.set_ylim(len(lay.tips), -1)
    ax.set_yticks([])
    ax.set_xlabel("substitutions per site")
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    return {"tips": len(lay.tips), "clades": {g: len(c.tips) for g, c in clades.items()},
            "cyclostome_clades": [len(c.tips) for c in cc], "supported_nodes": n_sup}


def _internal(root: Node) -> list[Node]:
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        if n.children and n is not root:
            out.append(n)
        stack.extend(n.children)
    return out


def _sup(n: Node) -> str:
    s = support(n)
    return "" if s is None else f"{s[0]:g}/{s[1]:g}"
