"""Paper 1's range drawn, and the exhibits of the range checks.

:func:`draw_range` draws one row per eukaryotic clade, grouped by supergroup.
A grey track is the proteomes swept and the coloured bar is the fraction with
an IP3 receptor call, on a fixed 0–1 scale. The prokaryotes are collapsed to
one row per domain. A cross marks a clade where the absence was confirmed in
controlled genome assemblies (a smaller one where it held for classes inside
the clade). Used by the Range tab and by the P1 exhibits, which draw only
from their check's own outcome.
"""

from __future__ import annotations

import numpy as np
from matplotlib.patches import Patch

from . import range_genomes as RG
from .range_table import SUPERGROUPS

__all__ = ["GROUP_COLORS", "draw_range", "draw_presence", "draw_relaxed",
           "draw_absences", "draw_copies", "draw_chase", "draw_genomes",
           "VERDICT_COLORS", "STATUS_COLORS", "COPIES_MAX"]

_TEXT, _TRACK, _RED, _GREY = "#d7dbe3", "#3a3e47", "#e05a5a", "#8a8f99"
#: One fixed colour per eukaryotic supergroup (prokaryotes never have a bar).
GROUP_COLORS = dict(zip(SUPERGROUPS, ("#5b9ef2", "#f28c4d", "#72cc80", "#cc80e6",
                                      "#f2cc4d", "#66d9d9", "#e6a0b4", _GREY, _GREY)))
_PROK = ("Archaea", "Bacteria")
_MISSING = "#6b7079"
#: Fixed, never auto-ranged: the most copies S23 counts is 18 (Macrostomum).
COPIES_MAX = 20
VERDICT_COLORS = {"controlled_cross_kingdom": "#5b9ef2", "controlled_by_target": "#66d9d9",
                  "controlled_partial": "#f28c4d", "no_control_bait": "#e05a5a"}
STATUS_COLORS = {"found_annotated": "#72cc80", "found_no_annotation": "#a8e0a0",
                 "found_unannotated": "#a8e0a0", "fragment_only": "#f2cc4d",
                 "assembly_gap": "#f28c4d", "tblastn_trace": "#cc80e6",
                 "no_locus": "#20232a"}
_SPANS = {True: "#72cc80", False: "#f28c4d"}


def _rows(clades, min_n: int, collapse: bool) -> list[tuple]:
    """(label, supergroup, n, present, clade names) in drawing order."""
    rows, prok = [], {}
    for c, g, n, k in clades:
        if collapse and g in _PROK:
            p = prok.setdefault(g, [0, 0, []])
            p[0] += n
            p[1] += k
            p[2].append(c)
        elif n >= min_n:
            rows.append((c, g, n, k, (c,)))
    for g, (n, k, names) in prok.items():
        rows.append((f"{g} ({len(names)} phyla)", g, n, k, tuple(names)))
    order = {s: i for i, s in enumerate(SUPERGROUPS)}
    return sorted(rows, key=lambda r: (order.get(r[1], len(order)), -r[2], r[0]))


def draw_range(ax, clades, absences: dict | None = None, min_n: int = 1,
               collapse: bool = True, legend: bool = True) -> dict:
    """``clades``: (clade, supergroup, proteomes, with a call) tuples;
    ``absences``: clade → (targets holding at the clade itself, targets
    holding for classes inside it). Returns the drawn rows (top first)."""
    rows = _rows(clades, min_n, collapse)
    y = np.arange(len(rows))
    ax.barh(y, [1.0] * len(rows), color=_TRACK, height=0.8)
    ax.barh(y, [k / n for _, _, n, k, _ in rows], height=0.8,
            color=[GROUP_COLORS.get(g, _GREY) for _, g, _, _, _ in rows])
    for i, (_, _, n, k, _) in enumerate(rows):
        ax.text(1.02, i, f"{k:,}/{n:,}", va="center", fontsize=5.5, color=_TEXT)
    absences = absences or {}
    for i, (label, *_rest) in enumerate(rows):
        own, inside = absences.get(label, (0, 0))
        if own:
            ax.plot(0.025, i, marker="x", color=_RED, ms=5, mew=1.5)
        elif inside:
            ax.plot(0.025, i, marker="x", color=_RED, ms=3, mew=0.8)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=6 if len(rows) <= 70 else 4)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of swept proteomes with an IP3R call (fixed 0–1)")
    for i in range(1, len(rows)):
        if rows[i][1] != rows[i - 1][1]:
            ax.axhline(i - 0.5, color=_TEXT, lw=0.5, alpha=0.5)
    if legend:
        groups = [g for g in SUPERGROUPS if g not in _PROK and any(r[1] == g for r in rows)]
        handles = [Patch(color=GROUP_COLORS[g], label=g) for g in groups]
        handles.append(Patch(color=_TRACK, label="proteomes swept"))
        if absences:
            handles.append(ax.plot([], [], "x", color=_RED, label="absence held in "
                                   "controlled genomes (small: in a class inside)")[0])
        ax.legend(handles=handles, fontsize=6, frameon=False, labelcolor=_TEXT,
                  loc="upper left", bbox_to_anchor=(1.16, 1.0))
    return {"rows": rows, "clades": len(rows),
            "with_call": sum(r[3] > 0 for r in rows),
            "proteomes": sum(r[2] for r in rows), "present": sum(r[3] for r in rows)}


# ------------------------------------------------------------- exhibits

def draw_presence(ax, d) -> None:
    """P1.presence_range / kingdom_absences / absence_targets: the clade
    bars, clades of fewer than two proteomes left out for legibility."""
    draw_range(ax, d["clades"], min_n=2, legend=False)


def draw_relaxed(ax, d) -> None:
    """P1.relaxed_controls: substantial matches per lineage × profile, the
    claimed absences above their in-kingdom controls."""
    lins = d["lineages"]
    profiles = ("PF08709", "itpr", "ryr", "PF01365", "PF08454", "PF02815")
    m = np.array([[d["table"].get(lin, {}).get(p, (0, 0))[1] for p in profiles]
                  for lin in lins], float)
    ax.imshow(np.log10(m + 1), cmap="cividis", vmin=0, vmax=4, aspect="auto")
    for i in range(len(lins)):
        for j in range(len(profiles)):
            v = int(m[i, j])
            ax.text(j, i, f"{v:,}", ha="center", va="center", fontsize=7,
                    color="#111" if v >= 300 else _TEXT)
    ax.set_xticks(range(len(profiles)), profiles, fontsize=7)
    ax.set_yticks(range(len(lins)), lins, fontsize=7)
    ax.set_xlabel("substantial matches (E ≤ 1e-5, ≥ half the model); colour log10, fixed 0–4")


def draw_absences(ax, d) -> None:
    """P1.absences: each target's genomes, controlled and with a gene."""
    rows = sorted(d["absences"], key=lambda r: (-r[1], r[0]))
    y = np.arange(len(rows))
    ax.barh(y, [g for _, g, _, _ in rows], color=_TRACK, label="genomes swept")
    ax.barh(y, [c for _, _, c, _ in rows], height=0.5, color="#5b9ef2",
            label="controlled")
    ax.barh(y, [f for *_, f in rows], height=0.5, color=_RED, label="with a gene")
    ax.set_yticks(y, [r[0] for r in rows], fontsize=5)
    ax.invert_yaxis()
    ax.set_xlabel("non-vertebrate genome assemblies")
    ax.legend(fontsize=6, frameon=False, labelcolor=_TEXT, loc="lower right")


def draw_copies(ax, d) -> None:
    """P1.copy_number: genomes by complete gene models carried."""
    c = np.asarray(d["copies"])
    ks = np.arange(0, c.max() + 1)
    ax.bar(ks, [int((c == k).sum()) for k in ks], color="#5b9ef2")
    ax.set_yscale("symlog", linthresh=5)
    ax.set_xlabel("complete IP3R gene models in the genome")
    ax.set_ylabel("genomes (symlog)")


def draw_chase(ax, d) -> None:
    """P1.record_chase: length against the best identity outside the
    kingdom, with the fragment floor and the contaminant bar."""
    style = {"real_gene": ("o", "real gene"), "fragment": ("^", "fragment")}
    for v, (mk, lab) in style.items():
        for k, col in (("Viridiplantae", "#72cc80"), ("Fungi", "#f28c4d")):
            p = [(n, i) for n, i, o, kk in d["points"] if o == v and kk == k]
            if p:
                ax.scatter(*zip(*p), marker=mk, s=12, color=col, label=f"{k}: {lab}")
    ax.axvline(d["floor"], color=_TEXT, lw=0.8, ls="--")
    ax.axhline(d["contaminant"], color=_RED, lw=0.8, ls="--")
    ax.set_ylim(0, 100)
    ax.set_xlabel("record length (aa); dashed: family floor")
    ax.set_ylabel("identity to nearest relative outside the kingdom (%)")
    ax.text(50, d["contaminant"] - 5, "contaminant bar", color=_RED, fontsize=6)
    ax.legend(fontsize=6, frameon=False, labelcolor=_TEXT, loc="center right")


def draw_genomes(ax, rows, title: str = "") -> dict:
    """S23's genomes of one clade: control verdict, copy-ledger status and
    contiguity against the genome's own bar as cells, complete gene models
    as bars on a fixed 0–:data:`COPIES_MAX` scale. Unknown values are grey."""
    from matplotlib.colors import to_rgb
    n = len(rows)
    img = np.empty((n, 3, 3))
    for i, r in enumerate(rows):
        span = _SPANS[r.spans_gene] if np.isfinite(r.contig_n50) else _MISSING
        for j, c in enumerate((VERDICT_COLORS.get(r.verdict, _MISSING),
                               STATUS_COLORS.get(r.status, _MISSING), span)):
            img[i, j] = to_rgb(c)
    ax.imshow(img, aspect="auto", interpolation="nearest",
              extent=(-0.5, 2.5, n - 0.5, -0.5))
    ax.set_xticks(range(3), ["control", "copy\nledger", "contig\n≥ bar"], fontsize=7)
    ax.set_yticks(range(n), [r.organism for r in rows] if n <= 60 else [""] * n,
                  fontsize=6 if n <= 30 else 4)
    bars = ax.inset_axes([1.04, 0.0, 0.3, 1.0], sharey=ax)
    bars.set_facecolor(ax.get_facecolor())          # the theme's, not matplotlib's white
    for side in bars.spines.values():
        side.set_color(ax.spines["left"].get_edgecolor())
    bars.tick_params(colors=ax.xaxis.label.get_color())
    bars.xaxis.label.set_color(ax.xaxis.label.get_color())
    bars.barh(range(n), [r.copies for r in rows], color=_GREY)
    for i, r in enumerate(rows):
        if r.copies:
            bars.text(r.copies + 0.3, i, str(r.copies), va="center", fontsize=6, color=_TEXT)
    bars.set_xlim(0, COPIES_MAX)
    bars.tick_params(labelleft=False, labelsize=6)
    bars.set_xlabel("complete gene models", fontsize=7)
    seen_v = {r.verdict for r in rows}
    seen_s = {r.status for r in rows}
    handles = ([Patch(color=VERDICT_COLORS[v], label=v.replace("_", " "))
                for v in RG.VERDICT_ORDER if v in seen_v]
               + [Patch(color=STATUS_COLORS[s], label=s.replace("_", " "))
                  for s in RG.STATUS_ORDER if s in seen_s]
               + [Patch(color=_SPANS[True], label="contig N50 ≥ own bar"),
                  Patch(color=_SPANS[False], label="contig N50 < own bar")])
    ax.legend(handles=handles, fontsize=6, frameon=False, labelcolor=_TEXT,
              loc="upper left", bbox_to_anchor=(0.0, -0.08), ncol=2)
    if title:
        ax.set_title(title, fontsize=9)
    return {"genomes": n, "controlled": sum(r.controlled for r in rows),
            "with_gene": sum(r.copies > 0 for r in rows),
            "copies": sum(r.copies for r in rows)}
