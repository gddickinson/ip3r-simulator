"""The genome × paralog grid drawn, and the exhibits of the grid checks.

:func:`draw_grid` paints one layer of a :class:`~.genome_grid.GenomeGrid`
as a heat-map: a genome per row, the three paralogue cells and the RyR
control per column, and a contig-N50 strip on a fixed log scale beside
them. Colours are fixed per category; a cell the table says nothing about
is grey. Used by the GUI's Genomes tab and by the P3/P4 exhibits.
"""

from __future__ import annotations

import math

import numpy as np
from matplotlib import colormaps
from matplotlib.colors import LogNorm, to_rgb
from matplotlib.patches import Patch

from . import genome_grid as GG

__all__ = ["LAYER_STYLE", "LAYER_TITLES", "MISSING", "N50_RANGE", "draw_grid",
           "draw_misses", "draw_logistic", "draw_recovery"]

MISSING = "#6b7079"
N50_RANGE = (1e3, 1e8)          # fixed, never auto-ranged (bp)
_TEXT = "#d7dbe3"
_BLUE, _GREEN, _YELLOW = "#5b9ef2", "#72cc80", "#f2cc4d"
_ORANGE, _PURPLE, _RED, _NONE = "#f28c4d", "#cc80e6", "#e05a5a", "#20232a"

#: layer → ((value, colour, legend label), ...), in legend order.
LAYER_STYLE = {
    "search": (("found", _GREEN, "found"), ("assembly_gap", _YELLOW, "assembly gap"),
               ("fragment", _ORANGE, "fragment"), ("trace", _PURPLE, "translated-search trace"),
               ("absent", _RED, "absent")),
    "miss": (("found", _BLUE, "found"), ("missed", _RED, "missed (gene known present)"),
             ("not a control", _NONE, "no known gene")),
    "state": (("present_single_locus", _GREEN, "single locus"),
              ("present_truncated", _YELLOW, "truncated"),
              ("present_partial", _ORANGE, "partial"),
              ("present_fragmented", _PURPLE, "fragmented"),
              ("paralog_unassignable", _RED, "paralogue unassignable"),
              ("absent", "#ffffff", "absent")),
    "recovery": ((GG.REACHABLE, _BLUE, "reachable from a protein record"),
                 (GG.NO_PROTEOME, _ORANGE, "no reference proteome"),
                 (GG.FRAGMENTS_ONLY, _YELLOW, "only fragmentary records"),
                 (GG.OTHER_PARALOG, _PURPLE, "no record resolves to this paralogue"),
                 (GG.NO_RECORD, _RED, "no family record for the species"),
                 (GG.NO_GENE, _NONE, "no gene")),
}

LAYER_TITLES = {"search": "What the sweep found", "miss": "Missed genes (method control)",
                "state": "S15a evidence state", "recovery": "Protein-record recovery"}


def _n50_rgb(values) -> np.ndarray:
    cmap, norm = colormaps["cividis"], LogNorm(*N50_RANGE, clip=True)
    out = np.empty((len(values), 3))
    for i, v in enumerate(values):
        out[i] = to_rgb(MISSING) if not math.isfinite(v) else cmap(norm(v))[:3]
    return out


def _bar_row(genomes) -> int | None:
    """Index of the first genome below the bar, if the rows are in N50
    order and the bar falls inside them."""
    above = [g.above_bar for g in genomes]
    if not any(above) or all(above):
        return None
    k = above.index(False)
    return k if all(above[:k]) and not any(above[k:]) else None


def draw_grid(ax, grid, layer: str, legend: bool = True) -> dict:
    """Paint ``layer`` of ``grid`` in its current row order; returns counts."""
    style = LAYER_STYLE[layer]
    colour = {v: to_rgb(c) for v, c, _ in style}
    m = grid.layers[layer]
    img = np.empty((grid.n, len(grid.cells) + 1, 3))
    img[:, 0] = _n50_rgb([g.contig_n50 for g in grid.genomes])
    for j in range(len(grid.cells)):
        img[:, j + 1] = [colour.get(v, to_rgb(MISSING)) for v in m[:, j]]
    ax.imshow(img, aspect="auto", interpolation="nearest",
              extent=(-0.5, len(grid.cells) + 0.5, grid.n - 0.5, -0.5))
    ax.axvline(0.5, color=_NONE, lw=2)
    ax.set_xticks(range(len(grid.cells) + 1),
                  ["contig\nN50"] + [c if c != "RYR" else "RyR\n(control)" for c in grid.cells],
                  fontsize=7)
    k = _bar_row(grid.genomes)
    if k is not None:
        ax.axhline(k - 0.5, color=_TEXT, lw=1.0, ls="--")
        ax.text(len(grid.cells) + 0.45, k - 0.8, "contiguity bar", color=_TEXT,
                fontsize=6, ha="right", va="bottom")
    _row_labels(ax, grid)
    counts = grid.counts(layer)
    if legend:
        handles = [Patch(color=c, label=f"{lab} ({counts.get(v, 0)})")
                   for v, c, lab in style if counts.get(v, 0)]
        if "" in set(m.ravel()):
            handles.append(Patch(color=MISSING, label="not in this table"))
        handles.append(Patch(color=colormaps["cividis"](1.0),
                             label="N50 strip: dark 1 kb → yellow 100 Mb (log)"))
        ax.legend(handles=handles, fontsize=6, frameon=False, labelcolor=_TEXT,
                  loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return {"genomes": grid.n, "counts": counts, "bar_row": k,
            "above_bar": sum(g.above_bar for g in grid.genomes)}


def _row_labels(ax, grid) -> None:
    if grid.n <= 60:
        ax.set_yticks(range(grid.n), [g.organism for g in grid.genomes], fontsize=6)
        return
    classes = [g.vclass for g in grid.genomes]
    blocks, start = [], 0
    for i in range(1, grid.n + 1):
        if i == grid.n or classes[i] != classes[start]:
            blocks.append((start, i, classes[start]))
            start = i
    if len(blocks) <= 16:                   # contiguous classes: label the blocks
        big = [(a, b, c) for a, b, c in blocks if b - a >= 4 and c]
        ax.set_yticks([(a + b - 1) / 2 for a, b, _ in big], [c for *_, c in big], fontsize=6)
        for a, _, _ in blocks[1:]:
            ax.axhline(a - 0.5, color=_NONE, lw=0.6)
    else:
        ax.set_yticks([])
        ax.set_ylabel(f"{grid.n} genomes", fontsize=8)


# ------------------------------------------------------------- exhibits

def draw_misses(ax, d) -> None:
    """P3.miss_by_contiguity: the miss layer in N50 order, from the
    outcome's own arrays."""
    genomes = [GG.Genome("", "", "", v) for v in d["n50"]]
    lay = np.array([["" if v is None else ("missed" if v else "found")
                     for v in (d["missed"][c][i] for c in GG.CELLS)]
                    for i in range(len(genomes))], dtype=object)
    lay[lay == ""] = "not a control"
    draw_grid(ax, GG.order(GG.GenomeGrid(genomes, {"miss": lay}), "contig N50"), "miss")


def draw_logistic(ax, d) -> None:
    """P3.contiguity_tests: found (1) or missed (0) against log10 N50, the
    fitted logistic curve of each control series, and the bar."""
    colours = {"itpr_present": _BLUE, "ryr_sister": _ORANGE}
    rng = np.random.default_rng(0)                 # jitter only, for legibility
    xs = np.linspace(math.log10(N50_RANGE[0]), math.log10(N50_RANGE[1]) + 0.5, 200)
    for s, f in d["fits"].items():
        x, y = np.asarray(f["x"]), np.asarray(f["y"])
        off = 0.04 if s == "ryr_sister" else -0.04
        ax.scatter(x, y + off + rng.uniform(-0.02, 0.02, len(y)), s=3, alpha=0.4,
                   color=colours[s])
        ax.plot(xs, 1 / (1 + np.exp(-(f["b0"] + f["b1"] * xs))), color=colours[s],
                lw=1.4, label=f"{s}: OR per 10× {math.exp(f['b1']):.2f}")
    ax.axvline(math.log10(d["bar"]), color=_TEXT, lw=0.8, ls="--")
    ax.set_xlabel("log10 contig N50 (bp); dashed: contiguity bar")
    ax.set_ylabel("P(gene found)")
    ax.legend(fontsize=7, frameon=False, labelcolor=_TEXT, loc="lower right")


def draw_recovery(ax, d) -> None:
    """P4.recovery_channels: each cell's genes split by recovery channel."""
    per = d["per_cell"]
    cells = list(per)
    left = np.zeros(len(cells))
    for v, c, lab in LAYER_STYLE["recovery"]:
        if v == GG.NO_GENE:
            continue
        w = np.array([per[k].get(v, 0) for k in cells], float)
        ax.barh(cells, w, left=left, color=c, label=lab)
        left += w
    ax.invert_yaxis()
    ax.set_xlabel("genes demonstrated in the genome")
    ax.legend(fontsize=6, frameon=False, labelcolor=_TEXT, loc="upper left",
              bbox_to_anchor=(0.0, -0.18), ncol=2)
