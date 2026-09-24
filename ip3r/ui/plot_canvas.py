"""A matplotlib canvas styled to sit beside the dark viewport.

Every 2-D plot in the application — pore profiles, bell curves, Ca2+
traces, puff rasters, check exhibits — is drawn on one of these, so they
share fonts, colours and the way a redraw is requested.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("QtAgg")

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

__all__ = ["PlotCanvas", "PALETTE"]

#: Categorical series colours (readable on the dark background).
PALETTE = ["#5b9ef2", "#f28c4d", "#72cc80", "#cc80e6", "#f2cc4d", "#66d9d9"]

_BG = "#12151c"
_FG = "#d7dbe3"


class PlotCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5.0, height=3.2, rows=1, cols=1):
        self.figure = Figure(figsize=(width, height), dpi=100, facecolor=_BG)
        super().__init__(self.figure)
        self.setParent(parent)
        self.reset(rows, cols)

    def reset(self, rows: int = 1, cols: int = 1):
        """Clear and return a fresh grid of styled axes."""
        self.figure.clear()
        axes = self.figure.subplots(rows, cols, squeeze=False)
        for ax in axes.ravel():
            self.style(ax)
        self.axes = axes
        return axes[0, 0] if rows == cols == 1 else axes

    def span_row(self, axes_row):
        """Replace a row of axes by one styled axis spanning it."""
        gs = axes_row[0].get_subplotspec()
        spec = gs.get_gridspec()[gs.rowspan.start, :]
        for ax in axes_row:
            ax.remove()
        return self.style(self.figure.add_subplot(spec))

    @staticmethod
    def style(ax):
        ax.set_facecolor(_BG)
        for sp in ax.spines.values():
            sp.set_color("#4a505c")
        ax.tick_params(colors=_FG, labelsize=8)
        ax.xaxis.label.set_color(_FG)
        ax.yaxis.label.set_color(_FG)
        ax.title.set_color(_FG)
        ax.title.set_fontsize(9)
        return ax

    def legend(self, ax, **kw):
        leg = ax.legend(fontsize=7, frameon=False, **kw)
        for t in leg.get_texts():
            t.set_color(_FG)
        return leg

    def draw_now(self):
        self.figure.tight_layout()
        self.draw_idle()
