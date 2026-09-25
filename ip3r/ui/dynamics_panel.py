"""The dynamics tab: what the channel does with IP3 and Ca2+.

Three views of IP3R gating, from one molecule to one cell:

* **Gating** (``gating_panel``) — steady-state open probability against
  Ca2+ at several IP3 levels, under De Young-Keizer, Mak et al. (1998) or
  park/drive, or all three side by side: the bell of Bezprozvanny et al.
  (1991), and IP3 moving its right flank (inhibition). Mak and park/drive
  leave the left flank (activation) where it is; DYK does not.
* **Oscillations** — the closed-cell model at a chosen IP3; a scan finds the
  IP3 window in which Ca2+ oscillates with no oscillating input.
* **Puffs** (``puffs_panel``) — a stochastic cluster of De Young-Keizer or
  park/drive receptors, with and without the Ca2+ coupling that lets one
  opening recruit the next.

Every constant comes from the parameter registry; the plots say which model
drew them.
"""

from __future__ import annotations

from PyQt6.QtWidgets import (QDoubleSpinBox, QFormLayout, QHBoxLayout, QLabel,
                             QPushButton, QTabWidget, QVBoxLayout, QWidget)

from ..physics.calcium import oscillation_metrics, oscillation_window, simulate
from .gating_panel import GatingPanel
from .plot_canvas import PALETTE, PlotCanvas
from .puffs_panel import PuffsPanel
from .workers import run_async

__all__ = ["DynamicsPanel"]


def _spin(value, lo, hi, step, suffix="", decimals=2):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setSingleStep(step)
    s.setDecimals(decimals)
    s.setValue(value)
    s.setSuffix(suffix)
    return s


class _Oscillation(QWidget):
    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.p = _spin(0.5, 0.0, 5.0, 0.05, " µM")
        self.t_end = _spin(100.0, 10.0, 1000.0, 10.0, " s", 0)
        form.addRow("IP3", self.p)
        form.addRow("Duration", self.t_end)
        lay.addLayout(form)
        row = QHBoxLayout()
        run = QPushButton("Simulate")
        run.clicked.connect(self.run)
        scan = QPushButton("Find oscillation window")
        scan.clicked.connect(self.scan)
        row.addWidget(run)
        row.addWidget(scan)
        lay.addLayout(row)
        self.canvas = PlotCanvas(self, height=3.4)
        lay.addWidget(self.canvas, 1)
        self.text = QLabel()
        self.text.setWordWrap(True)
        lay.addWidget(self.text)

    def run(self):
        tr = simulate(self.p.value(), t_end=self.t_end.value())
        m = oscillation_metrics(tr)
        ax = self.canvas.reset()
        ax.plot(tr.t, tr.c, color=PALETTE[0], label="cytosolic Ca²⁺")
        ax.plot(tr.t, tr.h, color=PALETTE[2], lw=0.8, label="h (inhibition gate)")
        ax.set_xlabel("time (s)")
        ax.set_ylabel("µM  /  fraction")
        ax.set_title(f"closed-cell Li-Rinzel model, IP3 = {self.p.value():g} µM")
        self.canvas.legend(ax, loc="upper right")
        self.canvas.draw_now()
        self.text.setText(
            f"Oscillating: period {m['period']:.1f} s, amplitude {m['amplitude']:.2f} µM."
            if m["oscillates"] else
            f"Not oscillating: settles at {m['mean']:.3f} µM.")

    def scan(self):
        self.text.setText("Scanning IP3 0.20-1.00 µM (81 simulations)…")
        run_async(oscillation_window, on_done=self._scanned,
                  on_error=lambda e: self.text.setText(e))

    def _scanned(self, out):
        lo, hi, rows = out
        ax = self.canvas.reset()
        p = [r["p"] for r in rows]
        ax.plot(p, [r["amplitude"] if r["oscillates"] else 0 for r in rows],
                color=PALETTE[1], label="sustained amplitude")
        ax.plot(p, [r["mean"] for r in rows], color=PALETTE[0], lw=0.8, label="mean Ca²⁺")
        ax.axvspan(lo, hi, color=PALETTE[1], alpha=0.12, lw=0)
        ax.set_xlabel("IP3 (µM)")
        ax.set_ylabel("Ca²⁺ (µM)")
        ax.set_title("oscillation window (measured by simulation)")
        self.canvas.legend(ax, loc="upper left")
        self.canvas.draw_now()
        self.text.setText(f"Ca²⁺ oscillates for IP3 between {lo:.2f} and {hi:.2f} µM "
                          "(grid 0.01 µM; damped spirals excluded).")


class DynamicsPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        tabs = QTabWidget()
        self.gating = GatingPanel()
        tabs.addTab(self.gating, "Gating")
        tabs.addTab(_Oscillation(), "Oscillations")
        self.puffs = PuffsPanel()
        tabs.addTab(self.puffs, "Puffs")
        lay.addWidget(tabs)
