"""The dynamics tab: what the channel does with IP3 and Ca2+.

Three views of the De Young-Keizer / Li-Rinzel model, from one molecule to
one cell:

* **Gating** — steady-state open probability against Ca2+ at several IP3
  levels: the bell of Bezprozvanny et al. (1991), and IP3 moving its right
  flank (inhibition) while the left (activation) stays put.
* **Oscillations** — the closed-cell model at a chosen IP3; a scan finds the
  IP3 window in which Ca2+ oscillates with no oscillating input.
* **Puffs** — a stochastic cluster of receptors, with and without the Ca2+
  coupling that lets one opening recruit the next.

Every constant comes from the parameter registry; the plots say which model
drew them.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QDoubleSpinBox, QFormLayout, QHBoxLayout, QLabel,
                             QPushButton, QTabWidget, QVBoxLayout, QWidget)

from ..physics.calcium import oscillation_metrics, oscillation_window, simulate
from ..physics.gating import bell_peak, h_inf, open_probability
from ..physics.puffs import PuffParams, fano, simulate_cluster
from .plot_canvas import PALETTE, PlotCanvas
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


class _Gating(QWidget):
    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        note = QLabel("Steady-state P_open = (m∞ n∞ h∞)³ at clamped Ca²⁺ and "
                      "IP3 (Li & Rinzel 1994; De Young & Keizer 1992 constants).")
        note.setWordWrap(True)
        lay.addWidget(note)
        self.canvas = PlotCanvas(self, height=3.4, cols=2)
        lay.addWidget(self.canvas, 1)
        self.text = QLabel()
        self.text.setWordWrap(True)
        lay.addWidget(self.text)
        self.draw()

    def draw(self):
        ax1, ax2 = self.canvas.reset(1, 2)[0]
        c = np.logspace(-2, 2, 400)
        lines = []
        for k, p in enumerate((0.1, 0.3, 1.0, 10.0)):
            ax1.semilogx(c, open_probability(c, p), color=PALETTE[k], label=f"IP3 {p:g} µM")
            cp, po = bell_peak(p)
            ax1.plot([cp], [po], "o", color=PALETTE[k], ms=3)
            ax2.semilogx(c, h_inf(c, p), color=PALETTE[k])
            lines.append(f"IP3 {p:g} µM: peak P_open {po:.3f} at {cp:.2f} µM Ca²⁺")
        ax1.set_xlabel("Ca²⁺ (µM)")
        ax1.set_ylabel("P_open (steady state)")
        ax1.set_title("the bell: activation then inhibition")
        ax2.set_xlabel("Ca²⁺ (µM)")
        ax2.set_ylabel("h∞ (inhibitory site free)")
        ax2.set_title("IP3 relieves Ca²⁺ inhibition")
        self.canvas.legend(ax1, loc="upper left")
        self.canvas.draw_now()
        self.text.setText("; ".join(lines) + ". More IP3 raises the peak and "
                          "moves it to higher Ca²⁺: IP3 relieves Ca²⁺ inhibition. "
                          "In this model the activating flank shifts too, less "
                          "than the inhibitory one; experiment (Mak et al. 1998) "
                          "finds inhibition alone is tuned.")


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


class _Puffs(QWidget):
    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        form = QFormLayout()
        pp = PuffParams()
        self.p = _spin(0.2, 0.0, 5.0, 0.05, " µM")
        self.n = _spin(pp.n_channels, 1, 200, 1, " channels", 0)
        self.coupling = _spin(pp.ca_per_open, 0.0, 20.0, 0.1, " µM / open channel")
        self.duration = _spin(10.0, 1.0, 60.0, 1.0, " s", 0)
        for label, w in (("IP3", self.p), ("Cluster size", self.n),
                         ("Ca²⁺ coupling", self.coupling), ("Duration", self.duration)):
            form.addRow(label, w)
        lay.addLayout(form)
        run = QPushButton("Simulate cluster (coupled vs uncoupled, same seed)")
        run.clicked.connect(self.run)
        lay.addWidget(run)
        self.canvas = PlotCanvas(self, height=3.4, rows=2)
        lay.addWidget(self.canvas, 1)
        self.text = QLabel()
        self.text.setWordWrap(True)
        self.text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.text)

    def _params(self, coupling):
        pp = PuffParams()
        pp.n_channels, pp.ca_per_open = self.n.value(), coupling
        return pp

    def run(self):
        self.text.setText("Simulating…")
        p, d = self.p.value(), self.duration.value()
        run_async(lambda: [simulate_cluster(p, d, 0, self._params(c))
                           for c in (self.coupling.value(), 0.0)],
                  on_done=self._show, on_error=self.text.setText)

    def _show(self, traces):
        axes = self.canvas.reset(2, 1)
        for ax, tr, label in zip(axes[:, 0], traces, ("coupled", "uncoupled")):
            ax.step(tr.t, tr.n_open, where="post", lw=0.6,
                    color=PALETTE[0 if label == "coupled" else 1])
            ax.set_ylabel("open")
            ax.set_title(f"{label}: Fano {fano(tr):.2f}, max {tr.n_open.max()} "
                         f"of {int(tr.params.n_channels)} open at once")
        axes[-1, 0].set_xlabel("time (s)")
        self.canvas.draw_now()
        self.text.setText(
            "Fano factor = variance / mean of the number open at once: ~1 when "
            "channels gate independently, above 1 when one opening recruits "
            "others through Ca²⁺-induced Ca²⁺ release. The DYK constants give a "
            "high resting activity, so puffs here are modest (see ROADMAP).")


class DynamicsPanel(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        tabs = QTabWidget()
        tabs.addTab(_Gating(), "Gating")
        tabs.addTab(_Oscillation(), "Oscillations")
        tabs.addTab(_Puffs(), "Puffs")
        lay.addWidget(tabs)
