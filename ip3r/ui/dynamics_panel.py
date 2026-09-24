"""The dynamics tab: what the channel does with IP3 and Ca2+.

Three views of IP3R gating, from one molecule to one cell:

* **Gating** — steady-state open probability against Ca2+ at several IP3
  levels, under De Young-Keizer or Mak et al. (1998): the bell of
  Bezprozvanny et al. (1991), and IP3 moving its right flank (inhibition).
  Only the Mak model leaves the left flank (activation) where it is.
* **Oscillations** — the closed-cell model at a chosen IP3; a scan finds the
  IP3 window in which Ca2+ oscillates with no oscillating input.
* **Puffs** (``puffs_panel``) — a stochastic cluster of De Young-Keizer or
  park/drive receptors, with and without the Ca2+ coupling that lets one
  opening recruit the next.

Every constant comes from the parameter registry; the plots say which model
drew them.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout, QHBoxLayout, QLabel,
                             QPushButton, QTabWidget, QVBoxLayout, QWidget)

from ..physics.calcium import oscillation_metrics, oscillation_window, simulate
from ..parameters import PARAMETERS as _P
from ..physics import gating_mak as mk
from ..physics.gating import bell_at, h_inf, open_probability
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


class _Gating(QWidget):
    """The steady-state bell under either gating model, with the flank test
    that tells them apart (``gating_mak.compare_flanks``)."""

    MODELS = ("De Young–Keizer (Li–Rinzel)", "Mak et al. 1998 (Hill)")
    LEVELS = ((0.1, 0.3, 1.0, 10.0), (0.01, 0.02, 0.033, 0.1, 10.0))

    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel("Model"))
        self.model = QComboBox()
        self.model.addItems(self.MODELS)
        self.model.currentIndexChanged.connect(self.draw)
        row.addWidget(self.model)
        row.addStretch(1)
        lay.addLayout(row)
        self.note = QLabel()
        self.note.setWordWrap(True)
        lay.addWidget(self.note)
        self.canvas = PlotCanvas(self, height=3.4, cols=2)
        lay.addWidget(self.canvas, 1)
        self.text = QLabel()
        self.text.setWordWrap(True)
        lay.addWidget(self.text)
        self.draw()

    def draw(self):
        mak = self.model.currentIndex() == 1
        ax1, ax2 = self.canvas.reset(1, 2)[0]
        c = np.logspace(-2, 2.5, 400)
        lines = []
        for k, p in enumerate(self.LEVELS[mak]):
            po = mk.open_probability(c, p) if mak else open_probability(c, p)
            b = mk.bell_at(p) if mak else bell_at(p)
            ax1.semilogx(c, po, color=PALETTE[k], label=f"IP3 {p:g} µM")
            ax1.plot([b.c_peak], [b.po_peak], "o", color=PALETTE[k], ms=3)
            if not mak:
                ax2.semilogx(c, h_inf(c, p), color=PALETTE[k])
            lines.append(f"IP3 {p:g} µM: peak {b.po_peak:.3f} at "
                         f"{b.c_peak:.2f} µM")
        ax1.set_xlabel("Ca²⁺ (µM)")
        ax1.set_ylabel("P_open (steady state)")
        ax1.set_title("the bell: activation then inhibition")
        if mak:
            p = np.logspace(-2.3, 1, 300)
            ax2.loglog(p, mk.k_inh(p), color=PALETTE[0])
            ax2.axhline(mk.MakParams().k_act, color="grey", ls=":", lw=1)
            ax2.set_xlabel("IP3 (µM)")
            ax2.set_ylabel("K_inh (µM)")
            ax2.set_title("IP3 tunes K_inh alone (K_act dotted)")
        else:
            ax2.set_xlabel("Ca²⁺ (µM)")
            ax2.set_ylabel("h∞ (inhibitory site free)")
            ax2.set_title("IP3 relieves Ca²⁺ inhibition")
        self.canvas.legend(ax1, loc="upper left")
        self.canvas.draw_now()
        self.note.setText(
            "P_open = P_max / [1 + (K_act/c)^H_act + (c/K_inh(IP3))^H_inh] "
            "(Mak, McBride & Foskett 1998, Eqs. 1–2; Xenopus IP3R-1). "
            "Steady state only: no kinetics." if mak else
            "Steady-state P_open = (m∞ n∞ h∞)³ at clamped Ca²⁺ and IP3 "
            "(Li & Rinzel 1994; De Young & Keizer 1992 constants).")
        self.text.setText("; ".join(lines) + ". " + self._flank_sentence())

    @staticmethod
    def _flank_sentence() -> str:
        r = mk.compare_flanks()
        lo, hi = _P.value("mak.compare_ip3_low"), _P.value("mak.compare_ip3_high")
        (da, di), (ma, mi) = r["DYK"], r["Mak"]
        return (f"IP3 {lo:g} → {hi:g} µM moves the half-inhibition point "
                f"{di:.1f}× and half-activation {da:.2f}× in De Young–Keizer; "
                f"{mi:.1f}× and {ma:.3f}× in Mak 1998, where IP3 tunes "
                f"inhibition alone, as measured.")


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
        self.gating = _Gating()
        tabs.addTab(self.gating, "Gating")
        tabs.addTab(_Oscillation(), "Oscillations")
        self.puffs = PuffsPanel()
        tabs.addTab(self.puffs, "Puffs")
        lay.addWidget(tabs)
