"""The Puffs view of the Dynamics tab: a stochastic cluster of receptors,
with and without the Ca2+ coupling that lets one opening recruit the next.

Two receptors can fill the same cluster: De Young-Keizer subunits and the
park/drive scheme (Siekmann; Cao et al. 2013). Both are read with one ruler
(``physics.puff_compare``): traces, the Fano factor, and the event-size
distribution, where a clean blip/puff split shows as a valley. "Scan
coupling" runs both receptors over the registered coupling range.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout, QHBoxLayout,
                             QLabel, QPushButton, QVBoxLayout, QWidget)

from ..physics import puff_compare as pc
from .plot_canvas import PALETTE, PlotCanvas
from .workers import run_async

__all__ = ["PuffsPanel"]

_NOTE = {
    "dyk": "Every subunit's activating site is already partly occupied at "
           "rest, so single channels open often and coupling adds only modest "
           "clustering: event sizes fall away with no gap.",
    "park-drive": "At rest nearly every receptor is parked, so the cluster "
                  "is quiet apart from brief park-mode flickers. Once one "
                  "receptor drives, its Ca²⁺ can pull the others into drive "
                  "mode. Event sizes split into blips and whole-cluster puffs "
                  "with a valley between.",
}


def _spin(value, lo, hi, step, suffix="", decimals=2):
    s = QDoubleSpinBox()
    s.setRange(lo, hi)
    s.setSingleStep(step)
    s.setDecimals(decimals)
    s.setValue(value)
    s.setSuffix(suffix)
    return s


class PuffsPanel(QWidget):
    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        form = QFormLayout()
        self.model = QComboBox()
        for m in pc.MODELS:
            self.model.addItem(pc.MODEL_LABELS[m], m)
        pp = pc.params_for(pc.MODELS[0])
        self.p = _spin(0.2, 0.0, 5.0, 0.05, " µM")
        self.n = _spin(pp.n_channels, 1, 200, 1, " channels", 0)
        self.coupling = _spin(pp.ca_per_open, 0.0, 20.0, 0.05, " µM / open channel")
        self.duration = _spin(10.0, 1.0, 60.0, 1.0, " s", 0)
        for label, w in (("Receptor", self.model), ("IP3", self.p),
                         ("Cluster size", self.n), ("Ca²⁺ coupling", self.coupling),
                         ("Duration", self.duration)):
            form.addRow(label, w)
        self.model.currentIndexChanged.connect(self._model_changed)
        lay.addLayout(form)
        row = QHBoxLayout()
        run = QPushButton("Simulate cluster (coupled vs uncoupled, same seed)")
        run.clicked.connect(self.run)
        scan = QPushButton("Scan coupling (both receptors)")
        scan.clicked.connect(self.scan)
        row.addWidget(run)
        row.addWidget(scan)
        lay.addLayout(row)
        self.canvas = PlotCanvas(self, height=4.2, rows=3)
        lay.addWidget(self.canvas, 1)
        self.text = QLabel()
        self.text.setWordWrap(True)
        self.text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.text)
        self.result: dict | None = None

    @property
    def model_key(self) -> str:
        return self.model.currentData()

    def _model_changed(self):
        self.coupling.setValue(pc.params_for(self.model_key).ca_per_open)

    def run(self):
        self.text.setText("Simulating…")
        m, p, d = self.model_key, self.p.value(), self.duration.value()
        n = self.n.value()
        run_async(lambda: [pc.simulate(m, p, d, 0, pc.params_for(m, c, n))
                           for c in (self.coupling.value(), 0.0)],
                  on_done=self._show, on_error=self.text.setText)

    def _show(self, traces):
        axes = self.canvas.reset(3, 1)[:, 0]
        stats = [pc.recruitment(tr) for tr in traces]
        for ax, tr, r, label, col in zip(axes[:2], traces, stats,
                                         ("coupled", "uncoupled"), PALETTE):
            ax.step(tr.t, tr.peaks, where="post", lw=0.6, color=col)
            ax.set_ylabel("open")
            ev = "event reaches" if r["large"] == 1 else "events reach"
            ax.set_title(f"{label}: Fano {r['fano']:.2f}, {r['large']} {ev} "
                         f"{r['large_at']} of {int(tr.params.n_channels)}",
                         fontsize=9)
        axes[1].set_xlabel("time (s)")
        ax = axes[2]
        top = max(int(tr.params.n_channels) for tr in traces)
        for tr, label, col in zip(traces, ("coupled", "uncoupled"), PALETTE):
            sizes = pc.event_sizes(tr)
            counts = np.bincount(sizes, minlength=top + 1)[1:]
            ax.step(np.arange(1, top + 1), np.where(counts > 0, counts, np.nan),
                    where="mid", color=col, label=label)
        ax.axvline(stats[0]["large_at"] - 0.5, color="grey", ls=":", lw=0.8)
        ax.set_yscale("log")
        ax.set_xlabel("event size (most channels open at once)")
        ax.set_ylabel("events")
        self.canvas.legend(ax)
        self.canvas.draw_now()
        c, u = stats
        self.result = {"model": self.model_key, "coupled": c, "uncoupled": u}
        self.text.setText(
            f"{pc.MODEL_LABELS[self.model_key]}. Uncoupled, {100 * u['open_fraction']:.1f} % "
            f"of channels are open at a time. Coupled: {c['blips']} blips, "
            f"{c['multi']} multi-channel events, {c['large']} of which reach "
            f"half the cluster (dotted line). {_NOTE[self.model_key]}")

    def scan(self):
        self.text.setText("Scanning both receptors over the coupling range…")
        p, d = self.p.value(), self.duration.value()
        run_async(lambda: pc.coupling_scan(p, d, 0), on_done=self._show_scan,
                  on_error=self.text.setText)

    def _show_scan(self, rows):
        axes = self.canvas.reset(2, 1)[:, 0]
        for (m, scan), col in zip(rows.items(), PALETTE):
            x = np.array([r["coupling"] for r in scan])
            axes[0].plot(x[1:], [r["fano"] for r in scan[1:]], "o-", color=col,
                         label=pc.MODEL_LABELS[m])
            axes[1].plot(x[1:], [r["large_per_s"] for r in scan[1:]], "o-", color=col)
        axes[0].axhline(1.0, color="grey", ls=":", lw=0.8)
        axes[0].set_ylabel("Fano factor")
        axes[1].set_ylabel("half-cluster events / s")
        for ax in axes:
            ax.set_xscale("log")
        axes[1].set_xlabel("Ca²⁺ coupling (µM per open channel)")
        self.canvas.legend(axes[0])
        self.canvas.draw_now()
        self.result = {"scan": rows}
        best = {m: max(scan, key=lambda r: r["large_per_s"]) for m, scan in rows.items()}
        self.text.setText("Same cluster, same seed, same couplings. The most "
                          "half-cluster events per second: " + "; ".join(
                              f"{pc.MODEL_LABELS[m]} {b['large_per_s']:.2f} at "
                              f"{b['coupling']:.2f} µM" for m, b in best.items()) + ".")
