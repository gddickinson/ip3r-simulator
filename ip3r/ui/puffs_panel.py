"""The Puffs view of the Dynamics tab: a stochastic cluster of receptors,
with and without the Ca2+ coupling that lets one opening recruit the next.

Two receptors can fill the same cluster: De Young-Keizer subunits and the
park/drive scheme (Siekmann; Cao et al. 2013). RyR1 sparks come in two
forms: one mean-field cluster Ca2+, or each channel's own Ca2+ in the
junctional cleft, or in the cleft with the gating fitted to Murayama's
measured bell. All are read with one ruler
(``physics.puff_compare``): traces, the Fano factor, and the event-size
distribution, where a clean blip/puff split shows as a valley. "Scan
coupling" runs both receptors over the registered coupling range.

The RyR1 receptors take free Mg2+ (``ryr_gating.with_mg``) at either
reading of the activation-site affinity (``spark_mg.READINGS``). Under the
fibre's Mg2+ no spark starts by itself, so "Triggered sparks vs Mg2+"
(cleft receptors) opens the available channels at t = 0 and times the array
to shut, over ``spark_mg.mg_values``.

The park/drive receptor can also be run in Cao 2013's microdomain with
fluo-4 (``puffs_domain_view``): F/F0, the number open, and the inter-puff
intervals against Thurley's refractory density.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout, QGridLayout,
                             QLabel, QPushButton, QVBoxLayout, QWidget)

from ..parameters import PARAMETERS as _P
from ..physics import puff_compare as pc
from ..physics import spark_mg as sm
from . import puffs_domain_view as dv
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
    "ryr1": "A RyR1 cluster (Stern 1997's unfitted scheme; IP3 is ignored). "
            "The coupling is the free-diffusion Ca²⁺ one open channel "
            "delivers at 30 nm (an unbuffered upper estimate). Sparks recruit "
            "most of the cluster, then settle at about 5 open (the mean-field "
            "point where 30 × P_open(cluster Ca²⁺) equals the number open) "
            "until they close by chance, so they last ~100 ms against a "
            "measured release of about 6 ms (frog).",
    "ryr1-cleft": "The same RyR1s in Stern 1997's junctional cleft: two rows, "
                  "C and V channels alternating, each C channel seeing its "
                  "own Ca²⁺ (its own release about 50 µM, a diagonal "
                  "neighbour's about 11 µM; the coupling above is the "
                  "nearest-neighbour value). Sparks now end by local "
                  "inactivation after about 20 ms, not ~130 ms, still about "
                  "3× the measured release (frog).",
    "ryr1-cleft-fit": "The cleft array with Stern's Ka and Ki fitted to "
                      "Murayama 2015's measured bell (Ka 4.9, Ki 249 µM, "
                      "against 7.1 and 10). Once a spark starts it never "
                      "ends: at the tens of µM a channel sees in the cleft, "
                      "a Ki of 249 µM inactivates too few channels. No "
                      "inactivation rate changes this (spark-termination "
                      "in the CLI). What terminates real sparks is missing "
                      "from a Ca²⁺-only scheme fitted without Mg²⁺.",
}


#: Receptors whose channels each see their own Ca2+ (the trigger needs one).
_CLEFT = (pc.SPARK_CLEFT, pc.SPARK_FIT)


def _mg_args(m: str, mg: float, reading: str) -> dict:
    """Free Mg2+ and K_Mg,A for receptor ``m`` (called on the worker: the
    fitted scheme is solved here); none for an IP3R receptor."""
    if m not in pc.SPARKS:
        return {}
    return {"mg": mg, "k_mg_a": sm.k_mg_a_reading(pc.params_for(m).gating, reading)}


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
        for m in pc.ALL_MODELS:
            self.model.addItem(pc.MODEL_LABELS[m], m)
        pp = pc.params_for(pc.MODELS[0])
        self.p = _spin(0.2, 0.0, 5.0, 0.05, " µM")
        self.n = _spin(pp.n_channels, 1, 200, 1, " channels", 0)
        self.coupling = _spin(pp.ca_per_open, 0.0, 20.0, 0.05, " µM / open channel")
        self.duration = _spin(10.0, 1.0, 60.0, 1.0, " s", 0)
        self.mg = _spin(0.0, 0.0, 10 * _P.value("ryr.mg_free"), 100.0, " µM free", 0)
        self.mg.setToolTip(f"Cytosolic free Mg²⁺ (RyR1 only); the fibre's is "
                           f"{_P.value('ryr.mg_free'):g} µM")
        self.reading = QComboBox()
        for key, label in sm.READINGS.items():
            self.reading.addItem(label, key)
        for label, w in (("Receptor", self.model), ("IP3", self.p),
                         ("Cluster size", self.n), ("Ca²⁺ coupling", self.coupling),
                         ("Duration", self.duration), ("Mg²⁺", self.mg),
                         ("Mg²⁺ at the activation site", self.reading)):
            form.addRow(label, w)
        self.model.currentIndexChanged.connect(self._model_changed)
        self.form = form
        lay.addLayout(form)
        row = QGridLayout()               # two columns: the dock stays narrow
        run = QPushButton("Simulate cluster (coupled vs uncoupled, same seed)")
        run.clicked.connect(self.run)
        scan = QPushButton("Scan coupling (both receptors)")
        scan.clicked.connect(self.scan)
        self.mg_scan_btn = QPushButton("Triggered sparks vs Mg²⁺")
        self.mg_scan_btn.clicked.connect(self.scan_mg)
        row.addWidget(run, 0, 0, 1, 2)
        row.addWidget(scan, 1, 0)
        row.addWidget(self.mg_scan_btn, 1, 1)
        lay.addLayout(row)
        self.domain = dv.DomainControls()
        self.domain.run_btn.clicked.connect(self.run_domain)
        lay.addWidget(self.domain)
        self.canvas = PlotCanvas(self, height=4.2, rows=3)
        lay.addWidget(self.canvas, 1)
        self.text = QLabel()
        self.text.setWordWrap(True)
        self.text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(self.text)
        self.result: dict | None = None
        self._model_changed()

    @property
    def model_key(self) -> str:
        return self.model.currentData()

    def _model_changed(self):
        pp = pc.params_for(self.model_key)
        self.coupling.setValue(pp.ca_per_open)
        self.n.setValue(pp.n_channels)
        ryr = self.model_key in pc.SPARKS
        # A control that cannot act on this receptor is hidden, not greyed:
        # RyR ignores IP3, and the IP3R models have no Mg2+ sites.
        self.p.setEnabled(not ryr)
        self.mg.setEnabled(ryr)
        self.reading.setEnabled(ryr)
        self.form.setRowVisible(self.p, not ryr)
        self.form.setRowVisible(self.mg, ryr)
        self.form.setRowVisible(self.reading, ryr)
        self.mg_scan_btn.setVisible(ryr)
        self.mg_scan_btn.setEnabled(self.model_key in _CLEFT)
        self.domain.setVisible(self.model_key == "park-drive")

    def run_domain(self):
        spec = self.domain.spec()
        p, n = self.p.value(), int(self.n.value())
        self.text.setText(f"Simulating {spec['duration']:g} s of the cluster in "
                          "its microdomain (about 0.3 s per simulated second)…")
        run_async(lambda: dv.simulate(p, n, **spec),
                  on_done=lambda out: self._show_domain(out, spec),
                  on_error=self.text.setText)

    def _show_domain(self, out, spec):
        text, self.result = dv.draw(self.canvas, *out, spec)
        self.text.setText(text)

    def run(self):
        self.text.setText("Simulating…")
        m, p, d = self.model_key, self.p.value(), self.duration.value()
        n, c = self.n.value(), self.coupling.value()
        mg, reading = self.mg.value(), self.reading.currentData()
        run_async(lambda: [pc.simulate(m, p, d, 0, pc.params_for(m, cc, n, **kw))
                           for kw in [_mg_args(m, mg, reading)]
                           for cc in (c, 0.0)],
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
        if any(r["n_events"] for r in stats):   # a silent cluster has nothing to log
            ax.set_yscale("log")
        else:
            ax.text(0.5, 0.5, "no events", transform=ax.transAxes, ha="center",
                    color="grey")
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
            f"half the cluster (dotted line). {_NOTE[self.model_key]}"
            + self._mg_sentence(traces[0]))

    def _mg_sentence(self, tr) -> str:
        sp = getattr(tr.params, "gating", None)
        if sp is None or sp.mg == 0:
            return ""
        return (f" Under {sp.mg:g} µM free Mg²⁺ (K_Mg,A {sp.k_mg_a:.0f} µM) "
                f"half-activation moves from {sp.k_a:.1f} to {sp.k_a_eff:.1f} µM"
                " and the inactivation gate sees the Mg²⁺ as Ca²⁺.")

    def scan_mg(self):
        m, reading = self.model_key, self.reading.currentData()
        self.text.setText(f"Triggering {pc.MODEL_LABELS[m]} over free Mg²⁺ "
                          f"({_P.value('spark.trigger_trials'):g} trials per point)…")

        def work():
            base = pc.params_for(m).gating
            k = sm.k_mg_a_reading(base, reading)
            return m, f"{k:.0f} µM ({sm.READINGS[reading]})", sm.mg_scan(base, k)
        run_async(work, on_done=self._show_mg, on_error=self.text.setText)

    def _show_mg(self, out):
        m, k_label, rows = out
        window = 1e3 * _P.value("spark.trigger_window")
        axes = self.canvas.reset(2, 1)[:, 0]
        mg = np.array([r.mg for r in rows])
        x = np.where(mg > 0, mg, np.nan)      # Mg2+ 0 is in the text, not on a log axis
        ended = np.array([r.ended == r.trials for r in rows])
        dur = np.array([r.duration_ms for r in rows])
        ax = axes[0]
        ax.plot(x[ended], dur[ended], "o-", color=PALETTE[0], label="every spark shut")
        ax.plot(x[~ended], np.full((~ended).sum(), window), "^", mfc="none",
                color=PALETTE[1], label="some never shut (window)")
        ax.set_yscale("log")
        ax.set_ylabel("median spark (ms)")
        ax.set_title("triggered: available channels opened at t = 0", fontsize=9)
        ax = axes[1]
        ax.plot(x, [r.opened for r in rows], "o-", color=PALETTE[2], label="opened")
        ax.plot(x, [r.inactivated_start for r in rows], "s-", color=PALETTE[3],
                label="inactivated at rest")
        ax.set_ylabel("channels")
        ax.set_xlabel("free Mg²⁺ (µM)")
        for a in axes:
            a.set_xscale("log")
            self.canvas.legend(a)
        self.canvas.draw_now()
        self.result = {"model": m, "mg_scan": rows}
        zero, top = rows[0], rows[-1]
        first = next((r for r in rows if r.ended == r.trials), None)
        self.text.setText(
            f"{pc.MODEL_LABELS[m]}, K_Mg,A {k_label}. "
            f"Without Mg²⁺, {zero.ended}/{zero.trials} triggered sparks shut "
            f"within {window:.0f} ms. "
            + (f"Every spark shuts from {first.mg:.3g} µM ({first.duration_ms:.1f} ms); "
               if first else "No Mg²⁺ in the scan shuts every spark; ")
            + f"at {top.mg:.0f} µM, {top.opened:.0f} channels open, "
            f"{top.inactivated_start:.0f} are inactivated before the trigger, and "
            f"the spark lasts {top.duration_ms:.1f} ms (frog release: about 6 ms). "
            "Competition at the activation site ends sparks by induction "
            "decay; Mg²⁺ at the inactivation site mostly removes channels "
            "before the spark (spark-mg in the CLI dissects the two).")

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
