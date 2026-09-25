"""The microdomain cluster in the Puffs panel (Round 7.3, from Round 7.2).

Park/drive receptors in Cao 2013's point domain with Cao 2014's pools and
fluo-4 (``physics.puffs_domain``). Puffs are read from F/F0, as the
experiments read them (``physics.puff_stats``), and the inter-puff intervals
are fitted by Thurley's refractory density (Cao 2013 Eq. 14) beside an
exponential with the same mean. The likelihood ratio says whether the
refractory period is significant in this run.

``DomainControls`` holds the controls that only this model has: the store
clamp, the release multiple and the run length. ``draw`` puts a trace on
the Puffs panel's canvas and returns the text and the result for the smoke
test.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtWidgets import QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox, QPushButton

from ..parameters import PARAMETERS as _P
from ..physics.puff_stats import fluorescence_events, puff_stats, thurley_pdf
from ..physics.puffs_domain import DomainPuffParams, simulate_cluster_domain
from .plot_canvas import PALETTE
from .view_state import Seeded, set_combo, set_spin

__all__ = ["DomainControls", "simulate", "draw", "CLAMP_LABELS"]

#: ``microdomain.CLAMPS`` in words, in the combo's order.
CLAMP_LABELS = {"store": "store held (Cao 2013)", "none": "store free",
                "bath": "store and cytosol held"}


class DomainControls(QGroupBox):
    """Clamp, release multiple, duration and the run button."""

    def __init__(self):
        super().__init__("In the microdomain (fluo-4, puffs from F/F0)")
        form = QFormLayout(self)
        self.clamp = QComboBox()
        for key, label in CLAMP_LABELS.items():
            self.clamp.addItem(label, key)
        self.scale = QDoubleSpinBox()
        self.scale.setRange(0.1, 10.0)
        self.scale.setSingleStep(0.5)
        self.scale.setValue(1.0)                 # Cao 2014's release, unscaled
        self.scale.setSuffix(" × release")
        self.duration = QDoubleSpinBox()
        self.duration.setRange(5.0, 3600.0)
        self.duration.setDecimals(0)
        self.duration.setSingleStep(30.0)
        self.duration.setSuffix(" s")
        self.duration.setObjectName("microdomain duration")
        self.seeded = Seeded()
        self.seeded.seed(self.duration, lambda: _P.value("domain.gui_duration"))
        self._tooltips()
        self.run_btn = QPushButton("Simulate in the microdomain (F/F0, IPIs)")
        form.addRow("Store", self.clamp)
        form.addRow("Release", self.scale)
        form.addRow("Duration", self.duration)
        form.addRow(self.run_btn)

    def _tooltips(self) -> None:
        self.scale.setToolTip(
            "Multiplies each open receptor's release k_ipr. At 1× the mean blip "
            "is dF/F0 0.6-0.7; Cao 2013's 1.6 is matched at "
            f"{_P.value('domain.blip_scale'):g}× (domain.blip_scale).")

    def follow(self) -> None:
        """A parameter changed: the duration follows unless typed over."""
        self._tooltips()
        self.seeded.follow()

    def restore(self, d: dict) -> list[str]:
        notes: list[str] = []
        if "clamp" in d:
            set_combo(self.clamp, d["clamp"], "store clamp", notes)
        for key, w, what in (("scale", self.scale, "release multiple"),
                             ("duration", self.duration, "microdomain duration")):
            if key in d:
                set_spin(w, d[key], what, notes)
        return notes

    def spec(self) -> dict:
        return {"clamp": self.clamp.currentData(), "scale": self.scale.value(),
                "duration": self.duration.value()}


def simulate(p: float, n: int, clamp: str, scale: float, duration: float,
             seed: int = 0):
    """One run and its statistics (called on a worker)."""
    dp = DomainPuffParams(n_channels=n, clamp=clamp)
    dp.domain = dp.domain.scaled(scale)
    tr = simulate_cluster_domain(p, duration, seed, dp)
    return tr, puff_stats(tr)


def draw(canvas, tr, s, spec: dict) -> tuple[str, dict]:
    """F/F0, number open and the IPI density on ``canvas``."""
    axes = canvas.reset(3, 1)[:, 0]
    cut = _P.value("domain.puff_threshold") * s.blip_mean
    puffs = [e for e in fluorescence_events(tr) if e["peak_df"] > cut]
    ax = axes[0]
    ax.plot(tr.t, tr.f_peak, lw=0.6, color=PALETTE[2])
    ax.axhline(1.0 + cut, color="grey", ls=":", lw=0.8)
    ax.plot([e["start"] for e in puffs], [1.0 + e["peak_df"] for e in puffs],
            "v", color=PALETTE[1], ms=3)
    ax.set_ylabel("F/F0")
    ax.set_title(f"fluo-4: {s.n_puffs} puffs (▼; dotted: "
                 f"{_P.value('domain.puff_threshold'):g} mean blips)")
    ax = axes[1]
    ax.step(tr.t, tr.peaks, where="post", lw=0.6, color=PALETTE[0])
    ax.set_ylabel("open")
    ax.set_xlabel("time (s)")
    _draw_ipis(axes[2], s)
    canvas.draw_now()
    result = {"model": "park-drive-domain", "stats": s, "spec": spec,
              "open_fraction": float(tr.n_open.mean() / tr.params.n_channels)}
    return _text(tr, s, spec), result


def _draw_ipis(ax, s) -> None:
    ax.set_xlabel("inter-puff interval (s)")
    ax.set_ylabel("density (1/s)")
    if len(s.ipis) < 5:
        ax.text(0.5, 0.5, f"{len(s.ipis)} intervals: too few to fit",
                transform=ax.transAxes, ha="center", color="grey")
        return
    ax.hist(s.ipis, bins=min(30, max(5, len(s.ipis) // 3)), density=True,
            color=PALETTE[0], alpha=0.5, label=f"{len(s.ipis)} intervals")
    t = np.linspace(0, s.ipis.max(), 300)
    ax.plot(t, thurley_pdf(t, s.lam, s.xi), color=PALETTE[1],
            label=f"Thurley: λ {s.lam:.2f}, ξ {s.xi:.2f} /s")
    rate = 1.0 / s.ipis.mean()
    ax.plot(t, rate * np.exp(-rate * t), color="#d7dbe3", ls="--", lw=1,
            label="exponential, same mean")
    ax.legend(fontsize=7, frameon=False, labelcolor="#d7dbe3")


def _text(tr, s, spec: dict) -> str:
    fit = ("too few intervals to fit Eq. 14" if len(s.ipis) < 5 else
           f"Thurley's Eq. 14 fits λ {s.lam:.2f} and ξ {s.xi:.2f} /s; the "
           f"refractory period is {'' if s.refractory_lr > 3.84 else 'not '}"
           f"significant (LR {s.refractory_lr:.1f} against 3.84)")
    return (f"Park/drive in the microdomain, {CLAMP_LABELS[spec['clamp']]}, "
            f"release ×{spec['scale']:g}, {spec['duration']:g} s. Rest "
            f"{tr.rest.c:.3f} µM, store {tr.rest.cs:.0f} µM; "
            f"{100 * tr.n_open.mean() / tr.params.n_channels:.1f} % open on average. "
            f"{s.n_blips} blips (mean dF/F0 {s.blip_mean:.2f}; Cao 2013: 1.6), "
            f"{s.n_puffs} puffs ({s.rate:.2f}/s, mean dF/F0 {s.amp_mean:.2f}, "
            f"{s.open_mean:.1f} open at the peak). Interval CV {s.cv:.2f} "
            f"(Cao: 0.65-0.95); {fit}. A single run: the CLI's "
            "'microdomain --scan ah42' scans h42 recovery.")
