"""The Gating view of the Dynamics tab: steady-state P_open against Ca2+.

Three IP3R models and one RyR1 scheme, each measured with one ruler
(``physics.bell``):

* **De Young-Keizer** (Li-Rinzel): the bell of Bezprozvanny et al. (1991),
  IP3 moving both flanks.
* **Mak et al. 1998**: a Hill steady state where IP3 tunes inhibition alone,
  as measured.
* **Park/drive** (Siekmann 2012; Cao 2013): Ca2+ and IP3 act only on the
  switch between a quiet park mode and an active drive mode.
* **The three side by side** at the two IP3 levels of Mak's flank test,
  each bell relative to its own peak, so the test can be seen.
* **RyR1**: Stern 1997's scheme against Murayama 2015's measured bell.

Steady state only (the gating variables at equilibrium at clamped Ca2+):
the single-channel condition.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QVBoxLayout, QWidget

from ..parameters import PARAMETERS as _P
from ..physics import gating as dyk
from ..physics import gating_mak as mk
from ..physics import park_drive as pd
from ..physics import ryr_gating as rg
from ..physics import spark_mg as sm
from .plot_canvas import PALETTE, PlotCanvas

__all__ = ["GatingPanel", "MODELS"]

#: key -> label, in the combo's order.
MODELS = {
    "dyk": "De Young–Keizer (Li–Rinzel)",
    "mak": "Mak et al. 1998 (Hill)",
    "pd": "Park/drive (Siekmann 2012; Cao 2013)",
    "compare": "IP3R: the three models side by side",
    "ryr1": "RyR1: Stern 1997 scheme vs Murayama 2015 bell",
}

#: The IP3R models: (module with open_probability / bell_at, IP3 levels, µM).
_IP3R = {"dyk": (dyk, (0.1, 0.3, 1.0, 10.0)),
         "mak": (mk, (0.01, 0.02, 0.033, 0.1, 10.0)),
         "pd": (pd, (0.03, 0.1, 0.3, 1.0, 10.0))}
_FLANK_KEY = {"dyk": "DYK", "mak": "Mak", "pd": "Park/drive"}

_NOTES = {
    "dyk": "Steady-state P_open = (m∞ n∞ h∞)³ at clamped Ca²⁺ and IP3 "
           "(Li & Rinzel 1994; De Young & Keizer 1992 constants).",
    "mak": "P_open = P_max / [1 + (K_act/c)^H_act + (c/K_inh(IP3))^H_inh] "
           "(Mak, McBride & Foskett 1998, Eqs. 1–2; Xenopus IP3R-1). "
           "Steady state only: no kinetics.",
    "pd": "Six states in two modes, fitted to IP3R-1 single-channel records "
          "(Siekmann et al. 2012). Inside a mode the rates are constant: a "
          "driving channel is open about 70 % of the time, a parked one almost "
          "never. Ca²⁺ and IP3 act only on the switch (Cao et al. 2013's "
          "gating variables at equilibrium).",
    "compare": "Each bell relative to its own peak (P_max differs: DYK "
               "follows IP3, Mak 0.81, park/drive's drive mode 0.70). "
               "Dotted: the half-peak flanks. Mak's test: raising IP3 must "
               "move the inhibitory flank and leave the activating one.",
}


class GatingPanel(QWidget):
    """The steady-state bell under a chosen model, with the flank test that
    tells the IP3R models apart (``gating_mak.compare_flanks``)."""

    def __init__(self):
        super().__init__()
        lay = QVBoxLayout(self)
        row = QHBoxLayout()
        row.addWidget(QLabel("Model"))
        self.model = QComboBox()
        for key, label in MODELS.items():
            self.model.addItem(label, key)
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

    def select(self, key: str) -> None:
        self.model.setCurrentIndex(self.model.findData(key))

    @property
    def model_key(self) -> str:
        return self.model.currentData()

    def draw(self):
        key = self.model_key
        if key == "ryr1":
            return self._draw_ryr()
        if key == "compare":
            return self._draw_compare()
        mod, levels = _IP3R[key]
        ax1, ax2 = self.canvas.reset(1, 2)[0]
        c = np.logspace(-2, 2.5, 400)
        lines = []
        for k, p in enumerate(levels):
            b = mod.bell_at(p)
            ax1.semilogx(c, mod.open_probability(c, p), color=PALETTE[k % 6],
                         label=f"IP3 {p:g} µM")
            ax1.plot([b.c_peak], [b.po_peak], "o", color=PALETTE[k % 6], ms=3)
            if key == "dyk":
                ax2.semilogx(c, dyk.h_inf(c, p), color=PALETTE[k % 6])
            elif key == "pd":
                ax2.semilogx(c, [pd.park_fraction(x, p) for x in c],
                             color=PALETTE[k % 6])
            lines.append(f"IP3 {p:g} µM: peak {b.po_peak:.3f} at {b.c_peak:.2f} µM")
        ax1.set_xlabel("Ca²⁺ (µM)")
        ax1.set_ylabel("P_open (steady state)")
        ax1.set_title("the bell: activation then inhibition")
        self._right_panel(key, ax2)
        self.canvas.legend(ax1, loc="upper left")
        self.canvas.draw_now()
        self.note.setText(_NOTES[key])
        self.text.setText("; ".join(lines) + ". " + self._flank_sentence())

    @staticmethod
    def _right_panel(key: str, ax) -> None:
        if key == "mak":
            p = np.logspace(-2.3, 1, 300)
            ax.loglog(p, mk.k_inh(p), color=PALETTE[0])
            ax.axhline(mk.MakParams().k_act, color="grey", ls=":", lw=1)
            ax.set_xlabel("IP3 (µM)")
            ax.set_ylabel("K_inh (µM)")
            ax.set_title("IP3 tunes K_inh alone (K_act dotted)")
        elif key == "pd":
            ax.set_xlabel("Ca²⁺ (µM)")
            ax.set_ylabel("fraction parked (C4 + O5)")
            ax.set_title("Ca²⁺ and IP3 act on the mode alone")
            ax.set_ylim(0, 1.02)
        else:
            ax.set_xlabel("Ca²⁺ (µM)")
            ax.set_ylabel("h∞ (inhibitory site free)")
            ax.set_title("IP3 relieves Ca²⁺ inhibition")

    def _draw_compare(self):
        """The three IP3R bells at Mak's two IP3 levels, one panel each."""
        lo, hi = _P.value("mak.compare_ip3_low"), _P.value("mak.compare_ip3_high")
        axes = self.canvas.reset(1, 2)[0]
        c = np.logspace(-2, 2.5, 400)
        for ax, p in zip(axes, (lo, hi)):
            for k, key in enumerate(_IP3R):
                mod = _IP3R[key][0]
                b = mod.bell_at(p)
                ax.semilogx(c, mod.open_probability(c, p) / b.po_peak,
                            color=PALETTE[k], label=MODELS[key].split(" (")[0])
                for x in (b.c_half_act, b.c_half_inh):
                    if np.isfinite(x):
                        ax.axvline(x, color=PALETTE[k], ls=":", lw=0.8)
            ax.set_xlabel("Ca²⁺ (µM)")
            ax.set_title(f"IP3 {p:g} µM")
        axes[0].set_ylabel("P_open / own peak")
        self.canvas.legend(axes[0], loc="upper right")    # empty above 10 µM
        self.canvas.draw_now()
        self.note.setText(_NOTES["compare"])
        peaks = "; ".join(
            f"{MODELS[k].split(' (')[0]} {m.bell_at(lo).po_peak:.2g} → "
            f"{m.bell_at(hi).po_peak:.2g}" for k, (m, _) in _IP3R.items())
        self.text.setText(self._flank_sentence() + f" P_open at the peak, IP3 "
                          f"{lo:g} → {hi:g} µM: {peaks}.")

    def _draw_ryr(self):
        """RyR1 has no IP3: one kinetic scheme's bell against one measured
        bell, each normalised to its own peak, flanks on one ruler."""
        ax1, ax2 = self.canvas.reset(1, 2)[0]
        c = np.logspace(-2, 4, 500)
        bells = rg.compare_bells()
        curves = (rg.open_probability(c), rg.murayama_activity(c))
        for k, ((name, b), y) in enumerate(zip(bells.items(), curves)):
            ax1.semilogx(c, y / b.po_peak, color=PALETTE[k], label=name)
            for x in (b.c_half_act, b.c_half_inh):
                ax1.axvline(x, color=PALETTE[k], ls=":", lw=0.8)
        fit = rg.fit_to_bell()
        y = rg.open_probability(c, fit)
        ax1.semilogx(c, y / y.max(), color=PALETTE[0], ls="--",
                     label="Stern, fitted to it")
        mg = rg.with_mg(fit)                   # the fibre's free Mg2+, K_Mg,A measured
        ax1.semilogx(c, rg.open_probability(c, mg) / y.max(), color=PALETTE[0],
                     ls="-.", label=f"fitted, {mg.mg / 1e3:g} mM Mg²⁺")
        ax1.set_xlabel("Ca²⁺ (µM)")
        ax1.set_ylabel("activity / own peak")
        ax1.set_title("RyR1 bells; dotted: half-peak flanks")
        occ = np.array([rg.stationary(x) for x in c])
        for k, s in enumerate(rg.STATES):
            ax2.semilogx(c, occ[:, k], color=PALETTE[k + 2], label=s)
        ax2.set_xlabel("Ca²⁺ (µM)")
        ax2.set_ylabel("stationary occupancy")
        ax2.set_title("Stern 1997: two gates in series")
        self.canvas.legend(ax1, loc="upper left")       # the empty low-Ca2+ side
        self.canvas.legend(ax2, loc="center left")
        self.canvas.draw_now()
        s, m = bells.values()
        self.note.setText(
            "Stern, Pizarro & Ríos 1997: activation by two Ca²⁺ (k_o c², k_o−), "
            "inactivation by one (k_i c, k_i−); the authors did not fit it to "
            "data. Murayama et al. 2015: [³H]ryanodine binding of rabbit RyR1, "
            "A = Amax fA (1 − fI) — an activity index, not P_open, so both are "
            "drawn relative to their own peaks.")
        self.text.setText(
            f"Half-activation {s.c_half_act:.1f} µM (scheme) vs "
            f"{m.c_half_act:.1f} µM (measured); half-inhibition "
            f"{s.c_half_inh:.0f} vs {m.c_half_inh:.0f} µM. The scheme activates "
            f"where RyR1 does but inactivates {m.c_half_inh / s.c_half_inh:.1f}× "
            "too readily, as its authors said of their inactivation site. "
            f"Fitted to both flanks (dashed; off rates moved, on rates kept): "
            f"Ka {fit.k_a:.1f} µM and Ki {fit.k_i:.0f} µM, against "
            f"{rg.SternParams().k_a:.1f} and {rg.SternParams().k_i:.0f}. "
            "In the cleft, sparks under the fitted scheme never end. "
            + self._mg_sentence(fit, mg, y.max()))

    @staticmethod
    def _mg_sentence(fit, mg, peak) -> str:
        """What the fibre's Mg2+ does to the fitted bell, at each K_Mg,A
        reading (``spark_mg.READINGS``)."""
        sel = rg.with_mg(fit, k_mg_a=sm.k_mg_a_reading(fit, "selectivity"))
        mei = rg.with_mg(fit, k_mg_a=sm.k_mg_a_reading(fit, "meissner"))
        b, bs, bm = rg.bell_at(mg), rg.bell_at(sel), rg.bell_at(mei)
        return (f"Under {mg.mg:g} µM free Mg²⁺ (dash-dot, drawn relative "
                f"to the Mg²⁺-free peak; Laver 2004's two sites) the peak falls to {100 * b.po_peak / peak:.0f} % of "
                f"its Mg²⁺-free height and half-activation moves to "
                f"{b.c_half_act:.0f} µM with K_Mg,A {mg.k_mg_a:.0f} µM "
                f"(measured), or {100 * bs.po_peak / peak:.0f} % and "
                f"{bs.c_half_act:.0f} µM with {sel.k_mg_a:.0f} µM (Laver's "
                f"selectivity on the fitted Ka), or {100 * bm.po_peak / peak:.0f} % "
                f"and {bm.c_half_act:.0f} µM with {mei.k_mg_a:.0f} µM (Meissner "
                "1997, same assay, Na⁺ competing). Under Mg²⁺ triggered sparks "
                "do end (Puffs: Triggered sparks vs Mg²⁺).")

    @staticmethod
    def _flank_sentence() -> str:
        r = mk.compare_flanks()
        lo, hi = _P.value("mak.compare_ip3_low"), _P.value("mak.compare_ip3_high")
        parts = [f"{a:.2f}× and {i:.1f}× in {MODELS[k].split(' (')[0]}"
                 for k, (a, i) in ((k, r[_FLANK_KEY[k]]) for k in _IP3R)]
        return (f"IP3 {lo:g} → {hi:g} µM moves half-activation and "
                f"half-inhibition " + "; ".join(parts) + ". Mak 1998 measured the "
                "activating flank fixed while IP3 moved the inhibitory one.")
