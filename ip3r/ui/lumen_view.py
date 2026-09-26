"""The Channel panel's lumen box: its controls and its plot (Rounds 7.10,
7.12).

Controls: draw the lumen; the wall charge's placement (none = the neutral
pore of Round 7.10, one of Round 7.11's closures, or Round 7.13's
dielectric one) and whether salt bridges are paired; what the surface is coloured by. The solve itself is
:class:`~ip3r.ui.lumen_controller.LumenController`'s.

The plot sets the 3-D reading beside the 1-D model's on S0's window: the
lumen's area; where the (K+) drop falls, neutral and charged; and, with a
charge, the wall potential, 3-D plane means against the 1-D Donnan
potential of the same charge on the inscribed circle.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QHBoxLayout, QLabel,
                             QVBoxLayout, QWidget)

from ..physics.lumen_charge import CLOSURE_LABELS, LUMEN_CLOSURES
from ..render.lumen_mesh import COLOURINGS
from .plot_canvas import PALETTE
from .view_state import set_check, set_combo

__all__ = ["LumenControls", "draw_lumen", "NO_CHARGE"]

#: The wall-charge choice meaning "the neutral pore".
NO_CHARGE = "none"


class LumenControls(QWidget):
    toggled = pyqtSignal(bool)
    #: The wall charge's placement or pairing changed (re-solve).
    charge_changed = pyqtSignal()
    #: Only the colouring changed (recolour, no solve).
    colour_changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self.show = QCheckBox("Draw the lumen (3-D, coloured by the potential) "
                              "and plot where the voltage falls")
        self.show.setToolTip(
            "The ion-accessible volume of Round 7.6's 3-D solve, K+ (Cl- not "
            "drawn). Solved on the deposit, so it is hidden on a morph or mode "
            "frame.")
        self.show.toggled.connect(self._toggled)
        lay.addWidget(self.show)
        row = QHBoxLayout()
        row.addWidget(QLabel("Wall charge"))
        self.charge = QComboBox()
        self.charge.addItem("none (neutral pore)", NO_CHARGE)
        for c in LUMEN_CLOSURES:
            self.charge.addItem(f"{c}: {CLOSURE_LABELS[c]}", c)
        self.charge.setToolTip(
            "Round 7.11's placements of the lining side chains' charge on the "
            "lumen; dielectric (Round 7.13) puts every charged group in the "
            "box at its own centre, behind a protein of ε "
            "dielectric.eps_protein (two or three minutes). With a charge, "
            "the drop is the K+ electrochemical potential in linear response "
            "(conductivity σ e^−u).")
        self.charge.currentIndexChanged.connect(self._charge)
        row.addWidget(self.charge, 1)
        self.pairs = QCheckBox("salt bridges paired")
        self.pairs.setToolTip("Drop the lining groups that form a salt bridge "
                              "(Barlow & Thornton), both partners. Under "
                              "dielectric, unpaired is the bridge as a dipole.")
        self.pairs.toggled.connect(self._charge)
        row.addWidget(self.pairs)
        lay.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(QLabel("Colour by"))
        self.colour = QComboBox()
        for k, label in COLOURINGS.items():
            self.colour.addItem(label, k)
        self.colour.setToolTip(
            "Drop: 0 at the luminal bath (blue) to 1 at the cytosolic (red). "
            "Wall potential: a fixed ± display.lumen_potential_range kT/e, "
            "cation wells blue; zero everywhere in the neutral pore.")
        self.colour.currentIndexChanged.connect(lambda _: self.colour_changed.emit())
        row.addWidget(self.colour, 1)
        lay.addLayout(row)
        self.info = QLabel("")
        self.info.setWordWrap(True)
        lay.addWidget(self.info)
        self._enable()

    # -------------------------------------------------------------- state

    @property
    def closure(self) -> str | None:
        c = self.charge.currentData()
        return None if c == NO_CHARGE else c

    @property
    def colouring(self) -> str:
        return self.colour.currentData()

    def _toggled(self, on: bool) -> None:
        self._enable()
        self.toggled.emit(on)

    def _charge(self, *_) -> None:
        self._enable()
        self.charge_changed.emit()

    def _enable(self) -> None:
        on = self.show.isChecked()
        self.charge.setEnabled(on)
        self.colour.setEnabled(on)
        self.pairs.setEnabled(on and self.closure is not None)

    def view_state(self) -> dict:
        return {"charge": self.charge.currentData(), "pairs": self.pairs.isChecked(),
                "colour": self.colour.currentData()}

    def restore(self, d: dict) -> list[str]:
        """Set the controls without emitting; the caller re-requests."""
        notes: list[str] = []
        widgets = (self.charge, self.pairs, self.colour)
        for w in widgets:
            w.blockSignals(True)
        try:
            if "charge" in d:
                set_combo(self.charge, d["charge"], "lumen wall charge", notes)
            if "pairs" in d:
                set_check(self.pairs, d["pairs"], "lumen salt bridges", notes)
            if "colour" in d:
                set_combo(self.colour, d["colour"], "lumen colouring", notes)
        finally:
            for w in widgets:
                w.blockSignals(False)
        self._enable()
        return notes


def _mark(axes, s, bottom) -> None:
    for c in s.constrictions.values():
        for ax in axes:
            ax.axvline(c.z, color="#8a8f99", lw=0.6, ls=":")
        bottom.annotate(c.name, (c.z, 0.05), xytext=(3, 0),
                        textcoords="offset points", color="#d7dbe3", fontsize=7,
                        xycoords=("data", "axes fraction"))


def draw_lumen(canvas, f, s, charged=None) -> str:
    """Plot the lumen field ``f`` (and ``charged``, a
    :class:`~ip3r.physics.lumen_charge.ChargedLumen`) on ``canvas``; returns
    the panel's text."""
    from ..parameters import PARAMETERS as _P
    rows = 3 if charged is not None else 2
    axes = canvas.reset(rows, 1)
    area, drop = axes[0, 0], axes[1, 0]
    area.plot(f.z, f.area_3d, color=PALETTE[0], lw=1.4, label="3-D: lumen region")
    area.plot(f.z, f.area_1d, color=PALETTE[2], lw=1.0,
              label="1-D: π (r_free − r_ion)²")
    area.set_ylabel(f"{f.species} area (Å²)")
    area.set_title(f"{f.name}: where the voltage falls, S0's window", fontsize=8)
    drop.plot(f.z, f.drop_3d, color=PALETTE[0], lw=1.4, label="3-D (Laplace)")
    drop.plot(f.z, f.drop_1d, color=PALETTE[2], lw=1.0, label="1-D (∫dz/A)")
    if charged is not None:
        drop.plot(f.z, charged.drop_3d, color=PALETTE[1], lw=1.4,
                  label=f"3-D charged ({charged.closure}), K+")
        drop.plot(f.z, charged.drop_1d, color=PALETTE[3], lw=1.0, ls="--",
                  label="1-D charged (∫dz/(A·exp(−u))), K+")
    drop.set_ylabel("share of the\nwindow's drop")
    drop.set_ylim(-0.02, 1.02)
    bottom = drop
    if charged is not None:
        wall = axes[2, 0]
        wall.plot(f.z, charged.u_3d, color=PALETTE[1], lw=1.4,
                  label=f"3-D ({charged.closure}): plane mean")
        wall.plot(f.z, charged.u_1d, color=PALETTE[3], lw=1.0, ls="--",
                  label="1-D: local Donnan")
        wall.axhline(0.0, color="#8a8f99", lw=0.6)
        wall.set_ylabel("wall potential\n(kT/e)")
        bottom = wall
    bottom.set_xlabel("z along the four-fold axis (Å; luminal ← → cytosolic)")
    w = _P.value("lumen.constriction_half_width")
    lines = [f.summary() + "."]
    if charged is not None:
        lines.append(charged.summary() + ".")
    _mark([a for a in axes[:, 0]], s, drop)
    for c in s.constrictions.values():
        if not f.conducts:
            continue
        d3, d1 = f.drop_across(c.z, w)
        text = (f"{c.name} (z {c.z:+.1f} ± {w:.0f} Å): {d3:.0%} of the drop in "
                f"3-D, " + (f"{d1:.0%} in 1-D" if np.isfinite(d1) else "1-D shut"))
        if charged is not None:
            text += f"; {charged.drop_across(c.z, w):.0%} of the charged K+ drop"
        lines.append(text)
    for ax in axes[:, 0]:
        canvas.legend(ax, loc="upper left")
    canvas.draw_now()
    note = ("Neutral pore, the family's bath; each curve normalised to its own "
            "drop across the window." if charged is None else
            "The family's bath; drops normalised to their own window drop. "
            "Charged drops are K+'s electrochemical potential (linear "
            "response), not the electrical potential.")
    return "<br>".join(lines) + f"<br><i>{note}</i>"
