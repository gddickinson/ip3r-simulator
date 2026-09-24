"""Per-atom colours for every colouring the viewer offers.

Scales that carry a quantity are **fixed**, never auto-ranged: a
conservation map auto-ranged per structure would paint the least conserved
residue of one deposit the same colour as the least conserved of another,
and the two could not be compared. Where a residue has no value (unresolved
in the alignment, outside the reference numbering) it is drawn **grey, never
the low end of the scale** — the luminal loop is both the least conserved
element and the worst resolved, and those two must not look alike (the rule
``ip3r_genes`` S24 applied to its painted structures).
"""

from __future__ import annotations

import numpy as np

from ..core.annotations import ELEMENT_COLORS, element_array
from ..core.structure import Structure

__all__ = ["CHAIN_PALETTE", "MISSING", "chain_colors", "element_colors",
           "ramp", "constraint_colors", "bfactor_colors", "value_colors",
           "uniform_color", "CONSERVATION_RANGE", "displacement_colors",
           "SHELL_COLORS", "shell_colors"]

#: Four subunits, four distinguishable hues.
CHAIN_PALETTE = np.array([
    (0.36, 0.62, 0.95), (0.95, 0.55, 0.30), (0.45, 0.80, 0.50),
    (0.80, 0.50, 0.90), (0.95, 0.80, 0.30), (0.40, 0.85, 0.85),
    (0.90, 0.45, 0.55), (0.65, 0.65, 0.70)], dtype=np.float32)

#: Colour for residues with no value on a quantitative scale.
MISSING = np.array((0.42, 0.43, 0.46), np.float32)

#: The fixed JSD range the conservation ramp spans. The first choice, 0.30,
#: painted almost the whole receptor red (deep-layer residues sit between
#: ~0.6 and ~0.9); 0.50 sits just above the luminal-loop mean (0.45-0.49),
#: so that loop reads as the variable end and the ramp resolves everything
#: else. 0.95 is above every residue the deep layer calls invariant.
CONSERVATION_RANGE = (0.5, 0.95)

# A perceptually ordered blue -> pale -> red ramp (low -> high).
_RAMP = np.array([(0.23, 0.30, 0.75), (0.55, 0.69, 0.99), (0.87, 0.87, 0.87),
                  (0.96, 0.60, 0.48), (0.71, 0.02, 0.15)], np.float32)


def ramp(t: np.ndarray) -> np.ndarray:
    """Map ``t`` in [0, 1] onto the ramp; NaN becomes :data:`MISSING`."""
    t = np.asarray(t, np.float64)
    out = np.tile(MISSING, (len(t), 1)).astype(np.float32)
    ok = np.isfinite(t)
    x = np.clip(t[ok], 0.0, 1.0) * (len(_RAMP) - 1)
    lo = np.floor(x).astype(int)
    hi = np.minimum(lo + 1, len(_RAMP) - 1)
    f = (x - lo)[:, None]
    out[ok] = (_RAMP[lo] * (1 - f) + _RAMP[hi] * f).astype(np.float32)
    return out


def chain_colors(st: Structure) -> np.ndarray:
    idx = {c: i for i, c in enumerate(st.chains)}
    k = np.array([idx[str(c)] for c in st.chain]) % len(CHAIN_PALETTE)
    return CHAIN_PALETTE[k]


def element_colors(st: Structure, paralog: str | None) -> np.ndarray:
    """Colour by functional element; grey where numbering does not apply."""
    out = np.tile(MISSING, (st.n_atoms, 1))
    if paralog is None:
        return out
    names = element_array(paralog, st.res_seq.astype(int))
    for name, rgb in ELEMENT_COLORS.items():
        out[(names == name) & ~st.hetero] = rgb
    return out


def constraint_colors(values: np.ndarray, lo: float | None = None,
                      hi: float | None = None) -> np.ndarray:
    lo = CONSERVATION_RANGE[0] if lo is None else lo
    hi = CONSERVATION_RANGE[1] if hi is None else hi
    return ramp((np.asarray(values, float) - lo) / (hi - lo))


def bfactor_colors(st: Structure) -> np.ndarray:
    """B-factor (or pLDDT for a prediction), ranged over this structure."""
    b = st.b_factor.astype(float)
    lo, hi = np.percentile(b, [2, 98]) if len(b) else (0.0, 1.0)
    return ramp((b - lo) / max(hi - lo, 1e-6))


def value_colors(values: np.ndarray) -> np.ndarray:
    """An arbitrary per-atom scalar, auto-ranged — for mode amplitudes."""
    v = np.asarray(values, float)
    ok = np.isfinite(v)
    if not ok.any():
        return np.tile(MISSING, (len(v), 1))
    lo, hi = np.percentile(v[ok], [2, 98])
    return ramp((v - lo) / max(hi - lo, 1e-9))


def displacement_colors(values: np.ndarray) -> np.ndarray:
    """Displacement between two states on a fixed 0..``display.displacement_max``
    scale; NaN (not measured) is grey."""
    from ..parameters import PARAMETERS as _P
    return ramp(np.asarray(values, float) / _P.value("display.displacement_max"))


#: S22's four ligand shells, nearest hottest; the fixed ramp sampled at four
#: points so the colours are the conservation ramp's own.
SHELL_COLORS = ramp(np.array([1.0, 0.8, 0.62, 0.3]))


def shell_colors(distance: np.ndarray) -> np.ndarray:
    """Per atom: the colour of its ligand shell (edges from the registry);
    NaN or beyond the search radius is grey."""
    from ..parameters import PARAMETERS as _P
    d = np.asarray(distance, float)
    out = np.tile(MISSING, (len(d), 1)).astype(np.float32)
    lo = 0.0
    for rgb, key in zip(SHELL_COLORS, ("ligand.contact_cutoff", "ligand.shell_second_edge",
                                       "ligand.shell_third_edge", "ligand.shell_radius")):
        hi = _P.value(key)
        out[(d >= lo) & (d < hi)] = rgb
        lo = hi
    return out


def uniform_color(st: Structure, rgb=(0.55, 0.62, 0.75)) -> np.ndarray:
    return np.tile(np.asarray(rgb, np.float32), (st.n_atoms, 1))
