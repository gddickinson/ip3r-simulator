"""The Ca2+ around a puff site: Cao et al. 2014's pools, with fluo-4.

Three pools in cytosolic units, as in Cao et al. 2014's model code, and the
Ca2+-bound indicator of Cao et al. 2013's Eq. 12 placed in the microdomain,
where the puff is imaged:

    cs  = gamma2 (ct - (cb + b)/gamma1 - c)          store (ER units)
    dc  = Jdiff + Jleak - Jserca + Jin - Jpm          cytosol
    dcb = gamma1 (Jipr - Jdiff) - Jdye                microdomain
    dct = Jin - Jpm                                   total (cytosolic units)
    db  = Jdye = kon (B - b) cb - koff b              fluo-4 bound

    Jipr  = k_ipr N_open (cs - cb)     Jdiff = k_diff (cb - c)
    Jleak = k_leak (cs - c)            Jserca = V c^n / (c^n + K^n)
    Jin   = j_leak_in + v_rocc p + V_socc K^4 / (cs^4 + K^4)
    Jpm   = V_pm c^2 / (K_pm^2 + c^2)

The bound indicator is counted in the total, so the dye takes its Ca2+ from
the microdomain and not from the store. An open receptor sees its own mouth,
``mouth_per_store x cs`` (the code's ``120 (cs/100)``); a closed one sees
``cb``. With ``clamp_store`` the store is held at rest, Cao 2013's
assumption of a constant single-channel flux, and ``ct`` follows whatever the
other pools do. The ``"bath"`` clamp also holds the cytosol, which leaves the
microdomain exchanging with fixed surroundings (``CLAMPS``).

The two limits are what the tests hold the module to: at rest the right-hand
side is zero, and with the store clamped and the dye at equilibrium each
open receptor raises the microdomain by ``k_ipr (cs - cb) / k_diff`` (about
0.11 µM), the mean-field model's coupling.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.optimize import brentq

from ..parameters import PARAMETERS as _P

__all__ = ["DomainParams", "Rest", "start_store", "rest_state",
           "whole_cell_rest", "rhs", "rk4_step",
           "CLAMPS", "clamped", "store", "fluorescence", "coupling_per_open"]


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class DomainParams:
    gamma1: float = _v("domain.gamma1")
    gamma2: float = _v("domain.gamma2")
    k_ipr: float = _v("domain.k_ipr")
    k_diff: float = _v("domain.k_diff")
    k_leak: float = _v("domain.k_leak")
    v_serca: float = _v("domain.v_serca")
    k_serca: float = _v("domain.k_serca")
    n_serca: float = _v("domain.n_serca")
    v_pm: float = _v("domain.v_pm")
    k_pm: float = _v("domain.k_pm")
    j_leak_in: float = _v("domain.j_leak_in")
    v_rocc: float = _v("domain.v_rocc")
    v_socc: float = _v("domain.v_socc")
    k_socc: float = _v("domain.k_socc")
    mouth_per_store: float = _v("domain.mouth_per_store")
    ct0: float = _v("domain.ct0")
    c0: float = _v("puff.ca_rest")     # the code's starting c = cb
    fluo_total: float = _v("domain.fluo_total")
    fluo_kon: float = _v("domain.fluo_kon")
    fluo_koff: float = _v("domain.fluo_koff")
    clamp: str = "none"                 # "none" | "store" | "bath" (CLAMPS)
    cs_clamp: float = float("nan")      # set by ``clamped``

    def scaled(self, k_ipr_scale: float) -> "DomainParams":
        return replace(self, k_ipr=self.k_ipr * k_ipr_scale)


@dataclass(frozen=True)
class Rest:
    c: float      # cytosol = microdomain, µM
    cs: float     # store, µM (ER units)
    b: float      # bound fluo-4, µM
    ct: float     # total, cytosolic units
    y: np.ndarray  # (c, cb, ct, b)


def _serca(c, d):
    return d.v_serca * c ** d.n_serca / (c ** d.n_serca + d.k_serca ** d.n_serca)


def _pm(c, d):
    return d.v_pm * c * c / (d.k_pm * d.k_pm + c * c)


def _j_in(cs, p, d):
    k4 = d.k_socc ** 4
    return d.j_leak_in + d.v_rocc * p + d.v_socc * k4 / (max(cs, 0.0) ** 4 + k4)


def _bound(c, d):
    return d.fluo_total * c / (c + d.fluo_koff / d.fluo_kon)


def start_store(d: DomainParams) -> float:
    """The store the code starts from: ``gamma2 (ct0 - c0/gamma1 - c0)``
    with ``c = cb = c0`` (about 449 µM)."""
    return d.gamma2 * (d.ct0 - d.c0 / d.gamma1 - d.c0)


def _rest(c, cs, d):
    b = _bound(c, d)
    ct = c + (c + b) / d.gamma1 + cs / d.gamma2
    return Rest(c, cs, b, ct, np.array([c, c, ct, b]))


def rest_state(p: float, d: DomainParams) -> Rest:
    """Every receptor shut, the store at ``start_store``, at IP3 ``p`` (µM).

    The cytosol (= microdomain) is at its steady state against that store:
    ``k_leak (cs - c) + Jin = Jserca + Jpm``, one equation in ``c``. With the
    store clamped this is the steady state. With the store free it is not
    quite one: the total drifts at ``Jin - Jpm``, a few tens of nM/s, the
    code's own slow approach to ``whole_cell_rest``.
    """
    cs = start_store(d)
    c = brentq(lambda c: d.k_leak * (cs - c) + _j_in(cs, p, d) - _pm(c, d)
               - _serca(c, d), 1e-6, 50.0, xtol=1e-14)
    return _rest(c, cs, d)


def whole_cell_rest(p: float, d: DomainParams) -> Rest:
    """The true steady state with every receptor shut: the ER leak balances
    SERCA (``cs = c + Jserca/k_leak``) and entry balances extrusion."""
    def cs_of(c):
        return c + _serca(c, d) / d.k_leak

    c = brentq(lambda c: _j_in(cs_of(c), p, d) - _pm(c, d), 1e-6, 50.0,
               xtol=1e-14)
    return _rest(c, cs_of(c), d)


#: What a run holds fixed: nothing; the store (Cao 2013's constant flux);
#: or the store and the cytosol, so the microdomain exchanges with a fixed
#: bath (the mean-field cluster's assumption, with the microdomain's
#: kinetics and the dye kept).
CLAMPS = ("none", "store", "bath")


def clamped(d: DomainParams, rest: Rest, clamp: str = "store") -> DomainParams:
    """``d`` with ``clamp`` applied at ``rest`` (see ``CLAMPS``)."""
    if clamp not in CLAMPS:
        raise ValueError(f"clamp must be one of {CLAMPS}, not {clamp!r}")
    return replace(d, clamp=clamp, cs_clamp=rest.cs)


def store(y, d: DomainParams) -> float:
    """Store Ca2+ (ER units) of state ``y = (c, cb, ct, b)``."""
    if d.clamp != "none":
        return d.cs_clamp
    return d.gamma2 * (y[2] - (y[1] + y[3]) / d.gamma1 - y[0])


def rhs(y, n_open: int, p: float, d: DomainParams):
    """``dy/dt`` for ``y = (c, cb, ct, b)`` with ``n_open`` receptors open."""
    c, cb, _, b = y
    cs = store(y, d)
    j_diff = d.k_diff * (cb - c)
    j_in = _j_in(cs, p, d)
    j_pm = _pm(c, d)
    j_ipr = d.k_ipr * n_open * (cs - cb)
    j_dye = d.fluo_kon * (d.fluo_total - b) * cb - d.fluo_koff * b
    dc = j_diff + d.k_leak * (cs - c) - _serca(c, d) + j_in - j_pm
    if d.clamp == "bath":
        dc = 0.0
    dcb = d.gamma1 * (j_ipr - j_diff) - j_dye
    if d.clamp != "none":  # ct is bookkeeping: whatever keeps cs fixed
        dct = dc + (dcb + j_dye) / d.gamma1
    else:
        dct = j_in - j_pm
    return (dc, dcb, dct, j_dye)


def rk4_step(y, n_open: int, p: float, d: DomainParams, dt: float):
    """One fourth-order Runge-Kutta step (Cao et al.'s integrator)."""
    def add(a, k, h):
        return (a[0] + h * k[0], a[1] + h * k[1], a[2] + h * k[2], a[3] + h * k[3])

    k1 = rhs(y, n_open, p, d)
    k2 = rhs(add(y, k1, dt / 2), n_open, p, d)
    k3 = rhs(add(y, k2, dt / 2), n_open, p, d)
    k4 = rhs(add(y, k3, dt), n_open, p, d)
    return tuple(y[i] + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                 for i in range(4))


def fluorescence(b, rest: Rest):
    """F/F0: bound indicator over its resting value (Cao 2013's ratio)."""
    return np.asarray(b, float) / rest.b if rest.b > 0 else np.full_like(
        np.asarray(b, float), math.nan)


def coupling_per_open(rest: Rest, d: DomainParams) -> float:
    """Steady microdomain rise per open receptor with the store clamped,
    ``k_ipr (cs - cb) / k_diff`` to first order: the mean-field coupling."""
    return d.k_ipr * (rest.cs - rest.c) / d.k_diff
