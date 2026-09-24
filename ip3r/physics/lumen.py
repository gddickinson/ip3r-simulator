"""SR Ca2+ content under the couplon: depletion as a terminator of release.

Round 6.6 found that the C channels, fitted to Murayama's bell, give no
release peak under any Mg2+ arrangement, and asked for a terminator that
switches on during the pulse. Emptying of the SR is the candidate this
module tests, the way Stern et al. 1997 (Fig. 20) did: the SR is one well-
mixed pool (Shirokova & Rios 1996, as Stern argues), every unitary current
is proportional to the content left, and the content falls with the
ensemble's release and refills with uptake:

    dS/dt = -k I(t) + (S0 - S) / tau,     k = density / (2 F)

with ``I`` the mean release current per couplon, ``density`` the couplons
per unit volume, and ``S`` the content as a myoplasmic concentration. The
pool is shared, so ``S`` follows the ensemble mean, not one couplon; the
path ``S(t)`` and the ensemble that releases along it are made consistent
by damped fixed-point iteration with the same seeds every time (common
random numbers). The couplon itself is exact given the path.

What this leaves out, on purpose: Stern's global cytosolic Ca2+ (a µM
boundary against tens of µM in the cleft), any luminal Ca2+ sensor on the
channel (no source for RyR1), and gradients inside the SR (Stern: small).
Uptake is first order in the deficit, with the time constant of Stern's
recovery after the pulse; their pump runs faster during it, so this
depletes somewhat more than theirs, which is the generous direction for
depletion as a terminator.

The measured rulers: a 100-ms pulse to +20 mV empties 50-60 % of the SR
(Rios 1993, Fig. 2), and a spark lowers the local free SR Ca2+ by at most
7.4 % (Launikonis 2006). Release corrected for depletion (divided by the
fraction left, Schneider et al. 1987) still has a peak in the fibre, so a
peak made only by emptying the store would vanish on correction:
:attr:`Depleted.corrected_peak_to_plateau` is the number that decides.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from ..parameters import PARAMETERS as _P
from .couplon import CouplonParams
from .ec_release import Ensemble, ensemble, with_gating

__all__ = ["FARADAY", "LumenParams", "Depleted", "mm_per_pa_ms",
           "content_path", "depleted_ensemble", "fig20_couplon",
           "depletion_panel", "pool_scan"]

FARADAY = 96485.33212     # C/mol


def _v(key: str):
    return field(default_factory=lambda: _P.value(key))


@dataclass
class LumenParams:
    content: float = _v("lumen.content")                # mM
    density: float = _v("lumen.couplon_density")        # 1/um^3
    refill_tau: float = _v("lumen.refill_tau")          # s
    iterations: float = _v("lumen.iterations")
    tolerance: float = _v("lumen.tolerance")
    damping: float = _v("lumen.damping")


def mm_per_pa_ms(lp: LumenParams) -> float:
    """SR content (mM) one couplon's 1 pA for 1 ms removes at this density:
    density [um^-3] = density * 1e15 per litre; 1 pA ms = 1e-15 C."""
    return lp.density * 1e15 * 1e-15 / (2.0 * FARADAY) * 1e3


def content_path(t: np.ndarray, flux: np.ndarray, lp: LumenParams
                 ) -> np.ndarray:
    """Content (mM) at each bin start of ``t`` (s), given the actual mean
    current ``flux`` (pA per couplon) over each bin: release removes
    ``k I dt``, then uptake relaxes the deficit exactly over the bin."""
    k = mm_per_pa_ms(lp)
    dt = np.diff(t, append=t[-1] + (t[1] - t[0] if len(t) > 1 else 0.0))
    decay = np.exp(-dt / lp.refill_tau)
    s = np.empty(len(t))
    s[0] = lp.content
    for i in range(len(t) - 1):
        after = s[i] - k * flux[i] * 1e3 * dt[i]
        s[i + 1] = lp.content - (lp.content - after) * decay[i]
    return np.maximum(s, 0.0)


@dataclass(frozen=True)
class Depleted:
    label: str
    e: Ensemble                 # the ensemble along the consistent path
    full: Ensemble              # the same couplons with the SR never emptying
    content: np.ndarray         # mM at each bin start
    s0: float
    iterations: int
    change: float               # the last iteration's largest move

    @property
    def fraction(self) -> np.ndarray:
        return self.content / self.s0

    def _window(self):
        pulse, win = _P.value("ec.pulse"), _P.value("ec.plateau_window")
        during = self.e.t < pulse
        return during, during & (self.e.t >= pulse - win)

    @property
    def released(self) -> float:
        """Fraction of the content gone at the end of the pulse."""
        during, _ = self._window()
        return 1.0 - float(self.fraction[during][-1])

    @staticmethod
    def _ratio(y, during, plateau) -> float:
        p = float(np.nanmean(y[plateau]))
        return float(np.nanmax(y[during])) / p if p > 0 else float("inf")

    @property
    def peak_to_plateau(self) -> float:
        return self._ratio(self.e.flux, *self._window())

    @property
    def corrected_peak_to_plateau(self) -> float:
        """The same ratio after dividing by the fraction left (Schneider)."""
        return self._ratio(self.corrected(), *self._window())

    @property
    def po_peak_to_plateau(self) -> float:
        """C open probability, peak over plateau: the permeability's own
        ratio, with no division by the content."""
        return self._ratio(self.e.c_po, *self._window())

    @property
    def full_peak_to_plateau(self) -> float:
        return self._ratio(self.full.flux, *self._window())

    def corrected(self) -> np.ndarray:
        f = np.where(self.fraction > 0, self.fraction, np.nan)
        return self.e.flux / f

    @property
    def c_after(self) -> float:
        after = self.e.t >= _P.value("ec.pulse") + 0.5 * _P.value("ec.after")
        return float(self.e.c_po[after].mean()) if after.any() else float("nan")

    def at_ms(self, ms: float) -> float:
        """Content (mM) at ``ms`` after the step."""
        return float(np.interp(ms * 1e-3, self.e.t, self.content))

    def row(self) -> str:
        return (f"{self.label:34s} {self.e.v:5.0f} mV  content at 100 ms "
                f"{self.at_ms(100):.2f} mM (released {self.released:.2f})  "
                f"flux peak/plateau {self.peak_to_plateau:5.2f} (full SR "
                f"{self.full_peak_to_plateau:5.2f}, corrected "
                f"{self.corrected_peak_to_plateau:5.2f}; C Po "
                f"{self.po_peak_to_plateau:5.2f})  C after "
                f"{self.c_after:.3f}  [{self.iterations} it]")


def depleted_ensemble(v: float, pp: CouplonParams | None = None,
                      lp: LumenParams | None = None, trials: int | None = None,
                      label: str = "", seed0: int = 0) -> Depleted:
    """The ensemble at ``v`` mV with the SR emptying along a path that is
    consistent with the ensemble's own release."""
    pp = pp or CouplonParams()
    lp = lp or LumenParams()
    full = ensemble(v, pp, trials, seed0)
    t = full.t
    frac = content_path(t, full.flux, lp) / lp.content
    e, change, it = full, np.inf, 0
    while it < int(round(lp.iterations)):
        it += 1
        e = ensemble(v, pp, trials, seed0, scale=(t, frac))
        new = content_path(t, e.flux, lp) / lp.content
        change = float(np.max(np.abs(new - frac)))
        frac = frac + lp.damping * (new - frac)
        if change < lp.tolerance:
            break
    return Depleted(label, e, full, content_path(t, e.flux, lp), lp.content,
                    it, change)


def fig20_couplon(pp: CouplonParams | None = None) -> CouplonParams:
    """Stern's Fig. 20 couplon: ``lumen.fig20_channels`` C channels."""
    pp = pp or CouplonParams()
    n = _P.value("lumen.fig20_channels")
    return replace(pp, cleft=replace(pp.cleft, n_channels=n))


def depletion_panel(configs: dict, voltages, trials: int | None = None,
                    small: bool = False) -> list[Depleted]:
    """Every configuration (label -> SternParams) at every voltage; ``small``
    uses the Fig. 20 couplon instead of the 60-channel one."""
    rows = []
    for label, sp in configs.items():
        pp = with_gating(sp)
        pp = fig20_couplon(pp) if small else pp
        for v in voltages:
            rows.append(depleted_ensemble(v, pp, trials=trials, label=label))
    return rows



def pool_scan(configs: dict, v: float, factors, trials: int | None = None,
              small: bool = True) -> list[tuple[float, Depleted]]:
    """Each configuration with the pool ``factor`` times ``lumen.content``:
    a larger pool empties by a smaller fraction for the same release. The
    factor at which Stern's constants release Rios's measured 50-60 % is
    the pool the measurement allows."""
    out = []
    for label, sp in configs.items():
        pp = with_gating(sp)
        pp = fig20_couplon(pp) if small else pp
        for f in factors:
            lp = replace(LumenParams(), content=f * _P.value("lumen.content"))
            out.append((f, depleted_ensemble(v, pp, lp, trials, label=label)))
    return out
