"""Charge-space competition: the excess chemical potential of a crowded,
charged fluid (Round 7.17).

Every closure so far treats the ions as points in a mean field. Such a
Donnan partition fills a filter of 20 M fixed charge with 20 M K+, and
nothing penalises it: a point ion takes no room, and the field it feels is
the average one. Nonner, Catacuzzeno & Eisenberg 2000 and Gillespie 2008
explain the Ca2+ selectivity of the L-type and ryanodine-receptor filters
with the two terms that picture omits. Both are local functions of the
densities of every species in a slice:

* **excluded volume** (hard spheres): inserting an ion of diameter sigma
  into a dense fluid costs work. Here the Boublik-Mansoori-Carnahan-
  Starling-Leland (BMCSL) mixture, whose one-component limit is Carnahan-
  Starling, differentiated from its free energy (:func:`hard_sphere`);
* **screening** (the mean spherical approximation, Blum 1975): an ion
  lowers its energy by the arrangement of its neighbours, and a small,
  divalent ion does it best (:func:`msa`).

Their sum, per species and per slice, is :func:`excess` (in kT). The
filter's own carboxylates are fluid species too, as in Gillespie's model:
two half-charged oxygens (``csc.oxygen_diameter``) per acid, and one sphere
per base (``csc.base_diameter``), at the density the wall charge map gives
them. They screen (MSA) and, when ``csc.structural_volume`` is 1, take up
room (hard spheres). The deposit's profile already excludes their atoms,
so 1 counts their volume twice; 0 omits it. Both are reported.

Water is an uncharged hard sphere at its bath density
(``csc.water_concentration``, ``csc.water_diameter``). It is displaced
where ions crowd in, and that is where most of the room comes from.

:func:`partition` solves one slice's equilibrium with its reservoir under
local electroneutrality: c_i = a_i exp(-z_i psi - mu_i(c)), sum z_i c_i +
X = 0, with a_i the activity (bath concentration times exp(mu_bath)). With
both terms switched off it is exactly the Donnan closure.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["hard_sphere", "msa", "excess", "bjerrum", "Structural",
           "structural", "Fluid", "bath_activity", "partition", "Partition",
           "TO_A3", "hs_free_energy", "msa_free_energy"]

#: mol/m^3 -> particles per cubic angstrom.
TO_A3 = 6.02214076e23 * 1e-30
_EPS0 = 8.8541878128e-12
_E = 1.602176634e-19
_KB = 1.380649e-23


def bjerrum(eps: float | None = None, temperature: float | None = None
            ) -> float:
    """Bjerrum length e^2 / (4 pi eps0 eps kT), A."""
    eps = _P.value("csc.permittivity") if eps is None else eps
    t = _P.value("permeation.temperature") if temperature is None else temperature
    return _E ** 2 / (4 * np.pi * _EPS0 * eps * _KB * t) * 1e10


# ------------------------------------------------------------ hard spheres
def _xi(rho, sigma):
    """xi_n = (pi/6) sum_j rho_j sigma_j^n, n = 0..3 (rho (n_sp, m) A^-3)."""
    return [np.pi / 6.0 * (rho * sigma[:, None] ** n).sum(axis=0)
            for n in range(4)]


def hs_free_energy(rho, sigma) -> np.ndarray:
    """BMCSL excess free energy density, kT per A^3 (for the tests)."""
    x0, x1, x2, x3 = _xi(np.asarray(rho, float), np.asarray(sigma, float))
    d, lg = 1.0 - x3, np.log1p(-x3)
    return 6.0 / np.pi * ((x2 ** 3 / x3 ** 2 - x0) * lg + 3 * x1 * x2 / d
                          + x2 ** 3 / (x3 * d ** 2))


def hard_sphere(rho, sigma) -> np.ndarray:
    """BMCSL excess chemical potential of each species, kT, (n_sp, m).

    ``rho`` in A^-3, ``sigma`` (n_sp,) in A. The derivative of
    :func:`hs_free_energy` taken analytically through the xi_n.
    """
    rho = np.asarray(rho, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    x0, x1, x2, x3 = _xi(rho, sigma)
    out = np.zeros_like(rho)
    ok = x3 > 1e-12
    if not np.any(ok):
        return out
    x0, x1, x2, x3 = (v[ok] for v in (x0, x1, x2, x3))
    d, lg = 1.0 - x3, np.log1p(-x3)
    f0 = -lg
    f1 = 3 * x2 / d
    f2 = 3 * x2 ** 2 * lg / x3 ** 2 + 3 * x1 / d + 3 * x2 ** 2 / (x3 * d ** 2)
    f3 = (-2 * x2 ** 3 * lg / x3 ** 3 - (x2 ** 3 / x3 ** 2 - x0) / d
          + 3 * x1 * x2 / d ** 2
          + x2 ** 3 * (-1.0 / (x3 ** 2 * d ** 2) + 2.0 / (x3 * d ** 3)))
    s = sigma[:, None]
    out[:, ok] = f0 + f1 * s + f2 * s ** 2 + f3 * s ** 3
    return out


# ------------------------------------------------------ screening (MSA)
def _msa_gamma(rho, z, sigma, lb, tol=1e-13, max_iter=500):
    """Blum's screening parameter Gamma (1/A) and eta, per slice."""
    rho, z = np.asarray(rho, float), np.asarray(z, float)
    s = sigma[:, None]
    kappa2 = 4 * np.pi * lb * (rho * z ** 2).sum(axis=0)
    gamma = 0.5 * np.sqrt(kappa2)
    delta = 1.0 - np.pi / 6.0 * (rho * s ** 3).sum(axis=0)
    eta = np.zeros_like(gamma)
    for _ in range(max_iter):
        g = 1.0 + gamma * s
        omega = 1.0 + np.pi / (2 * delta) * (rho * s ** 3 / g).sum(axis=0)
        eta = np.pi / (2 * delta * omega) * (rho * s * z / g).sum(axis=0)
        new = np.sqrt(np.pi * lb * (rho * ((z - eta * s ** 2) / g) ** 2)
                      .sum(axis=0))
        change = np.max(np.abs(new - gamma)) if new.size else 0.0
        gamma = 0.5 * gamma + 0.5 * new
        if change < tol:
            break
    g = 1.0 + gamma * s
    omega = 1.0 + np.pi / (2 * delta) * (rho * s ** 3 / g).sum(axis=0)
    eta = np.pi / (2 * delta * omega) * (rho * s * z / g).sum(axis=0)
    return gamma, eta


def msa(rho, z, sigma, lb: float) -> np.ndarray:
    """MSA electrostatic excess chemical potential, kT, (n_sp, m).

    ``z`` is (n_sp, m) or (n_sp,): a structural species' valence may vary
    along the pore (a partly protonated ring). Nonner et al. 2000's form of
    Blum's result: -l_B [Gamma z^2/(1+Gamma s) + eta s ((2z - eta s^2)
    /(1+Gamma s) + eta s^2/3)].
    """
    rho = np.asarray(rho, float)
    z = np.broadcast_to(np.asarray(z, float).reshape(len(sigma), -1),
                        rho.shape)
    sigma = np.asarray(sigma, float)
    gamma, eta = _msa_gamma(rho, z, sigma, lb)
    s = sigma[:, None]
    g = 1.0 + gamma * s
    return -lb * (gamma * z ** 2 / g
                  + eta * s * ((2 * z - eta * s ** 2) / g + eta * s ** 2 / 3.0))


def msa_free_energy(rho, z, sigma, lb: float) -> np.ndarray:
    """MSA excess free energy density, kT per A^3 (Blum & Hoye 1977): the
    energy plus Gamma^3 / 3pi. For the tests: its derivative is :func:`msa`."""
    rho = np.asarray(rho, float)
    z = np.broadcast_to(np.asarray(z, float).reshape(len(sigma), -1),
                        rho.shape)
    gamma, eta = _msa_gamma(rho, z, np.asarray(sigma, float), lb)
    s = np.asarray(sigma, float)[:, None]
    g = 1.0 + gamma * s
    energy = -lb * (gamma * (rho * z ** 2 / g).sum(axis=0)
                    + eta * (rho * z * s / g).sum(axis=0))
    return energy + gamma ** 3 / (3 * np.pi)


# ----------------------------------------------------------- the fluid
@dataclass(frozen=True)
class Structural:
    """The wall's groups as fluid species on the slices: density (mol/m^3,
    (k, n)), valence ((k, n): charge / count, so a mean protonation state
    enters as a fractional valence) and diameter (A, (k,))."""

    density: np.ndarray
    valence: np.ndarray
    diameter: np.ndarray
    names: tuple[str, ...] = ()

    @staticmethod
    def empty(n: int) -> "Structural":
        return Structural(np.zeros((0, n)), np.zeros((0, n)), np.zeros(0))

    @staticmethod
    def from_fixed(fixed) -> "Structural":
        """The wall from a net charge map alone: acid oxygens where it is
        negative (two per charge, -1/2 each), bases where positive."""
        f = np.asarray(fixed, float)
        neg, pos = np.maximum(-f, 0.0), np.maximum(f, 0.0)
        return Structural(np.array([2.0 * neg, pos]),
                          np.array([np.full_like(f, -0.5), np.ones_like(f)]),
                          np.array([_P.value("csc.oxygen_diameter"),
                                    _P.value("csc.base_diameter")]),
                          ("O", "N"))

    def charge(self) -> np.ndarray:
        """Signed charge density, mol/m^3: must equal the fixed charge."""
        return (self.density * self.valence).sum(axis=0)

    def take(self, order) -> "Structural":
        return Structural(self.density[:, order], self.valence[:, order],
                          self.diameter, self.names)


def structural(groups, z_A, radius_A, smoothing: float | None = None
               ) -> Structural:
    """Carboxylate oxygens and base groups of ``groups`` on the slices,
    Gaussian along the axis and divided by the same floored area as
    :func:`ip3r.physics.pore_charge.map_charge`, so their charge is the
    fixed-charge map exactly."""
    from .pore_charge import AVOGADRO, charge_per_length
    smoothing = (_P.value("pore_charge.smoothing") if smoothing is None
                 else smoothing)
    floor = _P.value("permeation.radius_potassium")
    area = np.pi * np.maximum(np.asarray(radius_A, float), floor) ** 2
    to_molar = 1e30 / AVOGADRO / area
    rows, names = [], []
    for name, kinds, per_group, key in (
            ("O", ("ASP", "GLU"), 2.0, "csc.oxygen_diameter"),
            ("N", ("LYS", "ARG"), 1.0, "csc.base_diameter")):
        mine = [g for g in groups if g.res_name in kinds]
        if not mine:
            continue
        count = per_group * charge_per_length(
            [_unit(g) for g in mine], z_A, smoothing) * to_molar
        charge = charge_per_length(mine, z_A, smoothing) * to_molar
        valence = np.divide(charge, count, out=np.zeros_like(count),
                            where=count > 0)
        rows.append((count, valence, _P.value(key)))
        names.append(name)
    if not rows:
        return Structural.empty(len(z_A))
    return Structural(np.array([r[0] for r in rows]),
                      np.array([r[1] for r in rows]),
                      np.array([r[2] for r in rows]), tuple(names))


def _unit(g):
    from dataclasses import replace
    return replace(g, charge=1.0)


@dataclass(frozen=True)
class Fluid:
    """What :func:`excess` needs besides the mobile densities: the mobile
    species' valences and diameters, the wall's groups, water, l_B, and
    which terms are on."""

    valence: np.ndarray                # (n,)
    diameter: np.ndarray               # (n,), A
    wall: Structural
    water: float                       # bath water, mol/m^3 (0 = none)
    water_diameter: float
    lb: float                          # A
    hs: bool = True
    msa: bool = True
    wall_volume: float = 1.0           # 0..1: structural_volume

    @staticmethod
    def of(species, wall: Structural, hs: bool = True, use_msa: bool = True,
           wall_volume: float | None = None, eps: float | None = None
           ) -> "Fluid":
        return Fluid(np.array([s.valence for s in species], float),
                     np.array([2.0 * s.radius for s in species], float),
                     wall, _P.value("csc.water_concentration") * 1000.0,
                     _P.value("csc.water_diameter"), bjerrum(eps), hs,
                     use_msa, (_P.value("csc.structural_volume")
                               if wall_volume is None else wall_volume))


def excess(c, water, fluid: Fluid, wall: Structural | None = None
           ) -> tuple[np.ndarray, np.ndarray]:
    """(mobile mu (n, m), water mu (m,)) in kT, from mobile concentrations
    ``c`` (n, m) and water (m,), mol/m^3."""
    wall = fluid.wall if wall is None else wall
    n, m = c.shape
    mu = np.zeros((n, m))
    mu_w = np.zeros(m)
    if fluid.hs:
        rho = np.vstack([c, water[None, :], wall.density * fluid.wall_volume]
                        ) * TO_A3
        sig = np.concatenate([fluid.diameter, [fluid.water_diameter],
                              wall.diameter])
        h = hard_sphere(rho, sig)
        mu += h[:n]
        mu_w += h[n]
    if fluid.msa and fluid.lb > 0:
        rho = np.vstack([c, wall.density]) * TO_A3
        z = np.vstack([np.repeat(fluid.valence[:, None], m, axis=1),
                       wall.valence])
        mu += msa(rho, z, np.concatenate([fluid.diameter, wall.diameter]),
                  fluid.lb)[:n]
    return mu, mu_w


def bath_activity(bath, fluid: Fluid) -> tuple[np.ndarray, float]:
    """Activities (mol/m^3) of the mobile species and of water in a bath
    of concentrations ``bath`` (n,) mol/m^3 with no wall."""
    bath = np.asarray(bath, float).reshape(-1, 1)
    empty = Structural.empty(1)
    mu, mu_w = excess(bath, np.array([fluid.water]), fluid, empty)
    return (bath[:, 0] * np.exp(mu[:, 0]),
            float(fluid.water * np.exp(mu_w[0])))


#: Largest packing fraction a Newton step may reach; a step that would
#: cross it is halved in that slice (hard spheres diverge at 1).
_PACKING_CAP = 0.9
#: Floor on a reservoir activity, mol/m^3 (keeps ln a finite where the
#: Nernst-Planck solve leaves a species at zero).
_A_FLOOR = 1e-20


def packing(c, water, fluid: Fluid, wall: Structural) -> np.ndarray:
    """Packing fraction xi_3 of each slice (the volume the hard spheres fill)."""
    rho = np.vstack([c, water[None, :], wall.density * fluid.wall_volume]
                    ) * TO_A3
    sig = np.concatenate([fluid.diameter, [fluid.water_diameter],
                          wall.diameter])
    return np.pi / 6.0 * (rho * sig[:, None] ** 3).sum(axis=0)


@dataclass
class Partition:
    """One equilibrium over the slices."""

    c: np.ndarray                      # (n, m) mol/m^3
    water: np.ndarray                  # (m,) mol/m^3
    psi: np.ndarray                    # (m,) V
    mu: np.ndarray                     # (n, m) kT
    converged: bool
    iterations: int

    def offset(self, valence, thermal) -> np.ndarray:
        """Each species' potential offset, V: psi + mu kT / (z e)."""
        return self.psi[None, :] + self.mu * thermal / np.asarray(
            valence, float)[:, None]


def partition(activity, fixed, fluid: Fluid, thermal: float,
              water_activity: float, wall: Structural | None = None,
              initial: Partition | tuple | None = None,
              tol: float | None = None,
              max_iter: int | None = None) -> Partition:
    """Local equilibrium of every slice with its reservoir (activities
    ``activity`` (n, m), mol/m^3) under electroneutrality with ``fixed``.
    ``initial``: a previous :class:`Partition` (warm start) or a guess
    ``(c, water)``; with none the first guess is ``a exp(-mu)`` for
    the ions at infinite dilution in pure water, close to any bath.

    Newton on u = (ln c_i, ln c_water, psi/phi_T) in every slice at once:
    ln c_i + z_i psi + mu_i(c) = ln a_i for the ions and water, and
    sum z_i c_i + X = 0. The Jacobian of mu is taken by finite
    differences in ln c (n + 1 evaluations of :func:`excess`); a step is
    capped at ``csc.max_step`` and halved in any slice where it would
    over-pack or not lower the residual. With both terms off, mu = 0 and
    the root is the Donnan partition.
    """
    wall = fluid.wall if wall is None else wall
    tol = _P.value("csc.tolerance") if tol is None else tol
    max_iter = int(_P.value("csc.max_iterations") if max_iter is None
                   else max_iter)
    cap = _P.value("csc.max_step")
    a = np.maximum(np.asarray(activity, float), _A_FLOOR)
    fixed = np.asarray(fixed, float)
    z = fluid.valence
    n, m = a.shape
    if isinstance(initial, Partition):
        c, w = initial.c.copy(), initial.water.copy()
        psi = initial.psi / thermal
    else:
        c, w = ((np.asarray(initial[0], float).copy(),
                 np.asarray(initial[1], float).copy()) if initial is not None
                else (a * np.exp(-excess(np.zeros_like(a),
                                         np.full(m, fluid.water), fluid,
                                         Structural.empty(m))[0]),
                      np.full(m, float(fluid.water))))
        psi = np.zeros(m)
    has_water = fluid.water > 0 and fluid.hs
    u = np.vstack([np.log(np.maximum(c, _A_FLOOR)),
                   np.log(np.maximum(w, _A_FLOOR))[None, :], psi[None, :]])
    ln_a = np.vstack([np.log(a), np.full((1, m), np.log(max(water_activity,
                                                            _A_FLOOR)))])

    def unpack(v):
        return np.exp(v[:n]), np.exp(v[n]), v[n + 1]

    def residual(v):
        cc, ww, pp = unpack(v)
        mu, mu_w = excess(cc, ww, fluid, wall)
        r = np.empty((n + 2, v.shape[1]))
        r[:n] = v[:n] + z[:, None] * pp + mu - ln_a[:n]
        r[n] = v[n] + mu_w - ln_a[n] if has_water else 0.0
        scale = (np.abs(z)[:, None] * cc).sum(axis=0) + np.abs(fixed) + 1e-12
        r[n + 1] = ((z[:, None] * cc).sum(axis=0) + fixed) / scale
        return r, mu, scale

    r, mu, scale = residual(u)
    ok, used, h = False, 0, 1e-6
    for used in range(1, max_iter + 1):  # noqa: B007
        norm = np.max(np.abs(r), axis=0)
        if float(np.max(norm)) < tol:
            ok = True
            break
        cc, ww, _ = unpack(u)
        jac = np.zeros((m, n + 2, n + 2))
        for j in range(n + 1):
            v = u.copy()
            v[j] += h
            c2, w2, _ = unpack(v)
            mu2, mu_w2 = excess(c2, w2, fluid, wall)
            jac[:, :n, j] = ((mu2 - mu) / h).T
            if has_water:
                jac[:, n, j] = (mu_w2 - excess(cc, ww, fluid, wall)[1]) / h
        jac[:, np.arange(n + 1), np.arange(n + 1)] += 1.0
        if not has_water:
            jac[:, n, :] = 0.0
            jac[:, n, n] = 1.0
        jac[:, :n, n + 1] = z[None, :]
        jac[:, n + 1, :n] = (z[:, None] * cc / scale[None, :]).T
        step = -np.linalg.solve(jac, r.T[:, :, None])[:, :, 0].T
        step *= np.minimum(1.0, cap / np.maximum(np.max(np.abs(step), axis=0),
                                                 1e-300))[None, :]
        frac = np.ones(m)
        for _ in range(40):
            trial = u + frac[None, :] * step
            c2, w2, _ = unpack(trial)
            r2, mu2, scale2 = residual(trial)
            worse = ~(np.max(np.abs(r2), axis=0) < norm)
            if fluid.hs:
                worse |= ~(packing(c2, w2, fluid, wall) < _PACKING_CAP)
            worse &= norm >= tol
            if not np.any(worse):
                break
            frac = np.where(worse, 0.5 * frac, frac)
        u, r, mu, scale = trial, r2, mu2, scale2
    cc, ww, pp = unpack(u)
    return Partition(cc, ww, pp * thermal, mu, ok, used)
