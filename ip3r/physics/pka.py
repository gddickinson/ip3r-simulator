"""Protonation of the pore's ionisable groups: a Tanford-Kirkwood network.

The wall-charge model (:mod:`.pore_charge`) gives every lining Asp and Glu
a charge of -1 and every Lys and Arg +1. In a 4-6 A lumen that is the least
tested assumption left: eight carboxylates of the D2518/D2522 rings sit
within 10 A of each other, and like charges that close repel. Protonating
one ring member costs nothing in charge-charge energy and relieves the
other three.

**The model.** Every ionisable side chain within ``pka.site_radius`` of the
groups asked about is a site with two states: protonated (theta = 1) or
not. Its charge is ``q = theta + q0``, with ``q0 = -1`` for an acid and 0 for
a base. The free energy of a state, in kT, is

    G(theta) = sum_i theta_i ln10 (pH - pKa_i) + 1/2 sum_ij q_i W_ij q_j,

where pKa_i is the group's model-compound value (Thurlkill 2006's
pentapeptides; Arg from Fitch 2015), and

    W_ij = l_B / (eps(r_ij) r_ij) exp(-r_ij / lambda_D).

``l_B`` is the vacuum Bjerrum length, eps(r) is Mehler & Solmajer's
sigmoidal distance-dependent permittivity (1.3 at contact, water's 78.4 by
~20 A), and lambda_D is the bath's Debye length in water. An isolated site
therefore follows Henderson-Hasselbalch exactly. The sites are titrated
together by Metropolis Monte Carlo. For each site the estimator averages
the conditional probability of being protonated given the others at every
sweep, so a fraction of 1e-4 is resolved without waiting for 1e4 flips.

**What it leaves out** (and PROPKA, the second route in :mod:`.pka_propka`,
includes): desolvation, and hydrogen bonds to backbone and side chains.
Desolvation raises an acid's pKa and lowers a base's. So this network's
carboxylate fractions are an upper bound on the charge that
charge-charge coupling alone allows. The two routes answer different parts
of the question, and they are reported side by side.

Charge centres are :data:`.pore_charge.CENTRE_ATOMS` (His: the ring
nitrogens' midpoint). A residue whose side chain was not modelled is not a
site; it is counted, as in the wall charge.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .pore_charge import CENTRE_ATOMS, AVOGADRO

__all__ = ["Site", "Titration", "sites", "eps_sigmoidal", "debye_angstrom",
           "bjerrum_vacuum", "interactions", "exact", "monte_carlo", "ring_blocks",
           "titrate", "MODEL_KEYS", "Q0"]

#: Charge when deprotonated: acids -1, bases 0.
Q0 = {"ASP": -1.0, "GLU": -1.0, "HIS": 0.0, "LYS": 0.0, "ARG": 0.0}

#: The model-compound pKa of each group, as a registered parameter.
MODEL_KEYS = {"ASP": "pka.model_asp", "GLU": "pka.model_glu",
              "HIS": "pka.model_his", "LYS": "pka.model_lys",
              "ARG": "pka.model_arg"}

_CENTRES = {**CENTRE_ATOMS, "HIS": ("ND1", "NE2")}

_E = 1.602176634e-19          # C
_EPS0 = 8.8541878128e-12      # F/m
_KB = 1.380649e-23            # J/K


@dataclass(frozen=True)
class Site:
    """One titratable side chain."""

    chain: str
    res_seq: int
    res_name: str
    xyz: tuple[float, float, float]

    @property
    def key(self) -> tuple[str, int]:
        return (self.chain, self.res_seq)

    @property
    def q0(self) -> float:
        return Q0[self.res_name]

    def label(self) -> str:
        return f"{self.res_name}{self.res_seq}/{self.chain}"


@dataclass
class Titration:
    """Mean protonation of every site at one pH."""

    sites: list[Site]
    ph: float
    protonated: np.ndarray             # <theta_i>
    method: str
    meta: dict = field(default_factory=dict)

    @property
    def charge(self) -> np.ndarray:
        return self.protonated + np.array([s.q0 for s in self.sites])

    def apparent_pka(self) -> np.ndarray:
        """pH + log10(theta / (1 - theta)): the Henderson-Hasselbalch pKa
        that gives this site's protonation at this pH."""
        th = np.clip(self.protonated, 1e-12, 1.0 - 1e-12)
        return self.ph + np.log10(th / (1.0 - th))

    def charges(self) -> dict[tuple[str, int], float]:
        return {s.key: float(q) for s, q in zip(self.sites, self.charge)}

    def by_residue(self) -> list[tuple[str, float, float, float]]:
        """``(label, mean charge, lowest, highest apparent pKa)`` per residue
        number, over its copies."""
        rows: dict[str, list[int]] = {}
        for i, s in enumerate(self.sites):
            rows.setdefault(f"{s.res_name}{s.res_seq}", []).append(i)
        pk, q = self.apparent_pka(), self.charge
        return [(k, float(q[ii].mean()), float(pk[ii].min()),
                 float(pk[ii].max())) for k, ii in rows.items()]


# ------------------------------------------------------------------ sites
def sites(st: Structure, around: np.ndarray, radius: float | None = None
          ) -> tuple[list[Site], list[str]]:
    """Every modelled titratable side chain whose centre is within
    ``radius`` A of any point of ``around``; and the stubbed ionisable
    residues whose C-alpha is (``RESNAME RESNUM/chain``)."""
    radius = _P.value("pka.site_radius") if radius is None else radius
    around = np.atleast_2d(np.asarray(around, dtype=float))
    protein = ~st.hetero
    found, placed = [], set()
    for name, atoms in _CENTRES.items():
        m = protein & (st.res_name == name) & np.isin(st.atom_name, atoms)
        groups: dict[tuple[str, int], list[int]] = {}
        for i in np.flatnonzero(m):
            groups.setdefault((str(st.chain[i]), int(st.res_seq[i])), []).append(i)
        for (chain, seq), ii in groups.items():
            if len(ii) != len(atoms):
                continue
            placed.add((chain, seq))
            c = st.xyz[ii].astype(float).mean(axis=0)
            if np.min(np.linalg.norm(around - c, axis=1)) <= radius:
                found.append(Site(chain, seq, name, tuple(float(x) for x in c)))
    stubs = []
    ca = st.mask_ca() & np.isin(st.res_name, list(_CENTRES))
    for i in np.flatnonzero(ca):
        key = (str(st.chain[i]), int(st.res_seq[i]))
        if key not in placed and np.min(np.linalg.norm(
                around - st.xyz[i].astype(float), axis=1)) <= radius:
            stubs.append(f"{st.res_name[i]}{key[1]}/{key[0]}")
    found.sort(key=lambda s: (s.res_seq, s.chain))
    return found, sorted(stubs)


# ---------------------------------------------------------- electrostatics
def eps_sigmoidal(r: np.ndarray) -> np.ndarray:
    """Mehler & Solmajer's eps(r) = A + B / (1 + k exp(-lambda B r)),
    B = eps_water - A; r in A."""
    a = _P.value("pka.ms_a")
    b = _P.value("pka.eps_water") - a
    k = _P.value("pka.ms_k")
    lam = _P.value("pka.ms_lambda")
    return a + b / (1.0 + k * np.exp(-lam * b * np.asarray(r, dtype=float)))


def bjerrum_vacuum(temperature: float) -> float:
    """e^2 / (4 pi eps0 kT), A."""
    return _E ** 2 / (4.0 * np.pi * _EPS0 * _KB * temperature) * 1e10


def debye_angstrom(ionic_strength: float, temperature: float) -> float:
    """Debye length of a bath of ``ionic_strength`` M in water, A."""
    if ionic_strength <= 0.0:
        return np.inf
    eps = _P.value("pka.eps_water")
    kappa2 = (2.0 * AVOGADRO * _E ** 2 * ionic_strength * 1000.0
              / (eps * _EPS0 * _KB * temperature))
    return float(1e10 / np.sqrt(kappa2))


def interactions(xyz: np.ndarray, ionic_strength: float,
                 temperature: float | None = None,
                 eps: float | None = None) -> np.ndarray:
    """W_ij in kT per e^2 (zero diagonal); ``eps`` replaces the sigmoidal
    permittivity by a uniform one (the sweep that bounds it)."""
    temperature = (_P.value("permeation.temperature") if temperature is None
                   else temperature)
    xyz = np.asarray(xyz, dtype=float)
    r = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    np.fill_diagonal(r, np.inf)
    if np.any(r < 1.0):
        raise ValueError("two charge centres closer than 1 A")
    lb, ld = bjerrum_vacuum(temperature), debye_angstrom(ionic_strength,
                                                         temperature)
    e_r = eps_sigmoidal(r) if eps is None else float(eps)
    w = lb / (e_r * r) * np.exp(-r / ld)
    np.fill_diagonal(w, 0.0)
    return w


# -------------------------------------------------------------- titration
def _terms(pka, q0, ph):
    return np.log(10.0) * (ph - np.asarray(pka, dtype=float)), np.asarray(q0, float)


def exact(w: np.ndarray, pka, q0, ph: float) -> np.ndarray:
    """<theta_i> by enumerating all 2^n states (n <= 20): the calibration
    the Monte Carlo is held to."""
    n = len(q0)
    if n > 20:
        raise ValueError("exact enumeration is for n <= 20")
    a, q0 = _terms(pka, q0, ph)
    th = np.array(list(itertools.product((0.0, 1.0), repeat=n)))
    q = th + q0
    g = th @ a + 0.5 * np.einsum("si,ij,sj->s", q, w, q)
    p = np.exp(-(g - g.min()))
    return (p @ th) / p.sum()


def _block_energies(states, b, a, q0, h, th, w):
    """Energy of every state of block ``b`` given the rest, kT."""
    wbb = w[np.ix_(b, b)]
    h_ext = h[b] - wbb @ (th[b] + q0[b])
    q = states + q0[b]
    return states @ a[b] + q @ h_ext + 0.5 * np.einsum("si,ij,sj->s", q, wbb, q)


def monte_carlo(w: np.ndarray, pka, q0, ph: float, sweeps: int | None = None,
                seed: int | None = None, blocks=None) -> np.ndarray:
    """<theta_i> by Monte Carlo. Each sweep makes Metropolis single flips;
    pair flips for pairs coupled by more than ``pka.mc_pair_coupling`` kT;
    and heat-bath moves of each ``block`` (index arrays, <= 12 sites; the
    four copies of a ring), which draw the whole block from its 2^m states
    given the rest. The block move is what lets a strongly coupled ring pass
    between its degenerate states (two opposite members protonated, or the
    other two), which single and pair flips cannot do at eps 4.
    After the burn-in, each sweep adds every site's probability of being
    protonated given the rest (its block's, for a block member)."""
    sweeps = int(_P.value("pka.mc_sweeps") if sweeps is None else sweeps)
    burn = int(sweeps * _P.value("pka.mc_burn_fraction"))
    rng = np.random.default_rng(int(_P.value("pka.mc_seed") if seed is None
                                    else seed))
    a, q0 = _terms(pka, q0, ph)
    n = len(q0)
    blocks = [np.asarray(b, dtype=int) for b in (blocks or []) if len(b) > 1]
    if any(len(b) > 12 for b in blocks):
        raise ValueError("a block is enumerated; keep it to 12 sites")
    states = [np.array(list(itertools.product((0.0, 1.0), repeat=len(b))))
              for b in blocks]
    th = (a < 0).astype(float)                   # start from each site alone
    h = w @ (th + q0)
    iu, ju = np.nonzero(np.triu(np.abs(w) > _P.value("pka.mc_pair_coupling"), 1))
    acc = np.zeros(n)
    for sweep in range(sweeps):
        for i in rng.permutation(n):
            d = 1.0 - 2.0 * th[i]
            de = d * a[i] + d * h[i]
            if de <= 0.0 or rng.random() < np.exp(-de):
                th[i] += d
                h += d * w[:, i]
        for k in rng.permutation(len(iu)):
            i, j = iu[k], ju[k]
            di, dj = 1.0 - 2.0 * th[i], 1.0 - 2.0 * th[j]
            de = di * (a[i] + h[i]) + dj * (a[j] + h[j]) + di * dj * w[i, j]
            if de <= 0.0 or rng.random() < np.exp(-de):
                th[i] += di
                th[j] += dj
                h += di * w[:, i] + dj * w[:, j]
        marginals = []
        for k in rng.permutation(len(blocks)):
            b, st = blocks[k], states[k]
            g = _block_energies(st, b, a, q0, h, th, w)
            p = np.exp(-(g - g.min()))
            p /= p.sum()
            new = st[rng.choice(len(st), p=p)]
            h += w[:, b] @ (new - th[b])
            th[b] = new
            marginals.append((b, p @ st))
        if sweep >= burn:
            # G(theta_i = 1) - G(theta_i = 0) given the rest (W_ii = 0)
            cond = 1.0 / (1.0 + np.exp(np.clip(a + h, -700, 700)))
            for b, m in marginals:
                cond[b] = m
            acc += cond
    return acc / max(sweeps - burn, 1)


def ring_blocks(site_list: list[Site]) -> list[np.ndarray]:
    """The copies of each residue number (a C4 ring), as index arrays."""
    by: dict[tuple[str, int], list[int]] = {}
    for i, s in enumerate(site_list):
        by.setdefault((s.res_name, s.res_seq), []).append(i)
    return [np.array(v) for v in by.values() if 1 < len(v) <= 12]


def titrate(site_list: list[Site], ph: float, ionic_strength: float,
            method: str = "mc", eps: float | None = None) -> Titration:
    """The network's protonation at ``ph`` (``method`` ``mc`` or ``exact``;
    ``eps`` a uniform permittivity in place of the sigmoidal one)."""
    xyz = np.array([s.xyz for s in site_list])
    w = interactions(xyz, ionic_strength, eps=eps)
    pka = [_P.value(MODEL_KEYS[s.res_name]) for s in site_list]
    q0 = [s.q0 for s in site_list]
    theta = (exact(w, pka, q0, ph) if method == "exact" else
             monte_carlo(w, pka, q0, ph, blocks=ring_blocks(site_list)))
    return Titration(site_list, ph, theta, method,
                     meta={"ionic_strength_M": ionic_strength, "eps": eps,
                           "debye_A": debye_angstrom(
                               ionic_strength, _P.value("permeation.temperature")),
                           "sites": len(site_list)})
