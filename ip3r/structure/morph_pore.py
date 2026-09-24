"""The gate along a morph: every heavy atom interpolated, not only C-alpha.

:func:`ip3r.structure.transition.displaced_coords` carries side chains
rigidly with their C-alpha, so a frame's side chains are the *start*
deposit's rotamers wherever the backbone has gone. At the gate that is the
wrong answer by construction: the last frame would be 8TKF's backbone
lined by 8TKG's side chains, and its gate would not be 8TKF's gate.

Here each heavy atom of a basis residue is matched **by name** to the same
atom of the same residue in the end deposit (carried by the transition's
own superposition), and its offset from its C-alpha is interpolated
linearly while the C-alpha follows the morph:

    x_a(t) = site(t) + (1 - t) (x_a,start - CA_start) + t (x_a,end - CA_end)

so frame 0 is the start deposit's atoms and the last frame the end
deposit's, exactly. Between them the side chain is an interpolation too —
the same "path, not trajectory" caveat as the morph — and a linear offset
takes a chord when a side chain swings, which ``offset_error`` measures
(the worst shortening of an atom's distance from its C-alpha, Å).

Atoms with no partner (unresolved in one deposit, off the basis, ligands,
lipids) are left out: the profile is the protein that both deposits
resolve. The axis is re-found in every frame from the frame's own C-alphas,
and the gate is found by the deposit rule
(:func:`ip3r.structure.pore.constriction_indices`) on a span re-measured
per frame, so the endpoint frames are measured exactly as a deposit is.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.annotations import elements
from ..core.structure import Structure
from .morph import MorphTrajectory
from .pore import PROFILE_MARGIN, constriction_indices, radial_profile
from .symmetry import Frame, _basis, axis_by_superposition
from .transition import Transition

__all__ = ["AtomPath", "GatePath", "atom_path", "gate_path"]


@dataclass
class AtomPath:
    """The matched heavy atoms of a transition and their two offsets."""

    site: np.ndarray            # (k,) basis site each atom rides
    offset_start: np.ndarray    # (k, 3) atom - own C-alpha, start deposit
    offset_end: np.ndarray      # (k, 3) the same in the superposed end
    vdw: np.ndarray             # (k,) van der Waals radius, Å
    labels: np.ndarray          # (k,) residue label, e.g. "ILE2565"
    start_xyz: np.ndarray       # (k, 3) the start deposit's own atoms
    index: np.ndarray           # (k,) their atom indices in the start deposit
    meta: dict = field(default_factory=dict)

    def coords(self, sites: np.ndarray, t: float) -> np.ndarray:
        """Atoms at path fraction ``t`` on the C-alpha sites ``sites``."""
        return sites[self.site] + (1.0 - t) * self.offset_start + t * self.offset_end

    def rigid(self, sites: np.ndarray, start_sites: np.ndarray) -> np.ndarray:
        """The same atoms carried rigidly with their C-alpha (the viewer's
        :func:`displaced_coords`): start rotamers on the moved backbone."""
        return self.start_xyz + (sites - start_sites)[self.site]

    def offset_error(self, t: float, mask: np.ndarray | None = None) -> float:
        """Worst shortening of an atom–C-alpha distance below its
        interpolated length at ``t`` (Å): the chord a swinging side chain
        takes. ``mask`` restricts it to some atoms (the gate's lining)."""
        s, e = self.offset_start, self.offset_end
        if mask is not None:
            s, e = s[mask], e[mask]
        now = np.linalg.norm((1.0 - t) * s + t * e, axis=1)
        short = (1.0 - t) * np.linalg.norm(s, axis=1) + t * np.linalg.norm(e, axis=1) - now
        return float(short.max(initial=0.0))


def _atoms_by_name(st: Structure, chain: str, residues: set[int]) -> dict:
    sel = np.flatnonzero((st.chain == chain) & ~st.hetero & (st.element != "H")
                         & st.mask_protein())
    out = {}
    for i in sel:
        key = (int(st.res_seq[i]), str(st.atom_name[i]))
        if key[0] in residues:
            out.setdefault(key, i)          # first altloc wins
    return out


def atom_path(st_start: Structure, st_end: Structure, tr: Transition) -> AtomPath:
    """Match every heavy atom of the basis residues between the two deposits."""
    if st_start.name != tr.start_id or st_end.name != tr.end_id:
        raise ValueError(f"{st_start.name}->{st_end.name} is not "
                         f"{tr.start_id}->{tr.end_id}")
    r, t = tr.meta["end_transform"]
    m = len(tr.residues)
    pos = {int(res): i for i, res in enumerate(tr.residues)}
    basis = set(pos)
    vdw_all = st_start.vdw_radii()
    site, i_start, xyz_end, unmatched = [], [], [], 0
    for k, (ca, cb) in enumerate(zip(tr.chains_start, tr.chains_end)):
        a = _atoms_by_name(st_start, ca, basis)
        b = _atoms_by_name(st_end, cb, basis)
        unmatched += len(a.keys() ^ b.keys())
        for key in sorted(a.keys() & b.keys()):
            site.append(k * m + pos[key[0]])
            i_start.append(a[key])
            xyz_end.append(st_end.xyz[b[key]])
    site = np.array(site)
    i_start = np.array(i_start)
    start_xyz = st_start.xyz[i_start].astype(np.float64)
    end_xyz = np.asarray(xyz_end, np.float64) @ r.T + t
    labels = np.char.add(st_start.res_name[i_start].astype(str),
                         st_start.res_seq[i_start].astype(str))
    return AtomPath(site, start_xyz - tr.start[site], end_xyz - tr.end[site],
                    vdw_all[i_start], labels, start_xyz, i_start,
                    {"n_atoms": len(site), "unmatched": unmatched})


def _frame_of(sites: np.ndarray, tr: Transition) -> Frame:
    """The four-fold frame of one morph frame, from its own C-alphas."""
    m = len(tr.residues)
    ca = {c: dict(zip(tr.residues.tolist(), sites[k * m:(k + 1) * m]))
          for k, c in enumerate(tr.chains_start)}
    axis, centre, _ = axis_by_superposition(ca, tr.n_subunits)
    if axis @ tr.frame.axis < 0:
        axis = -axis
    return Frame(axis, centre, _basis(axis), list(tr.chains_start))


def _span(sites_f: np.ndarray, pore_sites: np.ndarray) -> tuple[float, float]:
    z = sites_f[pore_sites, 2]
    return float(np.percentile(z, 2)), float(np.percentile(z, 98))


@dataclass
class GatePath:
    """The gate and filter measured on every frame of a morph."""

    fraction: np.ndarray        # (n,) path fraction of each frame
    gate: np.ndarray            # (n,) gate radius, atoms interpolated (Å)
    gate_z: np.ndarray          # (n,) its axial position in the frame (Å)
    filter: np.ndarray          # (n,) filter radius, atoms interpolated (Å)
    rigid_gate: np.ndarray      # (n,) gate radius with side chains rigid (Å)
    lining: list                # per frame, residues at the gate
    offset_error: np.ndarray    # (n,) worst side-chain chord, all atoms (Å)
    lining_error: np.ndarray    # (n,) the same over the gate's lining residues
    profiles: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def half_open(self) -> float:
        """First path fraction at which the gate has covered half its
        change (linear between frames); NaN if it never does."""
        g0, g1 = self.gate[0], self.gate[-1]
        half = 0.5 * (g0 + g1)
        s = np.sign(g1 - g0) or 1.0
        above = s * (self.gate - half) >= 0
        i = int(np.argmax(above))
        if not above[i]:
            return float("nan")
        if i == 0:
            return float(self.fraction[0])
        f0, f1 = self.fraction[i - 1], self.fraction[i]
        a, b = self.gate[i - 1], self.gate[i]
        return float(f0 + (half - a) / (b - a) * (f1 - f0))

    def overshoot(self) -> float:
        """How far the gate leaves the band between its endpoints (Å)."""
        lo, hi = sorted((self.gate[0], self.gate[-1]))
        return float(max(lo - self.gate.min(), self.gate.max() - hi, 0.0))


def _lining(f: np.ndarray, labels: np.ndarray, z: float, radius: float) -> list[str]:
    from ..parameters import PARAMETERS as _P
    slab, tol = _P.value("pore.lining_slab"), _P.value("pore.lining_tol")
    r = np.hypot(f[:, 0], f[:, 1])
    sel = (np.abs(f[:, 2] - z) <= slab) & (r <= radius + tol)
    return sorted(set(labels[sel].tolist()),
                  key=lambda s: (int("".join(c for c in s if c.isdigit())), s))


def _measure(xyz, frame, span, vdw):
    f = frame.to_frame(xyz)
    prof = radial_profile(f, vdw, span[0] - PROFILE_MARGIN, span[1] + PROFILE_MARGIN)
    return f, prof, constriction_indices(prof, span)


def gate_path(path: AtomPath, tr: Transition, mt: MorphTrajectory,
              keep_profiles: bool = False) -> GatePath:
    """Measure the pore on every morph frame, both ways."""
    pore = [e for e in elements(tr.paralog) if e.name == "channel"][0]
    m = len(tr.residues)
    pore_sites = np.flatnonzero((tr.residues >= pore.start) & (tr.residues <= pore.end))
    if not len(pore_sites):
        raise ValueError(f"no {pore.name} residues in the basis")
    n = len(mt.frames)
    frac = np.linspace(0.0, 1.0, n)
    rows = {k: np.full(n, np.nan) for k in ("gate", "gate_z", "filter", "rigid", "chord",
                                         "lchord")}
    lining, profiles = [], []
    for i, sites in enumerate(mt.frames):
        frame = _frame_of(sites, tr)
        span = _span(frame.to_frame(sites[:m]), pore_sites)
        f, prof, idx = _measure(path.coords(sites, frac[i]), frame, span, path.vdw)
        _, rprof, ridx = _measure(path.rigid(sites, mt.frames[0]), frame, span, path.vdw)
        if "gate" in idx:
            g = idx["gate"]
            rows["gate"][i], rows["gate_z"][i] = prof.r_min[g], prof.z[g]
            lining.append(_lining(f, path.labels, prof.z[g], prof.r_min[g]))
            rows["lchord"][i] = path.offset_error(
                frac[i], np.isin(path.labels, lining[-1]))
        else:
            lining.append([])
        if "filter" in idx:
            rows["filter"][i] = prof.r_min[idx["filter"]]
        if "gate" in ridx:
            rows["rigid"][i] = rprof.r_min[ridx["gate"]]
        rows["chord"][i] = path.offset_error(frac[i])
        if keep_profiles:
            profiles.append(prof)
    return GatePath(frac, rows["gate"], rows["gate_z"], rows["filter"], rows["rigid"],
                    lining, rows["chord"], rows["lchord"], profiles,
                    {"start": tr.start_id, "end": tr.end_id, "method": mt.method,
                     **path.meta})

