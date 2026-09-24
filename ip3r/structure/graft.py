"""Unresolved stretches of a deposit, filled from an AlphaFold model, seams kept.

Every deposit leaves 344-576 residues per ITPR3 subunit unbuilt: the
disordered loops, the ends, stretches the map could not trace. This places the
AlphaFold model's residues there, **one stretch at a time**, and keeps the
seams visible.

**One stretch at a time, anchored locally.** Each stretch is placed by a rigid
fit of the prediction's C-alphas onto the deposit's on the residues either
side of it (:data:`graft.anchor_window`). A global fit would spread the two
models' disagreement elsewhere into every join. An internal gap is anchored on
both sides, so it is interpolated. A terminus has one side only, so it is
extrapolated, cantilevered from its seam, and filled only when asked
(``mode="full"``).

**Anchors must be evidence.** An anchor residue is resolved in the deposit,
is not a backbone+CB stub (6DQN's 1434-1546 segment was built without a
sequence register), and has the same amino acid in the prediction at the same
number. A stretch without enough anchors is skipped and the reason is kept.
It is never fitted on whatever happens to be nearby.

**What is measured, per fill.** The anchor RMSD. The C-alpha distance across
each seam, which is 3.8 A for a real peptide bond, so a broken join shows as
a number. The mean pLDDT. The residues whose heavy atoms interpenetrate the
deposit, *including neighbouring subunits*: the prediction is a monomer and
knows nothing about them.

**The prediction must be in the deposit's numbering.** It is chosen by
measurement (:func:`prediction_for`), never by name. AlphaFold DB holds
isoforms, not the canonical sequence, for ITPR1 and rat ITPR1, and 7LHF and
9YKK are refused with the identities that refuse them.

**The fill is not the deposit.** A :class:`FilledModel` holds only the
predicted atoms, as a separate structure. Every measurement in this
application still runs on the deposit alone; the fill is drawn beside it.
:meth:`FilledModel.place` re-fits each stretch to new deposit coordinates,
so a fill follows a morph or a mode frame on its own anchors.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .numbering import chain_residues
from .symmetry import kabsch

__all__ = ["FILL_MODES", "GraftRefusal", "Stretch", "Fill", "Skip",
           "FilledModel", "prediction_for", "unresolved", "fill_stretches",
           "fill_structure", "placed_ca"]

#: (key, label, description). The default is the deposit alone: a view that
#: silently contained prediction is the outcome to avoid.
FILL_MODES = (
    ("none", "Deposited only", "Exactly what the experiment resolved."),
    ("gaps", "+ AlphaFold gaps", "Unresolved stretches inside the deposited "
     "range, each anchored on resolved residues at both ends (interpolated)."),
    ("full", "+ AlphaFold gaps and ends", "Also the unresolved termini, "
     "anchored at one end only (extrapolated)."),
)


class GraftRefusal(ValueError):
    """No downloaded prediction is in this deposit's numbering."""


@dataclass(frozen=True)
class Stretch:
    chain: str
    first: int
    last: int
    kind: str                     # "gap" | "n_term" | "c_term"

    @property
    def n_residues(self) -> int:
        return self.last - self.first + 1

    def label(self) -> str:
        return f"{self.chain}:{self.first}-{self.last}"


@dataclass
class Fill:
    stretch: Stretch
    pred_atoms: np.ndarray        # indices into the prediction
    anchor_dep: np.ndarray        # deposit atom indices (C-alpha) of the anchors
    anchor_pred: np.ndarray       # prediction C-alpha coordinates of the same
    seam_dep: np.ndarray          # deposit C-alpha index at each seam
    seam_pred: np.ndarray         # prediction C-alpha coordinates at each seam
    anchor_rmsd: float
    joins: tuple                  # C-alpha distance across each seam, A
    plddt: float                  # mean over the fill's C-alphas
    confident: float              # fraction of C-alphas >= graft.plddt_confident
    clashes: int = 0              # residues interpenetrating the deposit
    seam_fill: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    """Index into :attr:`FilledModel.atoms` of the fill's C-alpha at each seam."""

    @property
    def joined(self) -> bool:
        tol = _P.value("graft.join_tolerance")
        return all(j <= tol for j in self.joins)


@dataclass(frozen=True)
class Skip:
    stretch: Stretch
    reason: str


@dataclass
class FilledModel:
    deposit: str
    prediction: str
    identity: float               # prediction vs deposit, by number
    mode: str
    fills: list[Fill]
    skipped: list[Skip]
    atoms: Structure              # the predicted atoms only, placed
    _rows: list = field(default_factory=list, repr=False)

    @property
    def n_residues(self) -> int:
        return sum(f.stretch.n_residues for f in self.fills)

    def place(self, deposit_xyz: np.ndarray) -> np.ndarray:
        """Coordinates of :attr:`atoms` re-fitted onto new deposit coordinates."""
        pred = self._pred_xyz
        out = np.empty((self.atoms.n_atoms, 3), np.float64)
        dep = np.asarray(deposit_xyz, np.float64)
        for f, (a, b) in zip(self.fills, self._rows):
            r, t = kabsch(f.anchor_pred, dep[f.anchor_dep])
            out[a:b] = pred[f.pred_atoms] @ r.T + t
        return out

    def seam_segments(self, deposit_xyz: np.ndarray, fill_xyz: np.ndarray | None = None):
        """``(deposit end, fill end, joined)`` for every seam, in given coordinates."""
        fill_xyz = self.atoms.xyz if fill_xyz is None else fill_xyz
        dep = np.array([deposit_xyz[d] for f in self.fills for d in f.seam_dep]).reshape(-1, 3)
        end = np.array([fill_xyz[i] for f in self.fills for i in f.seam_fill]).reshape(-1, 3)
        gap = np.linalg.norm(dep - end, axis=1) if len(dep) else np.zeros(0)
        return dep, end, gap <= _P.value("graft.join_tolerance")

    def summary(self) -> str:
        broken = sum(not f.joined for f in self.fills)
        clash = sum(f.clashes for f in self.fills)
        ca = self.atoms.mask_ca()
        plddt = float(self.atoms.b_factor[ca].mean()) if ca.any() else float("nan")
        return (f"{self.deposit} + {self.prediction} ({self.identity:.1%} identical "
                f"by number): {len(self.fills)} stretches, {self.n_residues} "
                f"residues filled, mean pLDDT {plddt:.0f}; {broken} broken "
                f"seams; {clash} residues clash with the deposit; "
                f"{len(self.skipped)} stretches not filled")

    def warnings(self) -> list[str]:
        out = ["the filled residues are an AlphaFold MODEL placed on local "
               "anchors; no experiment resolves them, and nothing in this "
               "application measures on them"]
        if any(not f.joined for f in self.fills):
            out.append("some seams do not close (drawn red): the prediction's "
                       "loop does not reach both anchors in this conformation")
        if any(f.clashes for f in self.fills):
            out.append("some filled residues pass through the deposit (often "
                       "a neighbouring subunit, which a monomer prediction "
                       "cannot know about)")
        return out


# ------------------------------------------------------------------ choosing

def _identity(dep: dict, pred: dict) -> tuple[float, int]:
    scored = [r for r, (aa, stub) in dep.items() if not stub and r in pred]
    same = sum(dep[r][0] == pred[r][0] for r in scored)
    return (same / len(scored) if scored else 0.0), len(scored)


def _main_chain(st: Structure) -> str:
    return max(st.chains, key=lambda c: int((st.mask_ca() & (st.chain == c)).sum()))


def prediction_for(st: Structure, paths=None) -> tuple[Structure, float]:
    """The downloaded model in this deposit's numbering, and its identity.

    Identity by residue number over unstubbed residues, against the bar the
    variant painting uses (``numbering.min_identity``). Raises
    :class:`GraftRefusal` naming every model's identity when none clears it.
    """
    from ..io.predictions import load_prediction, local_predictions

    paths = local_predictions() if paths is None else paths
    if not paths:
        raise GraftRefusal("no AlphaFold model is downloaded "
                           "(run `python -m ip3r fetch`)")
    dep = chain_residues(st, _main_chain(st))
    scored = []
    for p in paths:
        pred = load_prediction(str(p))
        ident, _ = _identity(dep, chain_residues(pred, pred.chains[0]))
        scored.append((ident, pred))
    ident, best = max(scored, key=lambda s: s[0])
    if ident < _P.value("numbering.min_identity"):
        table = ", ".join(f"{p.name} {i:.1%}" for i, p in scored)
        raise GraftRefusal(
            f"{st.name}: no AlphaFold model is in its numbering ({table}; "
            f"the bar is {_P.value('numbering.min_identity'):.0%})")
    return best, ident


# ----------------------------------------------------------------- stretches

def unresolved(st: Structure, chain: str, pred_last: int,
               termini: bool = False) -> list[Stretch]:
    """Stretches of ``chain`` with no C-alpha, internal and (optionally) ends."""
    res = sorted(chain_residues(st, chain))
    if not res:
        return []
    nums = np.array(res)
    out = [Stretch(chain, int(nums[i]) + 1, int(nums[i + 1]) - 1, "gap")
           for i in np.flatnonzero(np.diff(nums) > 1)]
    if termini:
        if nums[0] > 1:
            out.insert(0, Stretch(chain, 1, int(nums[0]) - 1, "n_term"))
        if nums[-1] < pred_last:
            out.append(Stretch(chain, int(nums[-1]) + 1, pred_last, "c_term"))
    return out


def _ca_index(st: Structure, chain: str) -> dict[int, int]:
    m = st.mask_ca() & (st.chain == chain) & ~st.hetero
    return {int(r): int(i) for r, i in zip(st.res_seq[m], np.flatnonzero(m))}


def _anchors(s: Stretch, dep_res: dict, pred_res: dict) -> tuple[list, list]:
    """Anchor residues before and after a stretch (a terminus looks twice as far)."""
    w = int(_P.value("graft.anchor_window"))

    def ok(r):
        return (r in dep_res and not dep_res[r][1] and r in pred_res
                and dep_res[r][0] == pred_res[r][0])
    if s.kind == "n_term":
        return [], [r for r in range(s.last + 1, s.last + 1 + 2 * w) if ok(r)]
    if s.kind == "c_term":
        return [r for r in range(s.first - 2 * w, s.first) if ok(r)], []
    return ([r for r in range(s.first - w, s.first) if ok(r)],
            [r for r in range(s.last + 1, s.last + 1 + w) if ok(r)])


def _seams(s: Stretch) -> list[tuple[int, int]]:
    """(deposit residue, prediction residue) across each seam."""
    out = []
    if s.kind != "n_term":
        out.append((s.first - 1, s.first))
    if s.kind != "c_term":
        out.append((s.last + 1, s.last))
    return out


# -------------------------------------------------------------------- filling

def fill_stretches(st: Structure, pred: Structure, stretches: list[Stretch],
                   identity: float = float("nan"), mode: str = "gaps") -> FilledModel:
    """Place the prediction's residues for each stretch; skip with a reason."""
    from scipy.spatial import cKDTree

    k = int(_P.value("graft.min_anchor"))
    pchain = pred.chains[0]
    pred_res = chain_residues(pred, pchain)
    pred_ca = _ca_index(pred, pchain)
    pxyz = pred.xyz.astype(np.float64)
    heavy = (st.element != "H") & ~np.isin(st.res_name, ["HOH", "WAT"])
    tree = cKDTree(st.xyz[heavy])
    heavy_idx = np.flatnonzero(heavy)
    clash_d = _P.value("graft.clash_distance")
    conf = _P.value("graft.plddt_confident")

    fills, skipped, parts, rows, start = [], [], [], [], 0
    by_chain: dict[str, tuple] = {}
    for s in stretches:
        if s.chain not in by_chain:
            by_chain[s.chain] = (chain_residues(st, s.chain), _ca_index(st, s.chain))
        dep_res, dep_ca = by_chain[s.chain]
        before, after = _anchors(s, dep_res, pred_res)
        need_b = 0 if s.kind == "n_term" else (2 * k if s.kind == "c_term" else k)
        need_a = 0 if s.kind == "c_term" else (2 * k if s.kind == "n_term" else k)
        if len(before) < need_b or len(after) < need_a:
            skipped.append(Skip(s, f"anchors {len(before)} before / {len(after)} "
                                   f"after, need {need_b} / {need_a}"))
            continue
        inside = (pred.chain == pchain) & (pred.res_seq >= s.first) \
            & (pred.res_seq <= s.last) & ~pred.hetero
        if not inside.any():
            skipped.append(Skip(s, "the prediction does not model these residues"))
            continue
        anchor = before + after
        a_dep = np.array([dep_ca[r] for r in anchor])
        a_pred = pxyz[[pred_ca[r] for r in anchor]]
        r, t = kabsch(a_pred, st.xyz[a_dep])
        rmsd = float(np.sqrt((((a_pred @ r.T + t) - st.xyz[a_dep]) ** 2).sum(1).mean()))
        seams = [(dep_ca[d], pxyz[pred_ca[p]]) for d, p in _seams(s)
                 if d in dep_ca and p in pred_ca]
        joins = tuple(float(np.linalg.norm(q @ r.T + t - st.xyz[d])) for d, q in seams)
        idx = np.flatnonzero(inside)
        placed = pxyz[idx] @ r.T + t
        # Clashes with the deposit, the residues next to its own seams excepted.
        own = (st.chain == s.chain) & (st.res_seq >= s.first - 1) & (st.res_seq <= s.last + 1)
        hv = pred.element[idx] != "H"
        clash_res = {int(rs) for rs, h in zip(pred.res_seq[idx][hv],
                                              tree.query_ball_point(placed[hv], clash_d))
                     if h and not own[heavy_idx[h]].all()}
        clashing = len(clash_res)
        ca = pred.atom_name[idx] == "CA"
        pl = pred.b_factor[idx][ca].astype(float)
        fills.append(Fill(s, idx, a_dep, a_pred, np.array([d for d, _ in seams]),
                          np.array([q for _, q in seams]), rmsd, joins,
                          float(pl.mean()), float((pl >= conf).mean()), clashing))
        part = pred.subset(inside).copy_with_coords(placed.astype(np.float32))
        part.chain = np.full(part.n_atoms, s.chain, dtype=st.chain.dtype)
        fills[-1].seam_fill = np.array(
            [start + int(np.flatnonzero((part.res_seq == q) & (part.atom_name == "CA"))[0])
             for (d, q) in _seams(s) if d in dep_ca and q in pred_ca], int)
        parts.append(part)
        rows.append((start, start + part.n_atoms))
        start += part.n_atoms

    atoms = _concat(parts, f"{st.name}+{pred.name}", pred)
    model = FilledModel(st.name, pred.name, identity, mode, fills, skipped, atoms, rows)
    model._pred_xyz = pxyz
    return model


def _concat(parts: list[Structure], name: str, template: Structure) -> Structure:
    if not parts:
        return template.subset(np.zeros(template.n_atoms, bool), name=name)
    arrays = {f: np.concatenate([getattr(p, f) for p in parts])
              for f in Structure._ARRAY_FIELDS}
    built = Structure(name=name, meta={"source": "predicted"}, **arrays)
    built._build_residue_index()
    return built


def fill_structure(st: Structure, mode: str = "gaps", prediction=None,
                   chains=None) -> FilledModel:
    """Fill every chain of a deposit (``mode`` from :data:`FILL_MODES`).

    Raises :class:`GraftRefusal` when no downloaded model is in the deposit's
    numbering (or ``mode`` is ``none``, which has nothing to build).
    """
    if mode not in {k for k, _, _ in FILL_MODES} or mode == "none":
        raise ValueError(f"mode must be 'gaps' or 'full', not {mode!r}")
    if prediction is None:
        prediction, identity = prediction_for(st)
    else:
        identity, _ = _identity(chain_residues(st, _main_chain(st)),
                                chain_residues(prediction, prediction.chains[0]))
    pred_res = chain_residues(prediction, prediction.chains[0])
    last = max(pred_res)
    # Every chain that is itself in the prediction's numbering (a ligand-only
    # or foreign chain is not a subunit to fill).
    bar = _P.value("numbering.min_identity")
    chains = chains or [c for c in st.chains
                        if _identity(chain_residues(st, c), pred_res)[0] >= bar]
    stretches = [s for c in chains
                 for s in unresolved(st, c, last, termini=(mode == "full"))]
    return fill_stretches(st, prediction, stretches, identity, mode)


def placed_ca(model: FilledModel, fill: Fill, deposit_xyz=None) -> dict[int, np.ndarray]:
    """``resnum -> placed C-alpha`` for one fill (on the deposit, or new coords)."""
    a, b = model._rows[model.fills.index(fill)]
    xyz = model.atoms.xyz if deposit_xyz is None else model.place(deposit_xyz)
    sub = slice(a, b)
    ca = model.atoms.atom_name[sub] == "CA"
    return {int(r): x for r, x in zip(model.atoms.res_seq[sub][ca], xyz[sub][ca])}
