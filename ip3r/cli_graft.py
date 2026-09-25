"""AlphaFold fills, ``register(sub)``:

    python -m ip3r graft 8TKG [--mode full]   # fill every unresolved stretch
    python -m ip3r graft 8TKG --calibrate     # hide what other deposits miss
    python -m ip3r graft 8TKG --long          # windows 10-60 and islands (Round 7.9)
    python -m ip3r graft 7LHF                 # rat, through an alignment
"""

from __future__ import annotations

import numpy as np

__all__ = ["register"]


def _fills(model, chain: str) -> None:
    from .parameters import PARAMETERS as _P
    tol = _P.value("graft.join_tolerance")
    print(f"  {'stretch':14s} {'kind':6s} {'n':>3s} {'anchor':>6s} {'seams (Å)':>12s} "
          f"{'pLDDT':>5s} {'clash':>5s}")
    for f in model.fills:
        if chain and f.stretch.chain != chain:
            continue
        seams = "/".join(f"{j:.1f}" + ("!" if j > tol else "") for j in f.joins)
        print(f"  {f.stretch.label():14s} {f.stretch.kind:6s} {f.stretch.n_residues:3d} "
              f"{f.anchor_rmsd:6.2f} {seams:>12s} {f.plddt:5.0f} {f.clashes:5d}")
    for sk in model.skipped:
        if not chain or sk.stretch.chain == chain:
            print(f"  not filled {sk.stretch.label()}: {sk.reason}")
    for w in model.warnings():
        print(f"  ! {w}")


def _medians(trials) -> str:
    med = [float(np.median([getattr(t, k) for t in trials]))
           for k in ("rmsd_fill", "rmsd_line", "rmsd_global")]
    return (f"median fill {med[0]:.2f}, line {med[1]:.2f}, global {med[2]:.2f}; "
            f"fill beats the line in {sum(t.beats_line for t in trials)}/{len(trials)}")


def _calibrate(st, pred) -> None:
    from .io import loader
    from .io.registry import load_registry
    from .structure.graft import GraftRefusal, prediction_for
    from .structure.graft_calibration import calibrate
    others = []
    for e in load_registry():
        if e.pdb_id == st.name or not loader.is_local(e.pdb_id):
            continue
        o = loader.load(e.pdb_id)
        try:
            if prediction_for(o)[0].name == pred.name:
                others.append(o)
        except GraftRefusal:
            continue
    trials = calibrate(st, others, pred)
    print(f"\ncalibration: {len(trials)} stretches another deposit leaves "
          f"unresolved, hidden in {st.name} and filled (C-alpha RMSD, Å)")
    print(f"  {'stretch':14s} {'fill':>5s} {'line':>5s} {'global':>6s} {'pLDDT':>5s}")
    for t in trials:
        print(f"  {t.stretch.label():14s} {t.rmsd_fill:5.2f} {t.rmsd_line:5.2f} "
              f"{t.rmsd_global:6.2f} {t.plddt:5.0f}")
    if trials:
        print(f"  {_medians(trials)}")


def _long(st, pred, nm) -> None:
    from .parameters import PARAMETERS as _P
    from .structure.graft_calibration import island_trials, window_trials
    tol = _P.value("graft.join_tolerance")
    trials = window_trials(st, pred, numbering=nm)
    print(f"\nwindows: resolved stretches of each length hidden in {st.name} "
          "and filled (C-alpha RMSD, Å)")
    print(f"  {'length':>6s} {'n':>3s} {'fill':>5s} {'line':>5s} {'pLDDT':>5s} "
          f"{'seam max':>8s}")
    for n in sorted({t.stretch.n_residues for t in trials}):
        ts = [t for t in trials if t.stretch.n_residues == n]
        print(f"  {n:6d} {len(ts):3d} {np.median([t.rmsd_fill for t in ts]):5.2f} "
              f"{np.median([t.rmsd_line for t in ts]):5.2f} "
              f"{np.median([t.plddt_scored for t in ts]):5.0f} "
              f"{max(max(t.joins) for t in ts):8.2f}")
    isl = island_trials(st, pred, numbering=nm)
    print(f"\nislands: resolved runs between two gaps, hidden, the whole span "
          f"filled, scored on the island (seams over {tol:g} Å marked !)")
    print(f"  {'island':12s} {'span':>12s} {'fill':>5s} {'line':>5s} "
          f"{'pLDDT island/span':>17s}")
    for (a, b), t in isl:
        seam = "!" if max(t.joins) > tol else ""
        print(f"  {a:5d}-{b:<6d} {t.stretch.first:5d}-{t.stretch.last:<6d} "
              f"{t.rmsd_fill:5.1f} {t.rmsd_line:5.1f} {t.plddt_scored:8.0f}/"
              f"{t.plddt:<8.0f}{seam}")
    if isl:
        print(f"  {_medians([t for _, t in isl])}")


def _graft(args) -> int:
    from .io import loader
    from .structure.graft import GraftRefusal, fill_structure, prediction_for
    loader.ALLOW_FETCH = args.fetch
    st = loader.load(args.pdb)
    try:
        pred, nm = prediction_for(st)
    except GraftRefusal as exc:
        print(f"refused: {exc}")
        return 1
    model = fill_structure(st, args.mode, prediction=pred, numbering=nm)
    print(model.summary())
    _fills(model, args.chain)
    if args.calibrate:
        _calibrate(st, pred)
    if args.long:
        _long(st, pred, nm)
    return 0


def register(sub) -> None:
    p = sub.add_parser("graft", help="fill unresolved stretches from AlphaFold")
    p.add_argument("pdb")
    p.add_argument("--mode", choices=["gaps", "full"], default="gaps")
    p.add_argument("--chain", default="A", help="rows for one chain ('' = all)")
    p.add_argument("--calibrate", action="store_true",
                   help="hide stretches other deposits miss, fill, and score")
    p.add_argument("--long", action="store_true",
                   help="hide resolved windows (10-60) and islands, fill, and score")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_graft)
