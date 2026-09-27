"""Selectivity in 3-D (Round 7.19), ``register(sub)``:

    python -m ip3r sel3d                 # 9HEO, then 8TKF
    python -m ip3r sel3d 9HEO --mutants  # Xu 2006's five mutants on 9HEO
    python -m ip3r sel3d 7T3T --spacing 1.0
    python -m ip3r gate 9HEO             # Round 7.21: the gate widened
    python -m ip3r reversal 9HEO 8TKF    # Round 7.23: at bi-ionic reversal
    python -m ip3r reversal 8TKF --lumen neutral "pb + csc"   # Round 7.24
"""

from __future__ import annotations

import math

__all__ = ["register"]

_DEFAULT = ("9HEO", "8TKF")


def _deposit(r) -> None:
    flag = "" if r.converged else "  (n.c.)"
    print(f"\n{r.name}: formal wall {r.net:+.0f} e, grid {r.spacing:g} A; "
          f"measured P_Ca:P_K {r.measured:g}{flag}")
    if r.unreached:
        print(f"  groups reaching no lumen voxel (left out): {', '.join(r.unreached)}")
    print("  1-D, same bath: " + ", ".join(f"{k} {v:.2f}" for k, v in r.one_d.items()))
    print(f"  {'reading':12s} {'P_Ca:P_K':>8s} {'P_Cl:P_K':>8s}   resistance at "
          f"filter / gate (K+, Ca2+)   highest Ca2+, K+ (M)")
    for label, x in r.readings.items():
        f, g = x.shares["filter"], x.shares["gate"]
        print(f"  {label:12s} {x.ratio:8.2f} {x.pcl_pk:8.2f}   "
              f"K+ {f['K+']:4.0%} / {g['K+']:4.0%}   Ca2+ {f['Ca2+']:4.0%} / "
              f"{g['Ca2+']:4.0%}      {x.peak['Ca2+']:7.2f} {x.peak['K+']:6.2f}"
              f"{'' if x.converged else '  (n.c.)'}")


def _sel3d(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import selectivity3d as s3
    loader.ALLOW_FETCH = args.fetch
    print("P_Ca:P_K in 3-D linear response: P_i = D_i x the Boltzmann-weighted "
          "Laplace conductance, symmetric family KCl + "
          f"{_P.value('selectivity.ryr1_cacl2_lumen') * 1e3:g} mM CaCl2. "
          "Not a reversal potential: the 1-D readings beside it give the "
          "translation.")
    kw = {} if args.wall_volume is None else {"wall_volume": args.wall_volume}
    for pdb in args.pdb or _DEFAULT:
        st = loader.load(pdb)
        if args.mutants:
            rows = s3.mutant_panel(st, spacing=args.spacing, progress=lambda i, n, name: print(
                f"  [{i + 1}/{n}] {name}", flush=True), **kw)
            _deposit(rows[0][2])
            heads = [h for h in s3.READINGS if h != "neutral"] + ["1-D csc"]
            wt = rows[0][2]
            print("\n  Xu 2006's mutants, P_Ca:P_K over the wild type:")
            print(f"  {'':10s}{'measured':>9s}" + "".join(f"{h:>13s}" for h in heads))
            err = dict.fromkeys(heads, 0.0)
            for name, pca, r in rows[1:]:
                cells = ""
                for h in heads:
                    x = ((r.one_d[h] / wt.one_d[h]) if h.startswith("1-D")
                         else r.ratio(h) / wt.ratio(h))
                    err[h] += abs(math.log(x / (pca / rows[0][1])))
                    cells += f"{'x':>8s}{x:5.2f}"
                print(f"  {name:10s}{'x':>5s}{pca / rows[0][1]:4.2f}{cells}"
                      f"{'' if r.converged else '  (n.c.)'}")
            print(f"  {'sum |ln|':19s}" + "".join(f"{err[h]:13.2f}" for h in heads))
            continue
        _deposit(s3.selectivity_3d(st, spacing=args.spacing, **kw))
    return 0


def _gate(args) -> int:
    from .io import loader
    from .physics.gate_geometry import READINGS, gate_scan
    loader.ALLOW_FETCH = args.fetch
    print("The gate widened radially (gate.widen_half_width taper), then "
          "Round 7.6's neutral K+ conductance and Round 7.19's P_Ca:P_K "
          "(linear response); shares are Ca2+'s, at the deposited gate and filter.")
    for pdb in args.pdb or ("9HEO", "8TKF"):
        st = loader.load(pdb)
        scan = gate_scan(st, widths=args.delta, spacing=args.spacing,
                         progress=lambda i, n, d: print(f"  [{i + 1}/{n}] +{d:g} A", flush=True))
        print(f"\n{pdb}: measured P_Ca:P_K {scan.measured:g}")
        print(f"  {'+A':>4s} {'gate r':>6s} {'narrowest':>9s} {'g K+ pS':>8s}   "
              + "   ".join(f"{k:>11s} (gate/filter)" for k in READINGS))
        for r in scan.rows:
            cells = "   ".join(f"{r.ratios[k]:11.2f} ({r.gate_share[k]:4.0%}/"
                               f"{r.filter_share[k]:4.0%})" for k in READINGS)
            print(f"  {r.delta:4g} {r.gate_radius:6.2f} {r.narrowest:>9s} "
                  f"{r.conductance:8.0f}   {cells}{'' if r.converged else '  (n.c.)'}")
    return 0


def _rev_deposit(r) -> None:
    ruler = ("Xu's Eq. 1 (P_Cl taken 0)" if r.ryr
             else "Vais's Eq. 1 (P_Cl:P_K from the KCl/NMDG experiment)")
    print(f"\n{r.name}: formal wall {r.net:+.0f} e, grid {r.spacing:g} A; "
          f"measured P_Ca:P_K {r.measured:g}; ruler {ruler}"
          f"{'' if r.converged else '  (n.c.)'}")
    if r.unreached:
        print(f"  groups reaching no lumen voxel (left out): {', '.join(r.unreached)}")
    if r.one_d:
        print("  1-D, same wall: " + ", ".join(f"{k} {v:.2f}" for k, v in r.one_d.items()))
    extra = "  own P_Cl" if r.ryr else ""
    print(f"  {'reading':10s} {'V_rev mV':>18s} {'P_Ca:P_K':>9s} {'linear':>7s} "
          f"{'P_Cl:P_K':>9s}{extra}   highest Ca2+ (M)   solves")
    for label, x in r.readings.items():
        v = " / ".join(f"{x.v[k] * 1e3:+6.2f}" for k in sorted(x.v, reverse=True))
        lin = r.linear.get(label, math.nan)
        own = f"  {x.pca_pk_model_cl:8.2f}" if r.ryr else ""
        print(f"  {label:10s} {v:>18s} {x.pca_pk:9.2f} {lin:7.2f} "
              f"{x.pcl_pk:9.3f}{own}   {x.peak_ca:16.2f}   {x.solves:6d}"
              f"{'' if x.converged else '  (n.c.)'}")


def _reversal(args) -> int:
    from .io import loader
    from .physics import reversal3d as r3
    loader.ALLOW_FETCH = args.fetch
    print("P_Ca:P_K at bi-ionic reversal in 3-D: steady Poisson-Nernst-Planck "
          "on the voxels (Gummel), the root of the net current, read by the "
          "family's GHK ruler; 'linear' is Round 7.19's reading on the same grid.")
    for pdb in args.pdb or _DEFAULT:
        st = loader.load(pdb)
        if args.mutants:
            rows = r3.mutant_panel(st, spacing=args.spacing, progress=lambda i, n, name: print(
                f"  [{i + 1}/{n}] {name}", flush=True))
            _rev_deposit(rows[0][2])
            wt, heads = rows[0], list(rows[0][2].readings)
            print("\n  Xu 2006's mutants at reversal, P_Ca:P_K over the wild type:")
            print(f"  {'':10s}{'measured':>9s}" + "".join(f"{h:>13s}" for h in heads))
            err = dict.fromkeys(heads, 0.0)
            for name, pca, r in rows[1:]:
                cells = ""
                for h in heads:
                    x = r.readings[h].pca_pk / wt[2].readings[h].pca_pk
                    err[h] += abs(math.log(x / (pca / wt[1])))
                    cells += f"{'x':>8s}{x:5.2f}"
                print(f"  {name:10s}{'x':>5s}{pca / wt[1]:4.2f}{cells}"
                      f"{'' if r.converged else '  (n.c.)'}")
            print(f"  {'sum |ln|':19s}" + "".join(f"{err[h]:13.2f}" for h in heads))
            continue
        if args.lumen is not None:
            args.lumen = args.lumen or ["pb + csc"]
            _rev_lumen(st, args)
            continue
        if args.scale:
            print(f"\n{pdb}: the pb wall's charge scaled (point ions), at the "
                  "family's reversal; IP3R measured P_Cl:P_K 0.27, P_Ca:P_K 15.2 "
                  "(Vais 2010); RyR1 P_Ca:P_K 7.0 (Xu 2006, P_Cl taken 0)")
            print(f"  {'scale':>6s} {'P_Cl:P_K':>9s} {'P_Ca:P_K':>9s}")
            for f, x in r3.charge_scan(st, args.scale, spacing=args.spacing):
                print(f"  {f:6g} {x.pcl_pk:9.3f} {x.pca_pk:9.2f}"
                      f"{'' if x.converged else '  (n.c.)'}")
            continue
        _rev_deposit(r3.reversal_3d(st, spacing=args.spacing))
    return 0


def _rev_lumen(st, args) -> None:
    """Round 7.24: each ion on the lumen at the experiment's reversal."""
    from .parameters import PARAMETERS as _P
    from .physics.lumen_reversal import reversal_lumen
    from .structure.channel import measure_channel
    s = measure_channel(st)
    w = _P.value("lumen.constriction_half_width")
    for reading in args.lumen:
        r = reversal_lumen(st, reading, args.experiment, s, spacing=args.spacing)
        print(f"\n{st.name}, {args.experiment} experiment, {reading}: V_rev "
              f"{r.v * 1e3:+.2f} mV{'' if r.converged else ' (n.c.)'}")
        heads = list(s.constrictions)
        print(f"  {'ion':5s} {'bath lum/cyt M':>15s} {'peak M':>8s} {'at z':>7s} "
              f"{'I pA':>7s} {'steepest z':>11s}"
              + "".join(f"{h + ' share':>13s}" for h in heads))
        for ion in r.species:
            c, z = r.peak(ion)
            lum, cyt = r.baths[ion]
            print(f"  {ion:5s} {lum:7.3f}/{cyt:<7.3f} {c:8.3g} {z:+7.1f} "
                  f"{r.currents[ion] * 1e12:+7.2f} {r.steepest_z(ion):+11.1f}"
                  + "".join(f"{r.drop_across(ion, s.constrictions[h].z, w):13.0%}"
                            for h in heads))
        print("  " + ", ".join(f"{h} z {c.z:+.1f} A" for h, c in
                               s.constrictions.items())
              + f"; shares over +/- {w:g} A")


def register(sub) -> None:
    p = sub.add_parser("sel3d", help="P_Ca:P_K from the 3-D charged lumen, "
                       "with the charge-space excess (Round 7.19)")
    p.add_argument("pdb", nargs="*", help=f"deposits (default {' '.join(_DEFAULT)})")
    p.add_argument("--mutants", action="store_true",
                   help="Xu 2006's RyR1 mutants (give a RyR1 deposit)")
    p.add_argument("--spacing", type=float, default=None,
                   help="voxel spacing, A (default pore3d.spacing)")
    p.add_argument("--wall-volume", type=float, default=None,
                   help="csc.structural_volume for this run (0..1)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_sel3d)
    p = sub.add_parser("gate", help="the gate widened: conductance and "
                       "P_Ca:P_K again (Round 7.21)")
    p.add_argument("pdb", nargs="*", help="deposits (default 9HEO 8TKF)")
    p.add_argument("--delta", type=float, action="append", default=None,
                   help="widening, A (repeatable; default 0..gate.widen_max)")
    p.add_argument("--spacing", type=float, default=None,
                   help="voxel spacing, A (default pore3d.spacing)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_gate)
    p = sub.add_parser("reversal", help="P_Ca:P_K at bi-ionic reversal in "
                       "the 3-D charged lumen (Round 7.23)")
    p.add_argument("pdb", nargs="*", help=f"deposits (default {' '.join(_DEFAULT)})")
    p.add_argument("--mutants", action="store_true",
                   help="Xu 2006's RyR1 mutants (give a RyR1 deposit)")
    p.add_argument("--spacing", type=float, default=None,
                   help="voxel spacing, A (default reversal3d.spacing)")
    p.add_argument("--scale", type=float, action="append", default=None,
                   help="scale every wall charge (repeatable; pb only)")
    p.add_argument("--lumen", nargs="*", default=None, metavar="READING",
                   help="Round 7.24: each ion's peak and drop on the lumen at "
                   "one experiment's reversal, per reading (default pb + csc)")
    p.add_argument("--experiment", default="Ca2+", choices=("Ca2+", "Cl-"),
                   help="with --lumen: which experiment (Cl- is IP3R only)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_reversal)
