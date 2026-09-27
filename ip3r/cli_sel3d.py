"""Selectivity in 3-D (Round 7.19), ``register(sub)``:

    python -m ip3r sel3d                 # 9HEO, then 8TKF
    python -m ip3r sel3d 9HEO --mutants  # Xu 2006's five mutants on 9HEO
    python -m ip3r sel3d 7T3T --spacing 1.0
    python -m ip3r gate 9HEO             # Round 7.21: the gate widened
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
