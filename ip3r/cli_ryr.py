"""The RyR1 commands of ``python -m ip3r`` (kept apart from ``cli.py``):

    python -m ip3r mutants [PDB]    # charge mutants: model vs Xu 2006
    python -m ip3r ryr-gating       # Stern 1997 scheme vs Murayama 2015 bell
    python -m ip3r sparks [--scan] [--cleft [--fit]]  # a RyR1 cluster,
                                              # read with the puff ruler
    python -m ip3r spark-termination [--scan fit|ki|rate]  # what ends a
                                              # cleft spark
"""

from __future__ import annotations

import json

__all__ = ["register"]


def _mutants(args) -> int:
    from .io import loader
    from .physics.ryr_mutants import mutant_panel, open_deposit
    loader.ALLOW_FETCH = args.fetch
    wt, rows = mutant_panel(loader.load(args.pdb or open_deposit()))
    print("RyR1 charge mutants, conductance / wild type (the model neutralises "
          "the residue on all four subunits)")
    print(wt.row())
    for r in rows:
        print(r.row())
    return 0


def _ryr_gating(args) -> int:
    from .physics import ryr_gating as rg
    for name, b in rg.compare_bells().items():
        print(f"{name:20s} peak {b.po_peak:.4f} at {b.c_peak:.1f} µM; "
              f"half-activation {b.c_half_act:.2f} µM, half-inhibition "
              f"{b.c_half_inh:.1f} µM ({b.width_decades:.2f} decades)")
    return 0


def _sparks(args) -> int:
    import numpy as np
    from .parameters import PARAMETERS as _P
    from .physics import puff_compare as pc
    model = (pc.SPARK_FIT if args.fit else pc.SPARK_CLEFT) if args.cleft else pc.SPARK
    if args.fit and not args.cleft:
        print("--fit applies to the cleft array; add --cleft")
        return 2
    if args.cleft:
        from .physics.cleft import coupling_matrix, nearest_coupling
        g = coupling_matrix()
        print(f"cleft coupling (Stern 1997 geometry): own release "
              f"{np.median(np.diag(g)):.1f} µM, nearest neighbour "
              f"{nearest_coupling(g):.2f} µM, all others open "
              f"{np.median(g.sum(axis=1) - np.diag(g)):.1f} µM (medians)")
    else:
        from .physics.sparks import diffusion_coupling
        print(f"derived coupling {diffusion_coupling():.2f} µM per open channel "
              "(free diffusion at the channel spacing)")
    if args.scan:
        print(f"{'µM/open':>8s} {'Fano':>5s} {'open':>6s} {'blips':>6s} "
              f"{'multi':>6s} {'large':>6s} {'/s':>5s}")
        for r in pc.spark_scan(args.duration, args.seed, model=model):
            print(f"{r['coupling']:8.3f} {r['fano']:5.2f} {r['open_fraction']:6.3f} "
                  f"{r['blips']:6d} {r['multi']:6d} {r['large']:6d} "
                  f"{r['large_per_s']:5.2f}")
        return 0
    tr = pc.simulate(model, 0.0, args.duration, args.seed,
                     pc.params_for(model, args.coupling))
    print(json.dumps(pc.recruitment(tr), indent=1))
    ends = pc.spark_ends(tr)
    if ends:
        med = {k: np.median([e[k] for e in ends])
               for k in ("duration", "inactivated_start", "inactivated_end")}
        print(f"spark duration median {1e3 * med['duration']:.0f} ms over "
              f"{len(ends)} sparks; measured release (frog) "
              f"{_P.value('spark.published_release_duration'):g} ms")
        unended = sum(e["unterminated"] for e in ends)
        if unended:
            print(f"{unended} spark(s) still running when the trace ended: "
                  "the median is a lower bound")
        print(f"channels inactivated: {med['inactivated_start']:.0f} at the "
              f"start, {med['inactivated_end']:.0f} at the end (medians of "
              f"{int(round(tr.params.n_channels))})")
    return 0


def _spark_termination(args) -> int:
    from .physics import ryr_gating as rg
    from .physics import spark_termination as st
    if args.scan == "fit":
        for name, b in rg.compare_bells().items():
            print(f"{name:20s} half-activation {b.c_half_act:6.2f} µM, "
                  f"half-inhibition {b.c_half_inh:6.1f} µM")
        rows = st.refit(args.duration, args.seeds)
    elif args.scan == "ki":
        rows = st.ki_scan(args.duration, args.seeds)
    else:
        base = rg.fit_to_bell() if args.fitted else None
        rows = st.rate_scan(base, args.duration, args.seeds)
    print(f"cleft array, {args.seeds} seeds x {args.duration:g} s; 'unended' = "
          "still running when the trace ends (its duration a lower bound)")
    for r in rows:
        print(r.row())
    return 0


def register(sub) -> None:
    p = sub.add_parser("mutants", help="RyR1 charge-neutralising mutants: "
                       "modelled vs measured conductance ratio")
    p.add_argument("pdb", nargs="?", default=None,
                   help="RyR1 deposit (default: the curated open state)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_mutants)
    p = sub.add_parser("ryr-gating", help="RyR1 bells: Stern 1997 scheme vs "
                       "Murayama 2015 fit, on one ruler")
    p.set_defaults(fn=_ryr_gating)
    p = sub.add_parser("sparks", help="a stochastic RyR1 cluster (Ca2+ sparks)")
    p.add_argument("--coupling", type=float, default=None,
                   help="µM per open channel; with --cleft, between nearest "
                   "neighbours (default: the derived coupling)")
    p.add_argument("--scan", action="store_true",
                   help="over the registered spark coupling band")
    p.add_argument("--cleft", action="store_true",
                   help="each channel sees its own Ca2+ in Stern's cleft "
                   "(default: one mean-field cluster Ca2+)")
    p.add_argument("--fit", action="store_true",
                   help="with --cleft: Ka and Ki fitted to Murayama's bell")
    p.add_argument("--duration", type=float, default=10.0)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(fn=_sparks)
    p = sub.add_parser("spark-termination", help="cleft spark duration as the "
                       "inactivation gate is refitted and scanned")
    p.add_argument("--scan", choices=("fit", "ki", "rate"), default="fit",
                   help="fit: Stern vs fitted to Murayama (25, 37 C); ki: Ki "
                   "scan; rate: inactivation rate scan at fixed Ki")
    p.add_argument("--fitted", action="store_true",
                   help="with --scan rate: scan at the fitted Ka and Ki")
    p.add_argument("--duration", type=float, default=10.0)
    p.add_argument("--seeds", type=int, default=4)
    p.set_defaults(fn=_spark_termination)
