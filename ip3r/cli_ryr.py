"""The RyR1 commands of ``python -m ip3r`` (kept apart from ``cli.py``):

    python -m ip3r mutants [PDB]    # charge mutants: model vs Xu 2006
    python -m ip3r ryr-gating       # Stern 1997 scheme vs Murayama 2015 bell
    python -m ip3r sparks [--scan] [--cleft [--fit]]  # a RyR1 cluster,
                                              # read with the puff ruler
    python -m ip3r spark-termination [--scan fit|ki|rate]  # what ends a
                                              # cleft spark
    python -m ip3r spark-mg [--scan] [--reading R] [--spontaneous] [--two-site]  # Mg2+ and
                                              # the triggered cleft spark
    python -m ip3r ec [--scan] [--reading R] [--trials N] [--two-site]  # the couplon under
                                              # voltage clamp: V channels
                                              # (Rios 1993) trigger C channels
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
    elif args.scan == "use":
        from .physics import ryr_use as ru
        b = rg.murayama_bell()
        print(f"Murayama 25 C: half-peak {b.c_half_act:.2f} - "
              f"{b.c_half_inh:.1f} µM, width {b.width_decades:.2f} decades. "
              "Below, the Ca2+ gate is refitted WITH the use gate present, "
              "so the bell is not counted twice. Rows scan the recovery "
              "ratio k_use-/k_use; the speed does not enter")
        _print_residual_bound()
        for r in ru.bell_panel():
            print(r.row())
        if args.bell:
            return 0
        rows = ([st.measure(rg.SternParams(), "Stern 1997", args.duration,
                            args.seeds),
                 st.measure(rg.fit_to_bell(), "fitted, Ca2+ gate only",
                            args.duration, args.seeds)]
                + st.use_scan(args.duration, args.seeds)
                + [st.no_inactivation(args.duration, args.seeds)])
    elif args.scan == "ratio":
        print("the use gate's recovery ratio k_use-/k_use at the registered "
              "0 mV rate. The Ca2+ gate is refitted at every point, so Ki "
              "moves with it: a larger ratio leaves the use gate less of the "
              "bell's descending limb and puts Ki back where Round 6.8 had it")
        _print_residual_bound()
        rows = st.ratio_scan(None, args.duration, args.seeds)
    elif args.scan == "fraction":
        print("the share of channels carrying the use gate (Laver & Lamb "
              "1998: 80 % of skeletal RyRs), at the registered ratio. The "
              "shared Ca2+ gate is refitted to the POPULATION bell at every "
              "point. The last two rows take the registered fraction apart")
        rows = (st.fraction_scan(None, args.duration, args.seeds)
                + st.fraction_controls(None, None, args.duration, args.seeds))
    elif args.scan == "low-activity":
        print("Copello 1997's low-activity channels (share ryr.la_fraction, "
              "Po <= 0.1) added to the population bell the high-activity "
              "channels' Ca2+ gate is fitted to, at each reading of their "
              "half points; they are not simulated in the cleft (best case)")
        rows = st.low_activity_scan(None, args.duration, args.seeds)
    else:
        base = rg.fit_to_bell() if args.fitted else None
        rows = st.rate_scan(base, args.duration, args.seeds)
    print(f"cleft array, {args.seeds} seeds x {args.duration:g} s; 'unended' = "
          "still running when the trace ends (its duration a lower bound)")
    for r in rows:
        print(r.row())
    return 0


def _print_residual_bound() -> None:
    from .parameters import PARAMETERS as _P
    from .physics.ryr_use import residual_bound
    lo, hi = (residual_bound(_P.value(f"ryr.use_residual_40mv_{k}"))
              for k in ("min", "max"))
    print(f"Laver & Lamb 1998 Fig. 8 bounds the ratio at +40 mV only: at most "
          f"{lo:.2g}-{hi:.2g} across the skeletal RyRs that inactivated; "
          "at 0 mV it is unmeasured")


def _spark_mg(args) -> int:
    from .physics import ryr_gating as rg
    from .physics import spark_mg as sm
    if args.two_site:
        from .physics.ryr_two_site import fit_two_site
        fit = fit_two_site()
        print(f"two-site inactivation: K1 {fit.k_i:.0f}, K2 {fit.k_i2:.0f} µM")
    else:
        fit = rg.fit_to_bell()
    reading = "selectivity" if args.ratio else args.reading
    k = sm.k_mg_a_reading(fit, reading)
    base = rg.with_mg(fit, 0.0, k)
    print(f"fitted to Murayama 25 C (Ka {fit.k_a:.2f}, Ki {fit.k_i:.0f} µM); "
          f"activation-site Mg2+ affinity {base.k_mg_a:.0f} µM "
          f"({sm.READINGS[reading]})")
    if args.spontaneous:
        rows = sm.spontaneous(base, k, args.duration, args.seeds)
    elif args.scan:
        rows = sm.mg_scan(base, k)
    else:
        rows = sm.dissect(base)
    if not args.spontaneous:
        print("triggered: every available channel opened at t = 0 and timed "
              "until none is open")
    for r in rows:
        print(r.row())
    return 0


def _ec(args) -> int:
    import numpy as np
    from .physics import ec_release as er
    from .physics.allosteric_v import open_probability
    if args.use and args.two_site:
        print("--use and --two-site are different schemes for the same gate; "
              "run them separately")
        return 1
    configs = er.configurations(args.reading, two_site=args.two_site,
                                use=args.use)
    if args.only:
        configs = {k: v for k, v in configs.items() if args.only.lower() in k.lower()}
    vs = er.voltages() if args.scan else np.array([0.0, -30.0, -50.0])
    if args.depletion:
        return _ec_depletion(args, configs)
    print("couplon: 30 V channels (Rios 1993 fiber 827, Stern's rates) and 30 "
          f"C channels; V steady Po {', '.join(f'{v:.0f} mV {open_probability(v):.4f}' for v in vs)}")
    print(f"activation-site Mg2+ reading: {args.reading}; 'after' = C open "
          "probability over the second half of the time after repolarisation")
    for label, sp in configs.items():
        for v in vs:
            e = er.ensemble(v, er.with_gating(sp), args.trials)
            print(er.summarise(e, label).row(), flush=True)
        print(er.event_stats(sp, label, trials=args.trials).row(), flush=True)
    return 0


def _ec_depletion(args, configs) -> int:
    from .parameters import PARAMETERS as P
    from .physics import lumen
    lp = lumen.LumenParams()
    n = 30 if args.large else int(P.value("lumen.fig20_channels"))
    print(f"SR content {lp.content} mM, {lp.density} couplons/um^3, refill tau "
          f"{lp.refill_tau} s; couplon {2 * n} channels ({n} C). Measured: a "
          f"100-ms pulse to +{P.value('lumen.rios_voltage'):.0f} mV releases "
          f"{P.value('lumen.rios_released_low'):.0%}-"
          f"{P.value('lumen.rios_released_high'):.0%} (Rios 1993); Stern's "
          "Fig. 20 (28 channels) leaves 0.66 mM at 0 mV and 1.34 at -30 mV")
    vs = (P.value("lumen.rios_voltage"), 0.0, -30.0)
    if args.pool_scan:
        for f, d in lumen.pool_scan(configs, vs[0], (1.0, 1.5, 2.0, 3.0),
                                    args.trials, small=not args.large):
            print(f"pool x{f:<4}", d.row(), flush=True)
        return 0
    for d in lumen.depletion_panel(configs, vs, args.trials, small=not args.large):
        print(d.row(), flush=True)
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
    p.add_argument("--scan", choices=("fit", "ki", "rate", "use", "ratio",
                                        "fraction", "low-activity"),
                   default="fit",
                   help="fit: Stern vs fitted to Murayama (25, 37 C); ki: Ki "
                   "scan; rate: inactivation rate scan at fixed Ki; use: the "
                   "use gate's speed; ratio: its recovery/inactivation ratio; "
                   "fraction: the share of channels that carry it; "
                   "low-activity: Copello 1997's LA channels in the bell")
    p.add_argument("--bell", action="store_true",
                   help="with --scan use: recovery ratios against the measured "
                        "bell only "
                        "(no simulation): where a Ca2+ gate still fits")
    p.add_argument("--fitted", action="store_true",
                   help="with --scan rate: scan at the fitted Ka and Ki")
    p.add_argument("--duration", type=float, default=10.0)
    p.add_argument("--seeds", type=int, default=4)
    p.set_defaults(fn=_spark_termination)
    p = sub.add_parser("spark-mg", help="cytosolic Mg2+ and the cleft spark: "
                       "the two sites dissected, or scanned")
    p.add_argument("--scan", action="store_true",
                   help="triggered sparks over free Mg2+ (default: the two "
                   "sites dissected at the fibre's Mg2+)")
    p.add_argument("--ratio", action="store_true",
                   help="activation-site Mg2+ affinity from Laver 2004's "
                   "Mg2+/Ca2+ selectivity instead of their absolute value")
    p.add_argument("--reading", choices=("measured", "selectivity", "meissner"),
                   default="measured", help="which activation-site Mg2+ "
                   "affinity (--ratio is --reading selectivity)")
    p.add_argument("--spontaneous", action="store_true",
                   help="untriggered runs over free Mg2+")
    p.add_argument("--two-site", action="store_true",
                   help="the two-site inactivation fit (slope matched too)")
    p.add_argument("--duration", type=float, default=10.0)
    p.add_argument("--seeds", type=int, default=4)
    p.set_defaults(fn=_spark_mg)
    p = sub.add_parser("ec", help="E-C coupling: the couplon under voltage "
                       "clamp, V channels triggering C channels")
    p.add_argument("--scan", action="store_true",
                   help="over the registered voltage scan (default 0, -30, -50 mV)")
    p.add_argument("--reading", choices=("measured", "selectivity", "meissner"),
                   default="meissner", help="activation-site Mg2+ affinity")
    p.add_argument("--only", default=None,
                   help="only the configurations whose label contains this")
    p.add_argument("--trials", type=int, default=None,
                   help="couplons per ensemble (default ec.trials)")
    p.add_argument("--two-site", action="store_true",
                   help="the C scheme fitted with a two-site inactivation gate "
                        "(the bell's inhibitory slope matched too)")
    p.add_argument("--use", action="store_true",
                   help="the use-dependent (flux-driven) inactivation gate, "
                        "with the Ca2+ gate refitted beside it at each of the "
                        "speeds the measured bell permits")
    p.add_argument("--depletion", action="store_true",
                   help="with the SR emptying (Stern's Fig. 20 pool)")
    p.add_argument("--large", action="store_true",
                   help="with --depletion: the 60-channel couplon instead of "
                   "Stern's 28-channel Fig. 20 one")
    p.add_argument("--pool-scan", action="store_true",
                   help="with --depletion: the pool x1-x3 at Rios's +20 mV")
    p.set_defaults(fn=_ec)
