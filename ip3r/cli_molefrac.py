"""The mole-fraction prediction (Round 7.29), ``register(sub)``:

    python -m ip3r molefrac                   # 8TKF, the crossing site
    python -m ip3r molefrac 7T3T --depth 3.98 --only site
"""

from __future__ import annotations

__all__ = ["register"]


def _progress(i, n, label):
    print(f"  [{i + 1}/{n}] {label}", flush=True)


def _report(m, published) -> None:
    from .physics.mole_fraction import half_point
    print(f"\n{m.name}: band '{m.region}', hold {m.voltage * 1e3:+.0f} mV "
          f"(cytosol - lumen); the pore's P_Cl:P_K {m.pcl_pk:.3f}")
    lo, hi = published["ica_range"]
    for c in m.curves.values():
        head = c.reading if c.site is None else (
            f"{c.reading}: {c.site.label}, K_d {c.kd:.3g} mM")
        print(f"\n  {head}; i_K with no Ca2+ {c.i_k0 * 1e12:.3f} pA")
        print(f"    {'Ca mM':>7s} {'V_rev mV':>9s} {'P_Ca:P_K':>8s} "
              f"{'i_Ca(0) pA':>10s} {'i_K pA':>8s} {'i_K/0':>6s} "
              f"{'i_tot/0':>7s} {'theta':>6s} {'held':>5s}")
        kf, tf = c.k_fraction, c.total_fraction
        for p, k, t in zip(c.points, kf, tf):
            print(f"    {p.ca * 1e3:7.3g} {p.v_rev * 1e3:+9.2f} "
                  f"{p.pca_pk:8.2f} {p.i_ca * 1e12:10.4f} "
                  f"{p.i_hold['K+'] * 1e12:8.3f} {k:6.3f} {t:7.3f} "
                  f"{p.occupancy:6.3f} {p.held:5.2f}"
                  f"{'' if p.converged else '  (n.c.)'}", flush=True)
        half = half_point(c.ca, kf)
        slope = c.ica_slope(lo, hi) * 1e12 / 1e3            # pA/mM
        print(f"    i_K halved at {half * 1e3:.3g} mM" if half == half
              else "    i_K never halved on the grid")
        print(f"    i_Ca slope at 0 mV over {lo * 1e3:g}-{hi * 1e3:g} mM: "
              f"{slope:.4f} pA/mM (measured {published['ica_slope']:g})")
        pca = [p.pca_pk for p in c.points]
        print(f"    P_Ca:P_K {min(pca):.2f}-{max(pca):.2f} over the sweep "
              f"(measured {published['pca_pk']:g} at "
              f"{published['cacl2'] * 1e3:g} mM)")


def _molefrac(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics.mole_fraction import READINGS, sweep
    from .physics.selectivity import published as _pub
    loader.ALLOW_FETCH = args.fetch
    readings = tuple(args.only or READINGS)
    published = dict(_pub(), cacl2=_P.value("selectivity.cacl2_lumen"),
                     ica_range=(_P.value("selectivity.ca_lumen_low"),
                                _P.value("selectivity.ca_lumen_high")))
    print("The mole-fraction prediction (Round 7.29): Vais 2010's 140 mM KCl "
          "both sides, luminal CaCl2 swept; V_rev -> P_Ca:P_K (Eq. 1), i_Ca at "
          "0 mV, i_K at the hold against the same pore without Ca2+. Site = "
          "Round 7.27's compensated site blocking K+ "
          f"(f {_P.value('casite.block'):g}), {_P.value('casite.sites'):g} "
          "sites over the band. Uncharged wall, point ions.")
    for pdb in args.pdb or ("8TKF",):
        st = loader.load(pdb)
        try:
            m = sweep(st, readings, depth=args.depth, region=args.region,
                      spacing=args.spacing, progress=_progress)
        except ValueError as e:
            print(f"\n{pdb}: refused: {e}")
            continue
        _report(m, published)
    return 0


def register(sub) -> None:
    from .physics.mole_fraction import READINGS
    p = sub.add_parser("molefrac", help="P_Ca:P_K, i_Ca and the K+ current "
                       "against luminal Ca2+ through the Ca2+ site (Round 7.29)")
    p.add_argument("pdb", nargs="*", help="IP3R deposits (default 8TKF)")
    p.add_argument("--depth", type=float, default=None,
                   help="the site's depth, kT (default casite.gui_depth)")
    p.add_argument("--region", default="span",
                   choices=("span", "filter", "gate", "filter to gate",
                            "vestibule"), help="the band the sites cover")
    p.add_argument("--only", action="append", default=None, choices=READINGS,
                   help="readings to run (repeatable; default all three)")
    p.add_argument("--spacing", type=float, default=None,
                   help="voxel spacing, A (default reversal3d.spacing)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_molefrac)
