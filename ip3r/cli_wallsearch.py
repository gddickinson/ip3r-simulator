"""The wall-model search (Round 7.25), ``register(sub)``:

    python -m ip3r wallsearch                  # 8TKF: ceilings, rings, reversal
    python -m ip3r wallsearch 8TKF --only well
    python -m ip3r wallsearch 7T3T --no-reversal
"""

from __future__ import annotations

__all__ = ["register"]


def _progress(i, n, label):
    print(f"  [{i + 1}/{n}] {label}", flush=True)


def _ceilings(st, spacing) -> None:
    from .parameters import PARAMETERS as _P
    from .physics.selectivity3d import prepare
    from .physics.selectivity_bound import well_ceiling
    from .physics.wall_search import regions
    pore = prepare(st, spacing=spacing or _P.value("reversal3d.spacing"))
    print("\n  linear response: the most a Ca2+-only well of any depth can "
          "raise P_Ca:P_K (x the uncharged pore; 1/x = Ca2+'s resistance "
          "left outside it)")
    for name, (lo, hi) in regions(pore.summary).items():
        print(f"    {name:15s} z {lo:+7.1f} .. {hi:+7.1f} A   "
              f"x{well_ceiling(pore, lo, hi):8.1f}", flush=True)


def _rings(st, spacing) -> None:
    from .physics.wall_search import ring_search
    rows = ring_search(st, spacing=spacing)
    print("\n  opposite-charge C4 rings (linear response, x the uncharged "
          "pore): the parallel-path route round B <= 1")
    print(f"    {'z A':>7s} {'e/site':>6s} {'P_Ca:P_K':>9s} {'P_Cl:P_K':>9s} "
          f"{'B':>7s}")
    for r in rows:
        print(f"    {r.z:+7.1f} {r.charge:6g} {r.pca_pk:9.3f} {r.pcl_pk:9.3f} "
              f"{r.bound:7.3f}{'' if r.converged else '  (n.c.)'}")
    top = max(rows, key=lambda r: r.bound)
    print(f"    highest B {top.bound:.3f} at z {top.z:+.1f} A, {top.charge:g} e")


def _search(st, families, spacing) -> None:
    from .physics.wall_search import search
    s = search(st, families, spacing=spacing, progress=_progress)
    m_cl, m_ca = s.measured
    n = s.neutral
    b_meas = (m_cl / n.pcl_pk) * (m_ca / n.pca_pk) ** 2
    print(f"\n  {s.name} at Vais 2010's reversal (point ions); measured "
          f"P_Cl:P_K {m_cl:g}, P_Ca:P_K {m_ca:g} (B {b_meas:.0f} over the "
          f"uncharged pore's {n.pcl_pk:.3f}, {n.pca_pk:.3f})")
    print(f"    {'candidate':38s} {'P_Cl:P_K':>9s} {'P_Ca:P_K':>9s} "
          f"{'B':>7s} {'score':>6s} {'V_Ca mV':>8s} {'peak Ca2+ M':>11s}")
    for r in [n] + s.results:
        print(f"    {r.candidate.label:38s} {r.pcl_pk:9.3f} {r.pca_pk:9.2f} "
              f"{r.bound:7.3f} {r.score:6.2f} {r.v['Ca2+'] * 1e3:+8.2f} "
              f"{r.peak_ca:11.3g}{'' if r.converged else '  (n.c.)'}")
    b = s.best
    print(f"    best: {b.candidate.label} (score {b.score:.2f}; the uncharged "
          f"pore scores {n.score:.2f})")


def _wallsearch(args) -> int:
    from .io import loader
    from .physics.wall_search import FAMILIES
    loader.ALLOW_FETCH = args.fetch
    print("Wall models scored on P_Cl:P_K and P_Ca:P_K at once (Round 7.25). "
          "B = (P_Cl:P_K)(P_Ca:P_K)^2, each over the uncharged pore's: at "
          "most 1 for any point-ion potential in series.")
    for pdb in args.pdb or ("8TKF",):
        st = loader.load(pdb)
        print(f"\n{pdb}")
        if not args.no_ceiling:
            _ceilings(st, args.spacing)
        if not args.no_rings:
            _rings(st, args.spacing)
        if not args.no_reversal:
            try:
                _search(st, tuple(args.only or FAMILIES), args.spacing)
            except ValueError as e:
                print(f"\n  search at reversal refused: {e}")
    return 0


def register(sub) -> None:
    p = sub.add_parser("wallsearch", help="wall models scored on P_Cl:P_K and "
                       "P_Ca:P_K at once (Round 7.25)")
    p.add_argument("pdb", nargs="*", help="IP3R deposits (default 8TKF)")
    p.add_argument("--only", action="append", default=None,
                   choices=("charge", "well", "well + charge"),
                   help="reversal families to run (repeatable; default all)")
    p.add_argument("--no-ceiling", action="store_true",
                   help="skip the linear-response well ceilings")
    p.add_argument("--no-rings", action="store_true",
                   help="skip the opposite-charge ring pairs")
    p.add_argument("--no-reversal", action="store_true",
                   help="skip the search at reversal (the slow part)")
    p.add_argument("--spacing", type=float, default=None,
                   help="voxel spacing, A (default reversal3d.spacing)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_wallsearch)
