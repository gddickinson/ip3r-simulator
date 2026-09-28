"""The saturable Ca²⁺ site and its K⁺ block (Round 7.27), ``register(sub)``:

    python -m ip3r casite                      # 8TKF, the site over the span
    python -m ip3r casite 8TKF 7T3T --region "filter to gate"
"""

from __future__ import annotations

__all__ = ["register"]

KINDS = {"uncompensated": (False, 0.0), "compensated": (True, 0.0),
         "uncompensated + block": (False, None),
         "compensated + block": (True, None)}


def _progress(i, n, label):
    print(f"  [{i + 1}/{n}] {label}", flush=True)


def _casite(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics.ca_site import search
    loader.ALLOW_FETCH = args.fetch
    kinds = tuple(KINDS[k] for k in (args.kind or KINDS))
    print("A saturable Ca2+ site at Vais 2010's reversal (Round 7.27): "
          f"{_P.value('casite.sites'):g} sites over the band, depth d = the "
          "empty site's pull on Ca2+; compensated = -2e fixed per bound "
          "Ca2+; block = K+ stopped by an occupied site with probability "
          f"{_P.value('casite.block'):g}. Uncharged wall, point ions.")
    for pdb in args.pdb or ("8TKF",):
        st = loader.load(pdb)
        try:
            s = search(st, spacing=args.spacing, region=args.region,
                       kinds=kinds, solve_required=not args.no_required,
                       progress=_progress)
        except ValueError as e:
            print(f"\n{pdb}: refused: {e}")
            continue
        n = s.neutral
        print(f"\n{pdb}: band '{args.region}', site density "
              f"{s.density / 1000:.3g} M; measured P_Cl:P_K {s.measured[0]:g}, "
              f"P_Ca:P_K {s.measured[1]:g}; uncharged pore {n.pcl_pk:.3f}, "
              f"{n.pca_pk:.3f} (score {n.score:.2f})")
        print(f"    {'site':44s} {'P_Cl:P_K':>8s} {'P_Ca:P_K':>8s} "
              f"{'B':>8s} {'score':>6s} {'V_Ca mV':>8s} {'theta':>6s} "
              f"{'held':>5s} {'peak M':>7s}")
        for r in s.results:
            print(f"    {r.site.label:44s} {r.pcl_pk:8.3f} {r.pca_pk:8.2f} "
                  f"{r.product:8.2f} {r.score:6.2f} {r.v['Ca2+'] * 1e3:+8.2f} "
                  f"{r.occupancy:6.3f} {r.bound:5.2f} {r.peak_ca:7.3g}"
                  f"{'' if r.converged else '  (n.c.)'}", flush=True)
        b = s.best
        print(f"    best: {b.site.label} (score {b.score:.2f})")
        for label, d in s.required.items():
            print(f"    P_Ca:P_K = {s.measured[1]:g} at d = {d:.2f} kT: {label}")
        if not args.no_required and not s.required:
            print(f"    no kind crosses P_Ca:P_K {s.measured[1]:g} on the grid")
    return 0


def register(sub) -> None:
    p = sub.add_parser("casite", help="a saturable Ca2+ site, compensated or "
                       "not, and its K+ block, at Vais's reversal (Round 7.27)")
    p.add_argument("pdb", nargs="*", help="IP3R deposits (default 8TKF)")
    p.add_argument("--region", default="span",
                   choices=("span", "filter", "gate", "filter to gate"),
                   help="the band the sites cover")
    p.add_argument("--kind", action="append", default=None,
                   choices=tuple(KINDS), help="kinds to run (repeatable; "
                   "default all four)")
    p.add_argument("--no-required", action="store_true",
                   help="skip the root-finding for the depth giving 15.2")
    p.add_argument("--spacing", type=float, default=None,
                   help="voxel spacing, A (default reversal3d.spacing)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_casite)
