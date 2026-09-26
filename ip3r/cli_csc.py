"""Charge-space competition (Round 7.17), ``register(sub)``:

    python -m ip3r csc                  # RyR1 9HEO vs Xu 2006, then 8TKF vs Vais 2010
    python -m ip3r csc 9HEO --scan      # wall volume, eps, water and oxygen size
"""

from __future__ import annotations

__all__ = ["register"]

_SCANS = {"csc.structural_volume": (0.0, 0.5, 1.0),
          "csc.permittivity": (40.0, 60.0, 78.4),
          "csc.water_diameter": (2.5, 2.8, 3.1),
          "csc.oxygen_diameter": (2.4, 2.8, 3.2)}

_DEFAULT = ("9HEO", "8TKF")


def _energetics(w) -> None:
    from .physics import csc_readings as cr
    print("  filter binding, Gillespie 2008's split (kT; + favours Ca2+):")
    for kcl, ca in ((0.15, 0.001), (0.25, 0.01)):
        e = cr.energetics(w, kcl, ca)
        print(f"    {kcl * 1e3:3.0f} mM KCl + {ca * 1e3:2.0f} mM Ca2+ at z "
              f"{e.z:6.1f} ({e.fixed:+.1f} M): K+ {e.k:5.2f} M, Ca2+ "
              f"{e.ca:5.2f} M, water {e.water:4.1f} M, packing {e.packing:.2f};"
              f" mean {e.mean:+.2f}, screening {e.screening:+.2f}, excluded "
              f"volume {e.excluded:+.2f}")


def _shares(w) -> None:
    from .physics import csc_readings as cr
    print("  linear-response resistance (symmetric mixed bath), filter / "
          "gate share, and the P_Ca:P_K it gives:")
    for closure in ("donnan", "csc"):
        s = cr.shares(w, closure)
        print(f"    {closure:6s} K+ {s.filter['K+']:4.0%} / {s.gate['K+']:4.0%}"
              f"   Ca2+ {s.filter['Ca2+']:4.0%} / {s.gate['Ca2+']:4.0%}"
              f"   P_Ca:P_K {s.ratio:5.2f}")
    print(f"    the gate window alone, uncharged, by hand: "
          f"{cr.gate_limit(w):.2f}")


def _csc(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import csc
    from .physics import csc_readings as cr
    loader.ALLOW_FETCH = args.fetch
    print("Charge-space competition: hard spheres (BMCSL) + screening (MSA, "
          f"eps {_P.value('csc.permittivity'):g}, l_B {csc.bjerrum():.2f} A), "
          f"water {_P.value('csc.water_diameter'):g} A at "
          f"{_P.value('csc.water_concentration'):g} M, carboxylate oxygens "
          f"{_P.value('csc.oxygen_diameter'):g} A.")
    for pdb in args.pdb or _DEFAULT:
        st = loader.load(pdb)
        w = cr.csc_wall(st)
        print(f"\n{pdb}: formal wall {w.net:+.0f} e")
        _energetics(w)
        _shares(w)
        if args.scan:
            if not w.ryr:
                print("  (--scan reads Xu's protocol: RyR1 deposits only)")
            else:
                for key, values in _SCANS.items():
                    print(f"  {key}: value -> wild type P_Ca:P_K, D4899Q x")
                    for v, p, x in cr.ryr1_scan(st, key, values):
                        print(f"    {v:5.2f}: {p:5.2f}  x{x:4.2f}")
            continue
        heads = tuple(cr.READINGS)
        if w.ryr:
            print("  Xu 2006 (250 mM KCl, 10 mM CaCl2 luminal, their Eq. 1); "
                  "x = over wild type")
            print("  " + " " * 21 + "".join(f"{h:>24s}" for h in heads))
            rows = cr.ryr1_panel(st, heads)
            wt = rows[0]
            for r in rows:
                cells = "".join(
                    f"{r.pca[h]:7.2f} x{r.pca[h] / wt.pca[h]:4.2f} "
                    f"g x{r.g[h] / wt.g[h]:4.2f}" for h in heads)
                flag = "" if r.converged else "  (n.c.)"
                print(f"  {r.name:9s} {r.measured_pca:4.1f} "
                      f"x{r.measured_pca / wt.measured_pca:4.2f} "
                      f"g x{r.measured_g / wt.measured_g:4.2f}"
                      f"{cells}{flag}")
            print("  wild-type g: " + ", ".join(
                f"{h} {wt.g[h]:.0f} pS" for h in heads)
                + f" (measured {wt.measured_g:.0f})")
        else:
            print("  Vais 2010 (measured P_Ca:P_K "
                  f"{_P.value('selectivity.published_pca_pk'):g}, P_Cl:P_K "
                  f"{_P.value('selectivity.published_pcl_pk'):g}):")
            for h, (pca, pcl, g, ok) in cr.vais_reading(st, heads).items():
                print(f"    {h:20s} P_Ca:P_K {pca:6.2f}   P_Cl:P_K {pcl:5.2f}"
                      f"   g {g:6.1f} pS{'' if ok else '  (n.c.)'}")
    return 0


def register(sub) -> None:
    p = sub.add_parser("csc", help="charge-space competition vs the measured "
                       "selectivities")
    p.add_argument("pdb", nargs="*", help=f"deposits (default {' '.join(_DEFAULT)})")
    p.add_argument("--scan", action="store_true",
                   help="wall volume, eps, water and oxygen diameters (RyR1)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_csc)
