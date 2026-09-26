"""The lining salt bridge (Round 7.13), ``register(sub)``:

    python -m ip3r bridge [PDB ...]      # pKas of the pair, then its field
    python -m ip3r bridge 8TKF --scan    # + ε protein, width, ε water, box
    python -m ip3r bridge --mutants      # Xu 2006's RyR1 mutants, dielectric
"""

from __future__ import annotations

__all__ = ["register"]

_SCANS = {"dielectric.eps_protein": (2.0, 4.0, 10.0, 20.0),
          "dielectric.charge_width": (0.5, 1.0, 2.0, 3.0),
          "permeation.permittivity_pore": (20.0, 40.0, 80.0),
          "pore3d.box_half_width": (30.0, 40.0),
          "pore3d.bath_margin": (25.0, 35.0)}

_DEFAULT = ("8TKF", "7T3T", "9HEO")


def _bridge(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import bridge_charge as bc
    from .physics.protonation import lining_wall
    loader.ALLOW_FETCH = args.fetch
    h = args.spacing
    if args.mutants:
        from .physics.charged3d import mutant_panel_3d
        wt, rows = mutant_panel_3d(closures=("pb", "dielectric"), spacing=h)
        print(f"RyR1 9HEO, h {h:g} A, eps protein "
              f"{_P.value('dielectric.eps_protein'):g}: {wt.row()}")
        for r in rows:
            parts = "  ".join(f"{c} x{v:4.2f}" for c, v in r.ratios.items())
            print(f"  {r.name:7s} measured x{r.measured_ratio:4.2f}   1-D "
                  f"x{r.one_d:4.2f}   3-D {parts}")
        return 0
    print("Each lining salt bridge: its pKas (network with and without the "
          "base among the sites; PROPKA), then K+ g/g_neutral in 3-D with the "
          "protein in the field (dielectric closure, eps protein "
          f"{_P.value('dielectric.eps_protein'):g}, water "
          f"{_P.value('permeation.permittivity_pore'):g}; h {h:g} A).")
    for pdb in args.pdb or _DEFAULT:
        st = loader.load(pdb)
        wall = lining_wall(st)
        for acid, base in bc.lining_bridges(st, wall):
            print(f"\n{pdb}: D{acid}-R{base}' (pH {_ph(wall):g})")
            for t in bc.titrate_pair(st, acid, base, wall,
                                     with_propka=not args.no_propka):
                print(f"  {t.row()}")
            for label, w in bc.pair_readings(st, acid, base, wall.summary,
                                             spacing=h).items():
                c = bc.READINGS[label][0]
                f = w.fields[c]
                print(f"  {label:13s} g x{w.ratio(c):5.2f}  (charge in box "
                      f"{f.placed:+5.1f} e){'' if w.converged else '  [n.c.]'}")
            if args.scan:
                for key, values in _SCANS.items():
                    print(f"  {key}")
                    for v, r in bc.pair_scan(st, acid, base, key, values, h):
                        print(f"    {v:5.1f}: " + "  ".join(
                            f"{k} x{x:4.2f}" for k, x in r.items()))
    return 0


def _ph(wall) -> float:
    from .physics.protonation import family_conditions
    return family_conditions(wall.ryr)[0]


def register(sub) -> None:
    p = sub.add_parser("bridge", help="the lining salt bridge: pKas and field")
    p.add_argument("pdb", nargs="*", help=f"deposits (default {' '.join(_DEFAULT)})")
    p.add_argument("--scan", action="store_true",
                   help="eps protein, charge width, eps water and box scans")
    p.add_argument("--mutants", action="store_true",
                   help="Xu 2006's RyR1 mutants under the dielectric closure")
    p.add_argument("--no-propka", action="store_true")
    p.add_argument("--spacing", type=float, default=1.0,
                   help="voxel edge, A (default 1.0: the whole box is solved)")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_bridge)
