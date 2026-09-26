"""The image (Born) cost of the low-ε wall (Round 7.15), ``register(sub)``:

    python -m ip3r born [PDB ...]        # W on the axis, then K+ g with it
    python -m ip3r born 8TKF --scan      # box, reach, eps protein
    python -m ip3r born --mutants        # Xu 2006's RyR1 mutants with it
"""

from __future__ import annotations

__all__ = ["register"]

_SCANS = {"born.box_half_width": (8.0, 12.0, 16.0),
          "born.reach": (6.0, 10.0, 14.0),
          "dielectric.eps_protein": (2.0, 4.0, 10.0)}

_DEFAULT = ("8TKF", "7T3T", "9HEO")


def _axis(r) -> str:
    return "  ".join(f"{k} {w:4.2f} kT ({z:6.1f} A)" for k, (z, w) in r.axis.items())


def _born(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import born_readings as br
    loader.ALLOW_FETCH = args.fetch
    h = args.spacing
    print("Image (Born) self-energy W of a unit charge, kT (x z^2 for an ion): "
          f"protein eps {_P.value('dielectric.eps_protein'):g}, water "
          f"{_P.value('permeation.permittivity_pore'):g}, box +-"
          f"{_P.value('born.box_half_width'):g} A, reach "
          f"{_P.value('born.reach'):g} A, h {h:g} A. K+ g with the dielectric "
          "closure, x = over the neutral pore without the image.")
    if args.mutants:
        for label, (wt, rows) in br.born_mutants(h, args.workers).items():
            print(f"\n{label}: wild type {wt.charged['dielectric'] * 1e12:.1f} pS")
            for r in rows:
                print(f"  {r.name:7s} measured x{r.measured_ratio:4.2f}   "
                      f"model x{r.ratios['dielectric']:4.2f}")
        return 0
    for pdb in args.pdb or _DEFAULT:
        st = loader.load(pdb)
        r = br.born_reading(st, h, args.workers, cache=not args.no_cache)
        b = r.born
        where = "cached" if b.cached else f"{b.seconds:.0f} s"
        meas = " / ".join(f"{v:.0f}" for v in r.measured.values())
        print(f"\n{pdb}: {b.solved} voxels solved ({where}); measured {meas} pS")
        print(f"  axis W: {_axis(r)}")
        if r.bridge:
            print(f"  lining pair D{r.bridge[0]}-R{r.bridge[1]}'")
        for line in r.rows():
            print(f"  {line}")
        if args.scan:
            for key, values in _SCANS.items():
                print(f"  {key}")
                for v, s in br.born_scan(st, key, values, h, args.workers):
                    print(f"    {v:5.1f}: filter W {s.axis['filter'][1]:4.2f}  "
                          f"neutral + image x{s.neutral_image / s.neutral:4.2f}  "
                          f"dipole, swept x{s.ratio('dipole, swept'):4.2f}  "
                          f"dipole + image x{s.ratio('dipole + image'):4.2f}")
    return 0


def register(sub) -> None:
    p = sub.add_parser("born", help="the image cost of the low-eps wall")
    p.add_argument("pdb", nargs="*", help=f"deposits (default {' '.join(_DEFAULT)})")
    p.add_argument("--scan", action="store_true",
                   help="box half-width, reach and protein eps scans")
    p.add_argument("--mutants", action="store_true",
                   help="Xu 2006's RyR1 mutants with the image cost")
    p.add_argument("--spacing", type=float, default=1.0,
                   help="voxel edge, A (default 1.0)")
    p.add_argument("--workers", type=int, default=None,
                   help="processes for the self-energy (default: every core)")
    p.add_argument("--no-cache", action="store_true",
                   help="solve W again even if data/cache/born has it")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_born)
