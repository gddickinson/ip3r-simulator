"""The wall charge in 3-D (Round 7.11), ``register(sub)``:

    python -m ip3r wall3d [PDB ...]        # every open deposit, three closures
    python -m ip3r wall3d 8TKF --scan      # + width, lining margin, ε, grid
    python -m ip3r wall3d --mutants        # Xu 2006's RyR1 mutants in 3-D
"""

from __future__ import annotations

__all__ = ["register"]

_SCANS = {"pore_charge.smoothing": (1.5, 3.0, 4.5, 6.0),
          "pore_charge.lining_margin": (1.0, 3.0, 5.0, 8.0),
          "permeation.permittivity_pore": (20.0, 40.0, 60.0, 80.0)}


def _one_d(st, summary) -> tuple[float, float, float]:
    """The 1-D zero-voltage readings: neutral, charged, paired (pS)."""
    import numpy as np

    from .physics.charged3d import linear_response_1d
    from .physics.permeation import potassium_species
    from .physics.unitary import unitary
    u = unitary(st, summary)
    sp = potassium_species(bath=u.bath)
    r = np.maximum(u.profile.r_free, 0.0)
    return tuple(linear_response_1d(u.profile.z, r, sp, fx) * 1e12
                 for fx in (None, u.charge.density, u.paired_charge.density))


def _wall3d(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import charged3d as c3
    from .physics.shortfall import open_entries
    from .structure.channel import measure_channel
    loader.ALLOW_FETCH = args.fetch
    h = args.spacing or _P.value("pore3d.spacing")
    print("K+ conductance at zero voltage (linear response: the ions' "
          "equilibrium, then one weighted Laplace solve per species), neutral "
          "and with the lining charges placed three ways: slice = the 1-D "
          "charge per length over the real cross-section, local = each group "
          "a 3-D Gaussian at its own centre (both local Donnan), pb = local "
          f"screened by Poisson-Boltzmann (eps {_P.value('permeation.permittivity_pore'):g}). "
          f"h {h:g} A. Ratios are to the same deposit's neutral 3-D reading.")
    if args.mutants:
        wt, rows = c3.mutant_panel_3d(spacing=h, progress=lambda i, n, m: print(
            f"  {m} ({i + 1}/{n})", flush=True))
        print(wt.row())
        for r in rows:
            parts = "  ".join(f"{c} x{v:4.2f}" for c, v in r.ratios.items())
            print(f"{r.name:7s} measured x{r.measured_ratio:4.2f}   1-D "
                  f"x{r.one_d:4.2f}   3-D {parts}")
        return 0
    for pdb in args.pdb or [e.pdb_id for e in open_entries()]:
        st = loader.load(pdb)
        s = measure_channel(st)
        n1, c1, p1 = _one_d(st, s)
        print(f"\n{pdb}: 1-D (zero voltage) neutral {n1:6.1f}  charged {c1:6.1f} "
              f"(x{c1 / n1:4.2f})  paired {p1:6.1f} (x{p1 / n1:4.2f}) pS")
        for paired in (False, True):
            w = c3.wall_3d(st, s, spacing=h, pair_bridges=paired)
            peaks = ", ".join(f"{c} {v:.1f}" for c, v in w.peaks().items())
            print(f"  {'paired ' if paired else 'charged'} {w.row()}")
            print(f"          wall {w.charge.net_charge:+.0f} e; peak K+ {peaks} M")
        if args.scan:
            for key, values in _SCANS.items():
                print(f"  {key} (h {max(h, 1.0):g} A)")
                for v, w in c3.wall_scan(st, key, values, spacing=max(h, 1.0)):
                    print(f"    {v:5.1f}: {w.row()}")
            print("  grid")
            for g in (1.0, 0.75, 0.5):
                print(f"    h {g:4.2f}: {c3.wall_3d(st, s, spacing=g).row()}")
    return 0


def register(sub) -> None:
    p = sub.add_parser("wall3d", help="the wall charge in the lumen's 3-D shape")
    p.add_argument("pdb", nargs="*", help="deposits (default: every open one)")
    p.add_argument("--scan", action="store_true",
                   help="width, lining margin, permittivity and grid scans")
    p.add_argument("--mutants", action="store_true",
                   help="Xu 2006's RyR1 charge mutants in 3-D")
    p.add_argument("--spacing", type=float, default=None, help="voxel edge, A")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_wall3d)
