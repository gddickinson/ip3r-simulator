"""The conductance shortfall (Round 7.6), ``register(sub)``:

    python -m ip3r shortfall            # every open deposit, three geometries
    python -m ip3r shortfall --scan     # + window, seal and grid scans on 8TKF
"""

from __future__ import annotations

__all__ = ["register"]


def _shortfall(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import shortfall as sf
    from .structure.pore import PROFILE_MARGIN
    loader.ALLOW_FETCH = args.fetch
    print("Neutral K+ conductance of every open deposit, one electrolyte, three "
          "geometries: 1-D inscribed circle (Round 4), the voxelised lumen in "
          f"3-D (Laplace; h {_P.value('pore3d.spacing'):g} A and extrapolated to "
          "h 0 from twice it), and that lumen as a slice-area integral. "
          f"In-pore diffusivity {_P.value('permeation.diffusion_scale'):g} x bulk.")
    rows = sf.open_panel(corner=True, progress=lambda i, n, p: print(
        f"  measuring {p} ({i + 1}/{n})", flush=True))
    for r in rows:
        print(r.row())
    print()
    for r in rows:
        print(f"{r.pdb_id}: 3-D / 1-D {r.gain:.2f}x; short of "
              f"{min(r.measured.values()):.0f} pS by {r.short_by():.1f}x "
              f"(of {max(r.measured.values()):.0f} by "
              f"{r.short_by(max(r.measured.values())):.1f}x); diffusivity needed "
              f"{r.needed_scale():.2f}-{r.needed_scale(max(r.measured.values())):.2f}"
              f" x bulk; favourable corner (bulk D, K+ radius "
              f"{_P.value('permeation.sweep_radius_low'):g} A) {r.best_pS:.0f} pS "
              f"(h 0), {min(r.measured.values()) / r.best_pS:.2f}x short")
    if args.scan:
        pdb = args.pdb
        print(f"\n{pdb}: 1-D window (S0's margin {PROFILE_MARGIN:g} A scaled)")
        for margin, g, share in sf.window_scan(pdb):
            print(f"  margin {margin:5.1f} A: {g:6.1f} pS, access {share:.1%} of R")
        print(f"{pdb}: 3-D seal radius (coarse grid)")
        for seal, g in sf.seal_scan(pdb):
            print(f"  seal {seal:5.1f} A: {g:8.1f} pS")
        print(f"{pdb}: 3-D grid")
        for h, g in sf.grid_scan(pdb):
            print(f"  h {h:4.2f} A: {g:6.1f} pS")
    return 0


def register(sub) -> None:
    p = sub.add_parser("shortfall", help="the open pore in 1-D and 3-D against "
                       "the measured conductance")
    p.add_argument("pdb", nargs="?", default="8TKF",
                   help="the deposit --scan runs on")
    p.add_argument("--scan", action="store_true",
                   help="window, seal and grid scans")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_shortfall)
