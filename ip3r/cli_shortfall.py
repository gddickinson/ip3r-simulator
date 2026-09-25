"""The conductance shortfall (Round 7.6), ``register(sub)``:

    python -m ip3r shortfall            # every open deposit, three geometries
    python -m ip3r shortfall --scan     # + window, seal and grid scans on 8TKF
    python -m ip3r lumen [PDB ...]      # where the voltage falls (Round 7.10)
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


def _lumen(args) -> int:
    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics.lumen_field import lumen_field
    from .physics.shortfall import open_entries
    from .structure.channel import measure_channel
    loader.ALLOW_FETCH = args.fetch
    w = _P.value("lumen.constriction_half_width")
    print("Where the applied voltage falls in the neutral pore: the 3-D Laplace "
          "potential (Round 7.6's lumen) against the 1-D model's inscribed "
          "circle, each normalised to its drop across S0's window; share of "
          f"that drop within ±{w:g} A of each constriction.")
    for pdb in args.pdb or [e.pdb_id for e in open_entries()]:
        st = loader.load(pdb)
        s = measure_channel(st)
        f = lumen_field(st, s)
        print(f.summary())
        for c in s.constrictions.values():
            d3, d1 = f.drop_across(c.z, w)
            print(f"  {c.name:6s} z {c.z:+7.1f} A, r_min {c.radius:.2f} A: "
                  f"3-D {d3:.1%}, 1-D {d1:.1%}")
    return 0


def register(sub) -> None:
    q = sub.add_parser("lumen", help="where the voltage falls: 3-D potential "
                       "against the 1-D model")
    q.add_argument("pdb", nargs="*", help="deposits (default: every open one)")
    q.add_argument("--fetch", action="store_true")
    q.set_defaults(fn=_lumen)
    p = sub.add_parser("shortfall", help="the open pore in 1-D and 3-D against "
                       "the measured conductance")
    p.add_argument("pdb", nargs="?", default="8TKF",
                   help="the deposit --scan runs on")
    p.add_argument("--scan", action="store_true",
                   help="window, seal and grid scans")
    p.add_argument("--fetch", action="store_true")
    p.set_defaults(fn=_shortfall)
