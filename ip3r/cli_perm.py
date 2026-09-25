"""Permeation commands beyond the K+ conductance, ``register(sub)``:

    python -m ip3r selectivity [8TKF]   # Vais 2010's P_Cl:P_K, P_Ca:P_K, i_Ca
    python -m ip3r selectivity --closure radial   # PB across each slice
    python -m ip3r protonation [8TKF|9HEO] [--corners] [--no-propka]
"""

from __future__ import annotations

__all__ = ["register"]


def _selectivity(args) -> int:
    import numpy as np

    from .io import loader
    from .parameters import PARAMETERS as _P
    from .physics import selectivity as sel
    from .physics.unitary import permeation_profile
    from .structure.channel import measure_channel
    loader.ALLOW_FETCH = args.fetch
    st = loader.load(args.pdb)
    summary = measure_channel(st)
    rows = sel.selectivity_panel(st, summary, closure=args.closure)
    pub = sel.published()
    print(f"{st.name}: Vais 2010's lum-out protocols through the drift-diffusion "
          f"pore, {args.closure} closure (V_rev: KCl gradient / luminal CaCl2)")
    for r in rows:
        print(r.row())
    prof = permeation_profile(st, summary)
    bound = sel.slow_nmdg_bound(prof.z, np.maximum(prof.r_free, 0.0))
    print(f"neutral P_Cl:P_K with NMDG+ inside the pore (slow) instead of "
          f"excluded at the mouth: {bound:.2f}")
    print(f"measured (Vais 2010): P_Cl:P_K {pub['pcl_pk']:.2f}   P_Ca:P_K "
          f"{pub['pca_pk']:.1f}   i_Ca {pub['ica_slope']:.2f} pA/mM   "
          f"g {pub['conductance']:.0f} pS")
    _, ghk = sel.ghk_calcium_permeability(pub["conductance"] * 1e-12,
                                          pub["pcl_pk"], pub["pca_pk"],
                                          _P.value("permeation.bath_concentration"))
    print(f"GHK from the measured g and ratios predicts {ghk * 1e9:.2f} pA/mM "
          f"({ghk * 1e9 / pub['ica_slope']:.1f}x the measured slope)")
    return 0


def _titration_rows(wall, detail) -> None:
    lining = {(g.chain, g.res_seq) for g in wall.groups}
    names = sorted({(g.res_seq, g.res_name) for g in wall.groups})
    heads = [k for k in detail if k != "propka"]
    print("mean charge per lining residue (apparent pKa range over copies):")
    print("  residue   " + "".join(f"{h:>24s}" for h in heads)
          + ("         PROPKA pKa" if "propka" in detail else ""))
    for seq, name in names:
        cells = []
        for h in heads:
            t = detail[h]
            ii = [i for i, s in enumerate(t.sites)
                  if s.res_seq == seq and s.key in lining]
            pk = t.apparent_pka()[ii]
            cells.append(f"{t.charge[ii].mean():+6.2f} ({pk.min():5.1f}-"
                         f"{pk.max():5.1f})")
        line = f"  {name}{seq}  " + "".join(f"{c:>24s}" for c in cells)
        if "propka" in detail:
            v = [p for k, (_, p) in detail["propka"].items()
                 if k[1] == seq and k in lining]
            line += f"   {min(v):5.2f}-{max(v):5.2f}"
        print(line)


def _protonation(args) -> int:
    from .io import loader
    from .physics import protonation as pr
    from .physics.ryr_mutants import selectivity_mutants
    loader.ALLOW_FETCH = args.fetch
    st = loader.load(args.pdb)
    wall = pr.lining_wall(st)
    ph, ionic = pr.family_conditions(wall.ryr)
    rule = ("Xu 2006: 250 mM KCl, 10 mM CaCl2 luminal, their Eq. 1"
            if wall.ryr else "Vais 2010: 140 mM KCl, both ratios")
    print(f"{st.name}: lining groups titrated at pH {ph:g}, {ionic * 1e3:.0f} mM"
          f" ({rule})")
    rows, detail = pr.readings(wall, with_propka=not args.no_propka)
    _titration_rows(wall, detail)
    print("selectivity under each reading:")
    for r in rows:
        print("  " + r.row())
    if wall.ryr:
        print("Xu 2006's mutants (formal wall), P_Ca:P_K:")
        mut = selectivity_mutants(st)
        for m in mut:
            print("  " + m.row(mut[0]))
    if args.corners:
        out = pr.corners(wall)
        print("corners (each ring formal or neutral), largest P_Ca:P_K first:")
        for _, r in out[:5]:
            print("  " + r.row())
        bases = {g.res_seq for g in wall.groups if g.charge > 0}
        charged = [r for off, r in out if not off & bases]
        print(f"  best with every base ring charged: {charged[0].label} "
              f"P_Ca:P_K {charged[0].pca_pk:.2f}")
        inner = pr.interior(wall)
        print(f"  {len(inner)} interior states: largest P_Ca:P_K "
              f"{max(r.pca_pk for _, r in inner):.2f}")
    return 0


def register(sub) -> None:
    q = sub.add_parser("protonation", help="pKa of the lining groups (network "
                       "and PROPKA) and the selectivity under each reading")
    q.add_argument("pdb", nargs="?", default="8TKF")
    q.add_argument("--fetch", action="store_true")
    q.add_argument("--corners", action="store_true", help="also every ring "
                   "formal/neutral combination (a pKa-independent bound)")
    q.add_argument("--no-propka", action="store_true")
    q.set_defaults(fn=_protonation)

    p = sub.add_parser("selectivity", help="P_Cl:P_K, P_Ca:P_K and i_Ca from "
                       "the pore, against Vais 2010")
    p.add_argument("pdb", nargs="?", default="8TKF")
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--closure", choices=("donnan", "radial"), default="donnan",
                   help="how a charged slice is neutralised: uniformly (local "
                   "Donnan) or by Poisson-Boltzmann across it")
    p.set_defaults(fn=_selectivity)
