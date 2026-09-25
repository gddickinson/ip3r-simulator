"""Permeation commands beyond the K+ conductance, ``register(sub)``:

    python -m ip3r selectivity [8TKF]   # Vais 2010's P_Cl:P_K, P_Ca:P_K, i_Ca
    python -m ip3r selectivity --closure radial   # PB across each slice
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


def register(sub) -> None:
    p = sub.add_parser("selectivity", help="P_Cl:P_K, P_Ca:P_K and i_Ca from "
                       "the pore, against Vais 2010")
    p.add_argument("pdb", nargs="?", default="8TKF")
    p.add_argument("--fetch", action="store_true")
    p.add_argument("--closure", choices=("donnan", "radial"), default="donnan",
                   help="how a charged slice is neutralised: uniformly (local "
                   "Donnan) or by Poisson-Boltzmann across it")
    p.set_defaults(fn=_selectivity)
