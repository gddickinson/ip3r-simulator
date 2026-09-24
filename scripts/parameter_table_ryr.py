"""RyR1: the recording condition and the measured K+ conductances the pore
model is compared with. Imported into ``parameter_table.P``.

Read from Xu et al. 2006 (PMC1367051) on 2026-09-24: Table 2 and Methods,
"Single channel recordings". Recombinant rabbit RyR1 in planar lipid
bilayers, symmetric 250 mM KCl, 20 mM KHepes pH 7.4. The paper states no
bilayer temperature, so ``permeation.temperature`` (room temperature) is
used unchanged.
"""

from param_entry import entry as _p

_XU = "Xu et al. 2006 (PMC1367051), Table 2, gamma K+ in symmetric 250 mM KCl"


def _mutant(name, value, sd, n):
    return _p(f"permeation.published_ryr1_{name.lower()}",
              f"RyR1 {name} conductance", value, "pS", "empirical",
              "permeation", "xu2006", f"Recombinant rabbit RyR1-{name}, a "
              "charge-neutralising pore mutant, symmetric 250 mM KCl.",
              f"{_XU}: {value:g} +/- {sd:g} pS (n = {n})", 1.0, 3000.0)


RYR = [
    _p("permeation.ryr1_bath_concentration", "RyR1 bath KCl concentration",
       0.25, "M", "convention", "permeation", "xu2006", "Symmetric bath salt "
       "of the RyR1 recordings the model is compared with.", "Xu 2006 "
       "Methods: 250 mM KCl cis, trans raised to 250 mM after fusion", 0.001,
       3.0),
    _p("permeation.published_ryr1", "RyR1 conductance (wild type)", 801.0,
       "pS", "empirical", "permeation", "xu2006", "Recombinant rabbit RyR1 in "
       "planar bilayers, symmetric 250 mM KCl.",
       f"{_XU}: 801 +/- 7 pS (n = 17); the text rounds it to ~800", 1.0,
       3000.0),
    _mutant("D4899Q", 164.0, 4, 10),
    _mutant("E4900N", 505.0, 3, 3),
    _mutant("D4938N", 520.0, 6, 9),
    _mutant("D4945N", 737.0, 11, 6),
    _mutant("E4955Q", 812.0, 7, 4),
]
