"""Charge-space competition (Round 7.17): the hard-sphere and MSA excess
chemical potential in the 1-D pore. Imported into ``parameter_table.P``.

Read on 2026-09-26: Gillespie 2008 (PMC2212702), Theory and Methods and
the Appendix; Nonner, Catacuzzeno & Eisenberg 2000 (PMC1301088), abstract.
The ion diameters are twice the registered ``permeation.radius_*`` (the
Shannon/Pauling radii), which are Gillespie's (K+ 2.76, Ca2+ 2.00,
Cl- 3.62 A), so they are not registered twice.
"""

from param_entry import entry as _p

_G08 = ("Gillespie 2008 (PMC2212702), Theory and Methods")

CSC = [
    _p("csc.permittivity", "Permittivity of the MSA fluid", 78.4, "",
       "physical", "csc", "gillespie2008", "Relative permittivity in the "
       "Bjerrum length of the screening (MSA) term.", _G08 + ": 'the "
       "dielectric coefficient was constant at 78.4 throughout the system'. "
       "The mean field of the 1-D model is local electroneutrality, which "
       "has no permittivity; `csc --scan` lowers this to 40 "
       "(permeation.permittivity_pore)", 2.0, 100.0),
    _p("csc.oxygen_diameter", "Carboxylate oxygen diameter", 2.8, "A",
       "physical", "csc", "gillespie2008", "Each lining Asp/Glu enters the "
       "fluid as two oxygens of this diameter, carrying half its charge "
       "each.", _G08 + ": 'two independent, half-charged oxygen ions (2.8 A "
       "diameter)'", 1.0, 5.0),
    _p("csc.base_diameter", "Base group diameter", 2.8, "A", "method",
       "csc", "method_choice", "Each lining Lys/Arg enters the fluid as one "
       "sphere of this diameter carrying its charge.", "Gillespie's filter "
       "has no base, so there is no calibrated value; the oxygen's diameter "
       "is taken (an NH3+ or guanidinium nitrogen with its hydrogens is "
       "of that size)", 1.0, 6.0),
    _p("csc.water_diameter", "Water diameter", 2.8, "A", "physical", "csc",
       "gillespie2008", "Water is an uncharged hard sphere of this diameter "
       "in every slice, displaced where ions crowd in.", _G08 + ": 'water "
       "as an uncharged, hard sphere'. The paper gives no diameter; 2.8 A "
       "(the oxygen's, and the O-O distance of liquid water) is the value "
       "the Boda/Gillespie models use", 1.0, 4.0),
    _p("csc.water_concentration", "Bath water concentration", 55.5, "M",
       "physical", "csc", "convention", "Bath water; its activity sets "
       "the water in each slice.", "Pure water at 25 C (997 g/L / 18.015 "
       "g/mol); packing fraction 0.38 at 2.8 A", 0.0, 60.0),
    _p("csc.structural_volume", "Wall groups' volume counted", 1.0, "",
       "method", "csc", "method_choice", "1: the wall's oxygens and bases "
       "take up room in the lumen (Gillespie's tethered groups); 0: they "
       "screen but take none.", "The deposit's free radius already excludes "
       "those atoms, so 1 counts their volume twice and 0 not at all; both "
       "are reported (`csc`)", 0.0, 1.0),
    _p("csc.filter_reach", "Filter slice search", 6.0, "A", "method", "csc",
       "method_choice", "The binding energetics are read at the most "
       "charged slice within this distance of the filter constriction.",
       "RyR1 9HEO's D4899/E4900 peak sits 4 A luminal of the filter "
       "minimum; 6 A reaches it without reaching the next ring", 0.5, 20.0),
    _p("csc.max_step", "Largest Newton step", 2.0, "", "method", "csc",
       "method_choice", "Largest change of any log concentration (or of "
       "psi in thermal voltages) one Newton step of the local equilibrium "
       "may take.", "Picard iteration on the excess ran away (100 M K+ "
       "with no water) at the RyR1 filter; capped Newton with backtracking "
       "converges there in tens of steps", 0.1, 20.0),
    _p("csc.tolerance", "Partition tolerance", 1e-10, "", "method", "csc",
       "method_choice", "Largest residual (log activity, or the relative "
       "charge imbalance) at which the local equilibrium stops.", "1e-12 moves no reported ratio in "
       "its third figure", 1e-14, 1e-4),
    _p("csc.max_iterations", "Partition iterations", 200, "", "method",
       "csc", "method_choice", "Newton steps allowed per local "
       "equilibrium.", "The Gummel loop warm-starts each pass from the "
       "last, so after the first few it takes tens", 10, 100000),
]
