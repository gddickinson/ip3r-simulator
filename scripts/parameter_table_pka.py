"""Protonation in the pore (Round 7.4): model-compound pKas, the network's
electrostatics and sampling, the recording pHs, and RyR1's measured
P_Ca:P_K (the control). Imported into ``parameter_table.P``.

Read on 2026-09-25: Thurlkill 2006 and Fitch 2015 from their abstracts
(PMC2242523, PMC4420524); Vais 2010 (PMC2995152) and Xu 2006 (PMC1367051)
from Methods and Table 2.
"""

from param_entry import entry as _p

_THUR = ("Thurlkill et al. 2006 (PMC2242523), abstract: alanine pentapeptides "
         "with blocked termini, 0.1 M KCl, 25 C")
_MS = ("Mehler & Solmajer 1991's sigmoidal eps(r). The paper itself was not "
       "reachable; the constants are the ones AutoDock 4 implements for it "
       "(A = -8.5525, eps0 = 78.4, k = 7.7839, lambda = 0.003627 /A), "
       "checked here only as eps(r) running from ~1.3 at contact to 78.4")
_XU = ("Xu et al. 2006 (PMC1367051), Table 2: symmetric 250 mM KCl with 10 mM "
       "Ca2+ trans, read with their Eq. 1 (GHK, no Cl- term)")


def _model(res, value, citation, note):
    return _p(f"pka.model_{res.lower()}", f"Model pKa, {res}", value, "",
              "empirical", "pka", citation, f"{res} side chain in an "
              "unstructured model compound: the pKa a site has before any "
              "interaction.", note, 0.0, 16.0)


def _xu(name, value):
    return _p(f"selectivity.published_ryr1_pca_pk_{name.lower()}",
              f"Measured P_Ca:P_K, RyR1-{name}", value, "", "empirical",
              "selectivity", "xu2006", f"Recombinant rabbit RyR1-{name} in "
              "planar bilayers.", _XU, 0.0, 100.0)


PKA = [
    _model("ASP", 3.67, "thurlkill2006", _THUR + ": 3.67 +/- 0.04"),
    _model("GLU", 4.25, "thurlkill2006", _THUR + ": 4.25 +/- 0.05"),
    _model("HIS", 6.54, "thurlkill2006", _THUR + ": 6.54 +/- 0.04"),
    _model("LYS", 10.40, "thurlkill2006", _THUR + ": 10.40 +/- 0.08"),
    _model("ARG", 13.8, "fitch2015", "Fitch et al. 2015 (PMC4420524), "
           "abstract: 13.8 +/- 0.1 by potentiometry and NMR (not the ~12 "
           "often used)"),
    _p("pka.ph_vais", "pH, Vais 2010 recordings", 7.3, "", "empirical",
       "pka", "vais2010", "pH of every solution in the IP3R-3 selectivity "
       "and i_Ca experiments.", "Vais 2010 Methods: '10 mM HEPES, pH to 7.3 "
       "with KOH' (NMDG, NaOH in the substitution solutions)", 0.0, 14.0),
    _p("pka.ph_xu", "pH, Xu 2006 recordings", 7.4, "", "empirical", "pka",
       "xu2006", "pH of the RyR1 bilayer solutions.", "Xu 2006 Methods: "
       "'20 mM KHepes, pH 7.4'", 0.0, 14.0),
    _p("pka.eps_water", "Water permittivity (sigmoidal limit)", 78.4, "",
       "physical", "pka", "mehler1991", "Long-range limit of the "
       "distance-dependent permittivity, and the Debye length's solvent.",
       _MS, 1.0, 100.0),
    _p("pka.ms_a", "Mehler-Solmajer A", -8.5525, "", "method", "pka",
       "mehler1991", "Offset of the sigmoidal permittivity.", _MS, -20.0, 0.0),
    _p("pka.ms_k", "Mehler-Solmajer k", 7.7839, "", "method", "pka",
       "mehler1991", "Shape constant of the sigmoidal permittivity.", _MS,
       0.1, 100.0),
    _p("pka.ms_lambda", "Mehler-Solmajer lambda", 0.003627, "1/A", "method",
       "pka", "mehler1991", "Rate constant of the sigmoidal permittivity.",
       _MS, 1e-4, 0.1),
    _p("pka.site_radius", "Network radius", 20.0, "A", "method", "pka",
       "method_choice", "Titratable sites within this distance of the "
       "groups asked about are titrated with them.", "At 20 A a unit pair "
       "interacts by 0.02-0.03 kT (eps 76, Debye 8-6 A); "
       "tests/test_protonation.py shows 25 A moves no lining charge by "
       "0.01 e. The uniform eps-4 bound is long-ranged (0.6 kT at 20 A) and "
       "moves by < 0.04 e, which leaves P_Ca:P_K unchanged at 0.01",
       5.0, 60.0),
    _p("pka.mc_sweeps", "Monte Carlo sweeps", 3000, "", "method", "pka",
       "method_choice", "Sweeps of the network titration.", "Held to exact "
       "enumeration to 0.01 in tests/test_pka.py; on 8TKF two seeds agree to "
       "0.02 e even in the eps-4 bound", 10, 1000000),
    _p("pka.mc_burn_fraction", "Monte Carlo burn-in", 0.1, "", "method",
       "pka", "method_choice", "Fraction of the sweeps discarded before "
       "averaging.", "The start (each site at its own Henderson-Hasselbalch "
       "state) is near equilibrium; 10 % is generous", 0.0, 0.9),
    _p("pka.mc_seed", "Monte Carlo seed", 1, "", "method", "pka",
       "method_choice", "Random seed of the network titration.",
       "Fixed so a run is reproducible", 0, 2 ** 31),
    _p("pka.mc_pair_coupling", "Pair-move coupling", 1.0, "kT", "method",
       "pka", "method_choice", "Pairs coupled more strongly than this also "
       "get joint flips.", "Changes only how fast the sampler mixes, not "
       "what it samples (the exact-enumeration test holds either way)",
       0.0, 100.0),
    _p("pka.propka_radius", "PROPKA context radius", 25.0, "A", "method",
       "pka", "method_choice", "Whole residues with any atom this close to a "
       "lining group are written out for PROPKA.", "PROPKA 3's longest "
       "range is its 15 A desolvation count; tests/test_pka.py shows 25 and "
       "35 A give the same pKa", 16.0, 80.0),
    _p("protonation.interior_samples", "Interior samples", 8, "", "method",
       "pka", "method_choice", "Fractional wall states drawn inside the "
       "corners' box to check the corner bound.", "A check, not an "
       "estimate: 8 states take ~30 s on 8TKF", 0, 1000),
    _p("selectivity.ryr1_cacl2_lumen", "RyR1 P_Ca:P_K luminal CaCl2", 0.010,
       "M", "empirical", "selectivity", "xu2006", "CaCl2 added to the trans "
       "(luminal) side for the RyR1 reversal potential.", _XU, 1e-4, 0.1),
    _p("selectivity.published_ryr1_pca_pk", "Measured P_Ca:P_K, RyR1",
       7.0, "", "empirical", "selectivity", "xu2006", "Recombinant rabbit "
       "RyR1 (wild type) in planar bilayers.", _XU + ": E_rev 9.5 +/- 0.2 mV "
       "(n = 6) gives 6.96", 0.0, 100.0),
    _xu("D4899Q", 1.0),
    _xu("E4900N", 4.5),
    _xu("D4938N", 3.3),
    _xu("D4945N", 6.5),
    _xu("E4955Q", 8.3),
]
