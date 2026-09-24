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

# ------------------------------------------------------------ RyR1 gating
_STERN = "Stern, Pizarro & Rios 1997 (PMC2229377), Table I"
_MUR = ("Murayama et al. 2015 (PMC4482644), S1 Table, wild type at 25 C "
        "([3H]ryanodine binding, 0.17 M NaCl, 1 mM AMP, no Mg2+)")
_MUR37 = _MUR.replace("at 25 C", "at 37 C") + "; read from the S1 Table .doc"

RYR += [
    _p("ryr.k_act_on", "RyR1 activation on rate", 10.0, "uM^-2 s^-1",
       "empirical", "ryr", "stern1997", "Two-Ca2+ activation of the C "
       "channel: rate k_o x c^2.", f"{_STERN}: 'k o Activating on rate 10^13 "
       "M-2 s-1' = 10 uM^-2 s^-1. The authors: the model 'has not been "
       "objectively fitted to data'", 0.01, 1e4),
    _p("ryr.k_act_off", "RyR1 activation off rate", 500.0, "s^-1",
       "empirical", "ryr", "stern1997", "Closing of the activation gate.",
       f"{_STERN}: 'k o- Activating off rate 500 s-1'. With k_o this gives "
       "Ka = 7.1 uM; the text says 'the activating site ... was 10 uM'",
       1.0, 1e5),
    _p("ryr.k_inact_on", "RyR1 inactivation on rate", 2.0, "uM^-1 s^-1",
       "empirical", "ryr", "stern1997", "One-Ca2+ inactivation: rate k_i "
       "x c.", f"{_STERN} prints 'k i Inactivating on rate 2 x 10-6 M-1s-1', "
       "a sign typo: the text gives 'the K d for inactivation, also 10 uM', "
       "and 20 s-1 / 10 uM = 2 x 10^6 M-1 s-1 = 2 uM^-1 s^-1", 0.001, 1e3),
    _p("ryr.k_inact_off", "RyR1 inactivation off rate", 20.0, "s^-1",
       "empirical", "ryr", "stern1997", "Recovery from inactivation.",
       f"{_STERN}: 'k i- Inactivating off rate 20 s-1'", 0.01, 1e4),
    _p("ryr.murayama_amax", "RyR1 bell: Amax", 0.031, "", "empirical", "ryr",
       "murayama2015", "Peak binding activity (B/Bmax) of the fitted bell.",
       f"{_MUR}: 'WT 25 ... 0.031 +/- 0.002'", 1e-4, 1.0),
    _p("ryr.murayama_ka", "RyR1 bell: KA", 5.5, "uM", "empirical", "ryr",
       "murayama2015", "Half-activation of the fitted bell.",
       f"{_MUR}: KA 5.5 +/- 0.9 uM", 0.01, 1000.0),
    _p("ryr.murayama_na", "RyR1 bell: nA", 1.2, "", "empirical", "ryr",
       "murayama2015", "Hill coefficient of activation, fixed in the fit.",
       "Murayama 2015 Methods: 'we used fixed values for n A (1.2) and n I "
       "(1.5) for WT and all the mutants'", 0.1, 10.0),
    _p("ryr.murayama_ki", "RyR1 bell: KI", 270.0, "uM", "empirical", "ryr",
       "murayama2015", "Half-inhibition of the fitted bell.",
       f"{_MUR}: KI 0.27 +/- 0.03 mM", 1.0, 1e5),
    _p("ryr.murayama_ni", "RyR1 bell: nI", 1.5, "", "empirical", "ryr",
       "murayama2015", "Hill coefficient of inhibition, fixed in the fit.",
       "Murayama 2015 Methods: 'fixed values for n A (1.2) and n I (1.5)'",
       0.1, 10.0),
    _p("ryr.murayama_amax_37", "RyR1 bell at 37 C: Amax", 0.107, "",
       "empirical", "ryr", "murayama2015", "Peak binding activity of the "
       "wild-type bell at 37 C.", f"{_MUR37}: 'WT 37 ... 0.107 +/- 0.006'",
       1e-4, 1.0),
    _p("ryr.murayama_ka_37", "RyR1 bell at 37 C: KA", 20.5, "uM", "empirical",
       "ryr", "murayama2015", "Half-activation of the wild-type bell at "
       "37 C.", f"{_MUR37}: KA 20.5 +/- 2.5 uM", 0.01, 1000.0),
    _p("ryr.murayama_ki_37", "RyR1 bell at 37 C: KI", 410.0, "uM",
       "empirical", "ryr", "murayama2015", "Half-inhibition of the wild-type "
       "bell at 37 C.", f"{_MUR37}: KI 0.41 +/- 0.04 mM", 1.0, 1e5),
    # -------------------------------------------------------- spark cluster
    _p("spark.n_channels", "Spark cluster size", 30.0, "", "empirical",
       "spark", "stern1997", "Ca2+-gated RyR1s in one cluster.", "Stern 1997: "
       "'the total number of channels, 2N, was 60', half V (voltage-coupled) "
       "and half C (Ca2+-gated) channels; the C channels are the ones a "
       "Ca2+-driven spark recruits", 1.0, 500.0),
    _p("spark.unitary_current", "RyR1 unitary Ca2+ current", 0.3, "pA",
       "empirical", "spark", "stern1997", "Current of one open channel, "
       "for the derived coupling.", f"{_STERN}: 'i C C channel unitary "
       "current 0.3 pA'", 0.01, 10.0),
    _p("spark.d_ca", "Ca2+ diffusion coefficient", 5e-6, "cm^2/s",
       "physical", "spark", "stern1997", "Free Ca2+ diffusivity, for the "
       "derived coupling.", f"{_STERN}: 'D ca Ca2+ diffusion coefficient 5 x "
       "10-6 cm2 s-1'", 1e-7, 1e-4),
    _p("spark.channel_spacing", "RyR1 channel spacing", 30.0, "nm",
       "empirical", "spark", "stern1997", "Distance at which the derived "
       "coupling is evaluated.", f"{_STERN}: 'Channel spacing 30 nm'", 5.0,
       500.0),
    _p("spark.dt", "Spark cluster Ca2+ update step", 2.5e-5, "s", "method",
       "spark", "method_choice", "How often the cluster Ca2+ is updated; "
       "within a step the transition matrix is exact.", "Measured here (8 "
       "seeds x 20 s): 1e-5 and 2.5e-5 s agree within noise (spark rate "
       "1.56/1.48 s-1, median duration 125/116 ms), while 1e-4 s lengthens "
       "sparks by 15-20 %", 1e-6, 1e-2),
    _p("spark.published_release_duration", "Spark release duration (frog)",
       6.3, "ms", "empirical", "spark", "rios1999", "Measured open time of "
       "the release underlying a Ca2+ spark, compared with the model's "
       "spark duration.", "Rios et al. 1999 (PMC2229636): 'a release "
       "current of 16.9 pA, coming from a source of 0.5 um, with an open "
       "time of 6.3 ms'. Frog skeletal muscle, not mammalian RyR1: a scale, "
       "not a target", 0.1, 1000.0),
    _p("spark.scan_low", "Spark scan: lowest coupling factor", 0.03, "",
       "method", "spark", "method_choice", "Low end of the spark coupling "
       "scan, as a multiple of the derived coupling.", "Buffers reduce the "
       "free-diffusion estimate; a 30-fold band below it covers strong "
       "buffering", 0.001, 1.0),
    _p("spark.scan_high", "Spark scan: highest coupling factor", 1.0, "",
       "method", "spark", "method_choice", "High end of the spark coupling "
       "scan, as a multiple of the derived coupling.", "The derived value "
       "is an unbuffered upper estimate, so the scan stops at it", 0.01,
       10.0),
]

# ------------------------------------------ the junctional cleft (spatial)
_CLEFT = ("Stern, Pizarro & Rios 1997 (PMC2229377), Table I and Fig. 7 B "
          "(double row, V and C channels alternating)")

RYR += [
    _p("spark.cleft_width", "Couplon width", 60.0, "nm", "empirical",
       "spark", "stern1997", "Width of the junctional cleft: two rows of "
       "channels, one channel spacing apart.", f"{_CLEFT}: 'Couplon width "
       "60 nm'", 10.0, 1000.0),
    _p("spark.cleft_height", "Junctional gap thickness", 15.0, "nm",
       "empirical", "spark", "stern1997", "Height of the cleft; diffusion "
       "is two-dimensional across it.", f"{_CLEFT}: 'Junctional gap "
       "thickness 15 nm'", 1.0, 200.0),
    _p("spark.source_diameter", "Release source diameter", 30.0, "nm",
       "empirical", "spark", "stern1997", "Each open channel releases "
       "uniformly over a disc of this diameter (the foot).", "Stern 1997 "
       "Appendix: 'we treated Ca2+ release as a diffuse source spread over "
       "a 30-nm-diameter disk representing the foot process'", 1.0, 100.0),
    _p("spark.r_max", "Edge boundary: background distance", 1000.0, "nm",
       "empirical", "spark", "stern1997", "Distance at which Ca2+ escaping "
       "the cleft edge reaches background; sets the edge transfer "
       "coefficient D_inf = 2 pi D / (h ln(R_max / h)).", "Stern 1997 "
       "Appendix, Eq. 13: 'R max, which was taken as the length scale of "
       "the junctional strip (1 um)'", 20.0, 1e5),
    _p("spark.cleft_grid", "Cleft finite-volume grid spacing", 1.0, "nm",
       "method", "spark", "method_choice", "Cell size of the steady-state "
       "diffusion solve that gives the coupling coefficients.", "Measured "
       "here: the coupling coefficients at 1 nm and 0.5 nm agree to < 1 %",
       0.1, 10.0),
    # ------------------------------------ what terminates a cleft spark
    _p("spark.ki_scan_min", "Termination scan: lowest Ki", 3.0, "uM",
       "method", "spark", "method_choice", "Low end of the inactivation "
       "constant scan (on-rate held at Stern's k_i, off-rate moved).",
       "Below Stern's 10 uM, to see whether a still more sensitive "
       "inactivation shortens sparks", 0.1, 1e4),
    _p("spark.ki_scan_max", "Termination scan: highest Ki", 400.0, "uM",
       "method", "spark", "method_choice", "High end of the inactivation "
       "constant scan.", "Covers Murayama's measured KI at 25 C (270 uM) "
       "and 37 C (410 uM)", 1.0, 1e5),
    _p("spark.ki_scan_points", "Termination scan: Ki points", 9.0, "",
       "method", "spark", "method_choice", "Geometric steps of the Ki scan.",
       "Two to three per decade", 2.0, 50.0),
    _p("spark.rate_scan_max", "Termination scan: fastest inactivation", 30.0,
       "", "method", "spark", "method_choice", "Largest multiple of Stern's "
       "inactivation rates (on and off together, Ki fixed) in the rate scan; "
       "the scan runs from its inverse to it.", "Steady state fixes only the "
       "ratio k_i-/k_i, so the rate is scanned over 3 decades", 1.0, 1e3),
    _p("spark.rate_scan_points", "Termination scan: rate points", 7.0, "",
       "method", "spark", "method_choice", "Geometric steps of the rate "
       "scan.", "Two to three per decade", 2.0, 50.0),
]
