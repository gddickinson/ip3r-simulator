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

_M97 = ("Meissner, Rios, Tripathy & Pasek 1997 (JBC 272:1628), Table IV: "
        "[3H]ryanodine binding, 0.5 M choline-Cl, fitted with their Eqs. 3-4 "
        "(K_eff^na = Ka^na (1 + ([I]/Ki)^ni)), mean +/- S.D.")

_L04 = ("Laver, O'Neill & Lamb 2004 (PMC2234024), Table I: RyR1 (rabbit) in "
        "bilayers, 2 mM ATP, luminal Ca2+ 1 mM, cytoplasmic Cs+ 250 mM")

_LL98_FULL = ("Laver & Lamb 1998, full text (rabbit skeletal and sheep "
              "cardiac RyRs in bilayers, 250 mM Cs+)")
_LL98 = ("Laver & Lamb 1998 (PMC1299578), abstract (the full text is a "
         "paywalled page scan): rabbit skeletal and sheep cardiac RyRs in "
         "bilayers, inactivation after steps in voltage and in cytosolic Ca2+")

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
    # ------------------------------------------------------------ Mg2+
    _p("ryr.mg_free", "Cytosolic free Mg2+ (muscle)", 1000.0, "uM",
       "empirical", "ryr", "laver2018", "Free Mg2+ in the resting fibre "
       "cytosol; the spark Mg2+ scan ends here.", "Laver 2018 review: "
       "'physiological concentrations of Mg2+ (1 mM (Godt and Maughan "
       "1988))'", 0.0, 1e5),
    _p("ryr.k_mg_a", "RyR1 activation site: Mg2+ affinity", 54.0, "uM",
       "empirical", "ryr", "laver2004", "Dissociation constant of Mg2+ at "
       "the Ca2+ activation (A-) site, competitive with Ca2+.", f"{_L04}: "
       "Kapp(Mg2+) 54 +/- 4 uM (DIDS-modified: 120 +/- 17 uM). Measured "
       "with ATP present, which Murayama's bell lacks; the Laver 2018 review "
       "gives '~50 uM' for the A-site", 1.0, 1e5),
    _p("ryr.mg_i_relative", "RyR1 inactivation site: Mg2+/Ca2+ affinity",
       1.0, "", "empirical", "ryr", "laver1997mh", "Mg2+ affinity of the "
       "low-affinity inactivation (I1-) site relative to Ca2+'s: the "
       "inactivation gate's ligand is c + ratio x Mg2+.", "Laver et al. "
       "1997 abstract: 'the inhibitory effects of Ca2+ and Mg2+ were "
       "virtually identical for the same conditions in any given channel'; "
       "Hill coefficient ~2 there, 1 in this gate", 0.0, 100.0),
    _p("ryr.k_ca_a_laver", "RyR1 activation site: Ca2+ affinity (Laver)",
       0.51, "uM", "empirical", "ryr", "laver2004", "The A-site Ca2+ "
       "affinity measured beside ryr.k_mg_a in the same condition; their "
       "ratio is the Mg2+/Ca2+ selectivity of the site.", f"{_L04}: "
       "Kapp(Ca2+)c 0.51 +/- 0.07 uM (ratio to Kapp(Mg2+) 106; DIDS 133)",
       0.01, 1e3),
    _p("ryr.meissner_ki_mg", "Mg2+ inhibition constant (Meissner)", 18.0,
       "uM", "empirical", "ryr", "meissner1997", "Competitive Mg2+ "
       "constant at the Ca2+ activation site, in the assay Murayama used.",
       f"{_M97}: Ca2+/Mg2+ with 5 mM AMP, Ki 0.018 +/- 0.009 mM (n = 4); "
       "without AMP 0.013 +/- 0.004 mM. The +AMP row because Murayama's "
       "bell had 1 mM AMP", 1.0, 1e4),
    _p("ryr.meissner_ni_mg", "Mg2+ Hill coefficient (Meissner)", 1.2, "",
       "empirical", "ryr", "meissner1997", "Hill coefficient of the Mg2+ "
       "term in Eq. 4.", f"{_M97}: +AMP 1.2 +/- 0.3; -AMP 1.1 +/- 0.1",
       0.1, 5.0),
    _p("ryr.meissner_na_ca", "Ca2+ activation Hill coefficient (Meissner)",
       2.0, "", "empirical", "ryr", "meissner1997", "The Ca2+ Hill "
       "coefficient na of the same fit; the shift of half-activation is the "
       "na-th root of Eq. 4's factor.", f"{_M97}: Ca2+/Mg2+ +AMP na 2.0 +/- "
       "0.2 (Ka 0.18 uM); -AMP 1.9 +/- 0.4 (Ka 0.39 uM)", 0.1, 5.0),
    _p("ryr.meissner_ki_na", "Na+ inhibition constant (Meissner)", 24000.0,
       "uM", "empirical", "ryr", "meissner1997", "Competitive Na+ constant "
       "at the same site: Murayama's NaCl occupies it too.", f"{_M97}: "
       "Ca2+/Na+ +AMP Ki 24 +/- 6 mM (n = 3); -AMP 27 +/- 7 mM", 1.0, 1e7),
    _p("ryr.meissner_ni_na", "Na+ Hill coefficient (Meissner)", 1.7, "",
       "empirical", "ryr", "meissner1997", "Hill coefficient of the Na+ "
       "term.", f"{_M97}: +AMP 1.7 +/- 0.2; -AMP 1.8 +/- 0.2", 0.1, 5.0),
    _p("ryr.murayama_sodium", "Na+ in Murayama's assay", 170000.0, "uM",
       "empirical", "ryr", "murayama2015", "The Na+ competing at the "
       "activation site in the condition the fitted bell was measured in.",
       _MUR, 0.0, 1e7),
    _p("spark.mg_scan_min", "Mg2+ scan: lowest non-zero Mg2+", 10.0, "uM",
       "method", "spark", "method_choice", "Low end of the free Mg2+ scan "
       "of cleft sparks (0 is always included).", "A fifth of the A-site "
       "affinity: where Mg2+ begins to matter", 0.1, 1e4),
    _p("spark.mg_scan_points", "Mg2+ scan: non-zero points", 6.0, "",
       "method", "spark", "method_choice", "Geometric steps from "
       "spark.mg_scan_min to ryr.mg_free.", "Two to three per decade",
       2.0, 50.0),
    _p("spark.trigger_trials", "Triggered spark: trials", 20.0, "",
       "method", "spark", "method_choice", "Seeds per triggered-spark "
       "measurement (every available channel opened at t = 0).", "Enough "
       "for a median; each trial is one spark", 1.0, 1000.0),
    _p("spark.trigger_window", "Triggered spark: window", 2.0, "s",
       "method", "spark", "method_choice", "How long a triggered array is "
       "followed; a spark still running then is counted as not ended.",
       "300x the 6.3 ms measured release, 10x the longest ended spark seen",
       0.01, 100.0),
    # ------------------------------- use-dependent (flux-driven) inactivation
    _p("ryr.k_use_on", "RyR1 use-dependent inactivation rate at 0 mV", 0.056,
       "s^-1", "empirical", "ryr", "laverlamb1998", "Rate at which a "
       "conducting RyR1 enters use-dependent inactivation, at the SR's "
       "resting potential (taken as 0 mV). Entered only from the open state, "
       "so it is Ca2+-independent. Only ryr.use_recovery_ratio decides the "
       "steady state, the fitted Ca2+ gate and spark termination; this rate "
       "sets how fast the gate relaxes, far slower than a spark.",
       f"{_LL98_FULL}, Fig. 4: log10 of the macroscopic inactivation rate "
       "against bilayer potential, cardiac RyRs at Po > 0.8, the fit at "
       "positive potentials having 'slope and intercept at V = 0' of 20 V^-1 "
       "and 0.056 s^-1. The intercept is an extrapolation: the data start "
       "near +20 mV. The abstract's 'tau approximately 1-3 s' is at +40 mV "
       "(the same line gives 0.35 s^-1 there, tau 2.8 s), which Round 6.9 "
       "took for this rate. Skeletal RyRs were not measured against voltage",
       1e-4, 1e4),
    _p("ryr.use_rate_slope", "Use-gate rate: voltage slope", 20.0, "V^-1",
       "empirical", "ryr", "laverlamb1998", "log10 of the inactivation rate "
       "per volt of bilayer potential (cytosol relative to lumen), positive "
       "limb. Used only to check the reading of Fig. 4 against the paper's "
       "own abstract and charge estimate.", f"{_LL98_FULL}, Fig. 4 legend; "
       "20 V^-1 in log10 is z delta 1.18, as the text's 1.14 +- 0.25",
       0.0, 100.0),
    _p("ryr.use_recovery_ratio", "RyR1 use-gate recovery / inactivation",
       0.2, "", "empirical", "ryr", "unverified", "k_use-/k_use: the use "
       "gate's recovery rate relative to its rate from the open state. The "
       "one number Round 6.9's result rests on: with both rates far slower "
       "than a spark, the steady state, the refitted Ca2+ gate and spark "
       "termination depend on this ratio and not on either rate.",
       f"{_LL98_FULL} report no recovery rate at a fixed potential: "
       "inactivated channels recovered only when the voltage was reversed "
       "('reactivated in milliseconds'), and at +40 mV 'remained closed "
       "until ... reactivated by a brief voltage pulse to -40 mV'. Their "
       "Fig. 8 bounds the ratio at +40 mV only (ryr.use_residual_40mv_*: at "
       "most 0.03-0.52). 0.2 lies inside that span; at 0 mV nothing is "
       "measured, so the value is unverified and the ratio is scanned",
       1e-4, 1e4),
    _p("ryr.use_residual_40mv_min", "Use gate: least residual at +40 mV",
       0.03, "", "empirical", "ryr", "laverlamb1998", "Ensemble activity 5 s "
       "after a step from -40 to +40 mV, relative to its peak: the most "
       "inactivated skeletal RyR. A residual R bounds the recovery ratio by "
       "R/(1-R) (Po <= 1, and 5 s is not quite steady state).",
       f"{_LL98_FULL}, Fig. 8 (skeletal, open circles, y axis), read from "
       "the figure to about +-0.02; three of the six skeletal channels did "
       "not inactivate (residual 1) and are the heterogeneity, not this "
       "bound", 0.0, 1.0),
    _p("ryr.use_residual_40mv_max", "Use gate: largest residual at +40 mV",
       0.34, "", "empirical", "ryr", "laverlamb1998", "As "
       "ryr.use_residual_40mv_min, for the least inactivated skeletal RyR "
       "that inactivated at all.", f"{_LL98_FULL}, Fig. 8, read to about "
       "+-0.02", 0.0, 1.0),
    _p("ryr.use_inactivating_fraction", "RyR1 share carrying the use gate",
       0.8, "", "empirical", "ryr", "laverlamb1998", "Share of a cluster's "
       "RyR1 channels that have use-dependent inactivation at all; the rest "
       "never enter it. A stable property of the channel: 'only those that "
       "showed inactivation after voltage steps inactivated after [Ca2+] "
       "steps, and vice versa'. Read by ryr_mixed.",
       f"{_LL98_FULL}, Results: of channels with Po above 0.2, '80% of "
       "skeletal RyRs (12 of 15)' inactivated within 30 s of voltage steps "
       "(cardiac 56 %, 33 of 59); the abstract's 'one-half to two-thirds' "
       "pools both. Laver & Curtis 1996 (Biophys J 71:732): 70 % of 25 "
       "cardiac and skeletal channels declined after Ca2+ steps; Ma 1995: "
       "50-70 % of skeletal; Sitsapesan et al. 1995 (cardiac, +40 mV): "
       "17 % (4 of 24) with Ca2+ alone, 56 % (5 of 9) with ATP or EMD "
       "41000, the Po dependence Laver & Lamb's Po > 0.2 selection removes. "
       "The skeletal count is used; the fraction is scanned from 0 to 1",
       0.0, 1.0),
    _p("ryr.sitsapesan_po_early_m40", "RyR2 Po early at -40 mV", 0.841, "",
       "empirical", "ryr", "sitsapesan1995", "Open probability just after a "
       "rapid step to activating ligands, at -40 mV: the start of a 5-s "
       "record in which no inactivation was seen. Checks the use-gate rate "
       "read from Laver & Lamb's Fig. 4.", "Sitsapesan et al. 1995, Results and Fig. 2: sheep cardiac RyRs activated by 100 uM Ca2+ plus 1 mM ATP, held at -40 mV cis (cytosol) relative to trans; ensemble Po in the first and last part of a 5-s exposure, mean +- SEM, n = 5: Po early 0.841 +- 0.087",
       0.0, 1.0),
    _p("ryr.sitsapesan_po_late_m40", "RyR2 Po late at -40 mV", 0.874, "",
       "empirical", "ryr", "sitsapesan1995", "As "
       "ryr.sitsapesan_po_early_m40, at the end of the 5-s exposure.",
       "Sitsapesan et al. 1995, Results and Fig. 2: sheep cardiac RyRs activated by 100 uM Ca2+ plus 1 mM ATP, held at -40 mV cis (cytosol) relative to trans; ensemble Po in the first and last part of a 5-s exposure, mean +- SEM, n = 5: Po late 0.874 +- 0.072", 0.0, 1.0),
    _p("ryr.sitsapesan_po_sem_m40", "RyR2 Po at -40 mV: larger SEM", 0.087,
       "", "empirical", "ryr", "sitsapesan1995", "The larger of the two "
       "SEMs, used for both when the late/early ratio's uncertainty is "
       "formed (as if independent: conservative, the pairs are one channel "
       "each).", "Sitsapesan et al. 1995, Results and Fig. 2: sheep cardiac RyRs activated by 100 uM Ca2+ plus 1 mM ATP, held at -40 mV cis (cytosol) relative to trans; ensemble Po in the first and last part of a 5-s exposure, mean +- SEM, n = 5: SEM 0.087 (early), 0.072 (late)", 0.0, 1.0),
    _p("ryr.sitsapesan_record_m40", "RyR2 record at -40 mV: duration", 5.0,
       "s", "empirical", "ryr", "sitsapesan1995", "How long the channels "
       "were held open at -40 mV without a decline in Po.",
       "Sitsapesan et al. 1995, Fig. 2: 'over several seconds', 5-s "
       "exposures to the activating solution", 0.0, 1e3),
    _p("spark.use_fraction_scan_points", "Inactivating-fraction scan: points",
       6.0, "", "method", "spark", "method_choice", "Evenly spaced "
       "inactivating fractions from 0 to 1 in the mixed-cluster scan.",
       "Steps of 0.2: the measured 0.5-0.8 falls on two of them", 2.0, 50.0),
    _p("spark.use_scan_min", "Use-gate scan: fastest time constant", 0.001,
       "s", "method", "spark", "method_choice", "Low end of the scan of the "
       "use gate's time constant from the open state; also the fastest use "
       "gate the couplon is run with.", "1 ms: below the 6.3 ms measured "
       "release, so a gate this fast could end a spark on its own", 1e-6,
       1e3),
    _p("spark.use_scan_max", "Use-gate scan: slowest time constant", 10.0,
       "s", "method", "spark", "method_choice", "High end of the scan.",
       "Beyond Laver & Lamb's measured 1-3 s, so the measured range lies "
       "inside the scan", 1e-4, 1e4),
    _p("spark.use_scan_points", "Use-gate scan: points", 9.0, "", "method",
       "spark", "method_choice", "Geometric steps of the use-gate scan.",
       "Two per decade over the four decades from 1 ms to 10 s", 2.0, 50.0),
    _p("spark.use_ratio_scan_max", "Use-gate ratio scan: widest factor",
       10.0, "", "method", "spark", "method_choice", "Largest multiple of "
       "ryr.use_recovery_ratio in the ratio scan; the scan runs from its "
       "inverse to it. The ratio is load-bearing and unmeasured at 0 mV, so "
       "its effect is measured rather than assumed.", "Two decades around "
       "the unverified value: brackets a cluster that never fires and one "
       "that never stops", 1.0, 1e3),
    _p("spark.use_ratio_scan_points", "Use-gate ratio scan: points", 9.0, "",
       "method", "spark", "method_choice", "Geometric steps of the ratio "
       "scan.", "Four per decade: the working band is under a decade wide",
       2.0, 50.0),
]
