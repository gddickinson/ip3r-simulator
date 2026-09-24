"""E-C coupling: the V channel (Rios et al. 1993) and the voltage protocol
(Stern et al. 1997). Imported into ``parameter_table.P``.

Read from the supplied PDF of Rios, Karhanek, Ma & Gonzalez 1993 (J Gen
Physiol 102:449), Table I, on 2026-09-24. Fiber 827 is the row used: it is
the one the paper's kinetic simulations use (Fig. 9 B), and it is the row
that reproduces Stern 1997's V-channel open probabilities (Figs. 11-12:
0.46 / 0.05 / 0.003 at 0 / -30 / -50 mV; the other three reference rows
give 0.24-0.71 at 0 mV). Stern et al. took "the V channel parameters ...
from the original model of Rios et al. (1993), with time-dimensioned rate
constants increased by a factor of 2"; that factor is its own entry, so the
unscaled Table I values stay as printed. Their couplon figure (Fig. 12)
follows the factor; their stand-alone V check (Fig. 11) does not.
"""

from param_entry import entry as _p

_RIOS = "Rios et al. 1993 (supplied PDF), Table I, fiber 827 (reference)"
_STERN = "Stern, Pizarro & Rios 1997 (PMC2229377)"

EC = [
    _p("ec.rios_k", "Voltage-sensor steepness K", 4.5, "mV", "empirical",
       "ec", "rios1993", "Boltzmann steepness of one sensor, written 4K: "
       "Pa/Pr = exp((V - V_bar)/4K).", f"{_RIOS}: K 4.5 mV", 0.5, 50.0),
    _p("ec.rios_v_bar", "Voltage-sensor midpoint V_bar", -20.0, "mV",
       "empirical", "ec", "rios1993", "Voltage at which a sensor on a closed "
       "channel is equally resting and activating.", f"{_RIOS}: V_bar -20 mV",
       -120.0, 60.0),
    _p("ec.rios_k_l", "V-channel opening rate k_L", 2.0, "1/s", "empirical",
       "ec", "rios1993", "Opening rate with no sensor activated; each "
       "activated sensor divides it by f.", f"{_RIOS}: k_L 0.002 /ms", 1e-4,
       1e4),
    _p("ec.rios_k_minus_l", "V-channel closing rate k_-L", 9.0e5, "1/s",
       "empirical", "ec", "rios1993", "Closing rate with no sensor "
       "activated; each activated sensor multiplies it by f.",
       f"{_RIOS}: k_-L 900 /ms", 1.0, 1e8),
    _p("ec.rios_f", "Allosteric factor f", 0.175, "", "empirical", "ec",
       "rios1993", "Per activated sensor, the opening rate x 1/f and the "
       "closing rate x f (closing equilibrium x f^2).", f"{_RIOS}: f 0.175",
       0.01, 1.0),
    _p("ec.rios_alpha", "Sensor rate alpha", 200.0, "1/s", "empirical", "ec",
       "rios1993", "Inverse time constant of one sensor at V_bar: k_c = "
       "alpha/2 exp((V - V_bar)/8K), k_-c = alpha/2 exp(-(V - V_bar)/8K) "
       "(Eqs. 2-3).", "Rios 1993 Table I legend: 'the kinetic constant of "
       "charge movement ... was always 0.2 ms-1'", 1.0, 1e5),
    _p("ec.stern_rate_scale", "Stern's V-channel rate factor", 2.0, "",
       "empirical", "ec", "stern1997", "Every time-dimensioned rate of the "
       "Rios model (alpha, k_L, k_-L) is multiplied by this; equilibria are "
       "unchanged.", f"{_STERN}: 'with time-dimensioned rate constants "
       "increased by a factor of 2'. Their couplon runs (Fig. 12, V at 0 mV: "
       "0.43 at 10 ms) follow x2; their stand-alone V check (Fig. 11) follows "
       "the printed rates (x1), so it was made before the factor "
       "(tests/test_allosteric_v.py holds both)", 0.01, 100.0),
    _p("ec.v_unitary_current", "V-channel unitary Ca2+ current", 0.1, "pA",
       "empirical", "ec", "stern1997", "Current of one open V channel, for "
       "its Ca2+ at the C channels.", f"{_STERN}, Table I: 'i V V channel "
       "unitary current 0.1 pA'", 0.0, 10.0),
    _p("ec.holding", "Holding potential", -90.0, "mV", "convention", "ec",
       "stern1997", "Membrane potential before and after the pulse; the "
       "array starts from its stationary state here.", f"{_STERN}: release "
       "'during 100-ms depolarizations from a holding potential of -90 mV'",
       -150.0, 0.0),
    _p("ec.pulse", "Pulse duration", 0.1, "s", "convention", "ec",
       "stern1997", "Length of the depolarising step.", f"{_STERN}: "
       "100-ms step depolarizations (Figs. 1, 11, 12)", 1e-3, 10.0),
    _p("ec.after", "Record after the pulse", 0.1, "s", "method", "ec",
       "method_choice", "How long the array is followed after "
       "repolarisation.", "As Stern's Figs. 11-12 (200 ms traces)", 0.0, 10.0),
    _p("ec.plateau_window", "Plateau window", 0.04, "s", "method", "ec",
       "method_choice", "The plateau is the mean over this last part of the "
       "pulse.", "Stern's Fig. 12 plateaus are flat over the last ~60 ms of "
       "the 100-ms pulse; 40 ms keeps clear of the peak's tail at -50 mV",
       1e-3, 10.0),
    _p("ec.trials", "Couplons per ensemble", 200.0, "", "method", "ec",
       "method_choice", "Independent couplons averaged for an open "
       "probability or flux (Stern averaged 1,000-10,000).", "Enough for "
       "a peak to ~5 %; the CLI can raise it", 1.0, 1e5),
    _p("ec.spark_voltage", "Voltage for single events", -50.0, "mV",
       "convention", "ec", "stern1997", "Where V openings are sparse and C "
       "events stand apart, so their durations can be read.", f"{_STERN}, "
       "Fig. 14: the couplon at -50 mV, one V-triggered wave at a time",
       -120.0, 60.0),
    _p("ec.scan_v_min", "Voltage scan: lowest", -60.0, "mV", "method", "ec",
       "method_choice", "First step potential of the voltage scan.",
       "Stern's Fig. 1 B runs from -60 mV", -150.0, 60.0),
    _p("ec.scan_v_max", "Voltage scan: highest", 0.0, "mV", "method", "ec",
       "method_choice", "Last step potential of the voltage scan.",
       "Stern's Figs. 11-12 top at 0 mV", -150.0, 60.0),
    _p("ec.scan_v_step", "Voltage scan: step", 10.0, "mV", "method", "ec",
       "method_choice", "Spacing of the voltage scan.", "Stern's Fig. 1 B "
       "steps by 10 mV", 1.0, 100.0),
    # ------------------------------------------------ SR content (Round 6.7)
    _p("lumen.content", "SR Ca2+ content at rest", 2.0, "mM", "empirical",
       "ec", "stern1997", "Releasable SR Ca2+ as a concentration in "
       "accessible myoplasmic water; the unitary currents scale with the "
       "fraction left.", f"{_STERN}, Fig. 20: 'starting from an intra-SR "
       "calcium content equivalent to 2 mM in accessible myoplasmic water'. "
       "Rios 1993 (Fig. 2) measured 1.5 and 1.2 mM in fiber 827", 0.01, 100.0),
    _p("lumen.couplon_density", "Couplon density", 4.8, "1/um^3",
       "empirical", "ec", "stern1997", "Couplons per unit fibre volume: turns "
       "one couplon's release current into a rate of SR emptying.",
       f"{_STERN}: 'The density of couplons was 4.8 um-3' (junctional "
       "t-tubule area per Z disk from Peachey & Eisenberg 1978, times Z disks "
       "per length)", 0.01, 1e3),
    _p("lumen.refill_tau", "SR refilling time constant", 0.43, "s",
       "empirical", "ec", "stern1997", "Uptake returns the content to rest "
       "exponentially with this time constant (first order in the deficit).",
       f"{_STERN}, Fig. 20 C, digitised: after the pulse the content goes "
       "0.655 -> 0.929 mM (0 mV) and 1.34 -> 1.467 mM (-30 mV) in 96 ms, "
       "tau 421 and 450 ms. Their uptake is the Brum 1988 pump driven by "
       "global Ca2+, which runs faster during the pulse, so first order "
       "under-refills there", 1e-3, 1e3),
    _p("lumen.fig20_channels", "Fig. 20 couplon: C channels", 14.0, "",
       "convention", "ec", "stern1997", "C channels in the couplon of the "
       "depletion calibration (as many V channels again).", f"{_STERN}, "
       "Fig. 20: 'The dynamics of the model (28 channels/couplon)'", 1.0,
       500.0),
    _p("lumen.iterations", "SR path: most iterations", 12.0, "", "method",
       "ec", "method_choice", "The content path and the ensemble's release "
       "are made consistent by fixed-point iteration, at most this often.",
       "Converged in 4-7 on Stern's constants and up to 10 on the fitted scheme (Round 6.7)", 1.0,
       100.0),
    _p("lumen.tolerance", "SR path: tolerance", 0.01, "", "method", "ec",
       "method_choice", "Iteration stops when no point of the content "
       "fraction moves by more than this.", "1 % of the content: below the "
       "ensemble's sampling noise at 200 couplons", 1e-6, 1.0),
    _p("lumen.damping", "SR path: damping", 0.5, "", "method", "ec",
       "method_choice", "Each iteration moves the path this fraction of the "
       "way to the new one (1 = plain substitution).", "Plain substitution "
       "oscillates when the store empties within the pulse", 0.01, 1.0),
    _p("lumen.rios_released_low", "Released by a conditioning pulse: low",
       0.5, "", "empirical", "ec", "rios1993", "Fraction of SR content lost "
       "in a 100-ms pulse to +20 mV (fiber 827, 10 C).", "Rios 1993 Fig. 2 "
       "text: 'This depletion is ~60% for the pulse shown early in the "
       "experiment and ~50% in the late record'", 0.0, 1.0),
    _p("lumen.rios_released_high", "Released by a conditioning pulse: high",
       0.6, "", "empirical", "ec", "rios1993", "As the low value, early in "
       "the experiment.", "Rios 1993 Fig. 2 text (see the low value)", 0.0,
       1.0),
    _p("lumen.rios_voltage", "Conditioning pulse voltage", 20.0, "mV",
       "convention", "ec", "rios1993", "Step of Rios's conditioning pulse "
       "(100 ms from -80 mV).", "Rios 1993 Figs. 1-2: conditioning pulse to "
       "+20 mV, 100 ms", -150.0, 100.0),
    _p("lumen.skrap_max", "Local depletion in a spark: most", 0.074, "",
       "empirical", "ec", "launikonis2006", "Largest fall of free SR Ca2+ "
       "measured during a skeletal spark (frog, SEER).", "Launikonis et al. "
       "2006 PNAS (PMC1413852): 'implying a reduction of [Ca2+]SR by at most "
       "7.4%' (sulfate, sparks above the median)", 0.0, 1.0),
]
