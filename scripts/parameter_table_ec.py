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
]
