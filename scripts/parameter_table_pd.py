"""The park/drive receptor (Siekmann 2012 as extended by Cao et al. 2013) and
the puff-comparison method parameters. Imported into ``parameter_table.P``.

Every ``pd.*`` constant was read from the authors' own model code, Cao et al.
2014 Text S1 (the stochastic 6-state script), on 2026-09-24. The main text of
Cao et al. 2013 (PMC3852038) gives the equations; its Table S1 was not
reachable. The IP3 dependences are printed there as, e.g.,
``V24 = 62 + 880/(p^2 + 4)``. They are registered in the one Hill form the
module evaluates, ``base + amp * H(p)`` with ``H = p^h/(p^h + K^h)``
(rising) or ``K^h/(p^h + K^h)`` (falling). So ``880/(p^2+4)`` is
``amp = 880/4 = 220``, ``K = 2``, ``h = 2``, falling, and
``5/(p^2 + 0.5^2)`` is ``amp = 5/0.25 = 20``, ``K = 0.5``, ``h = 2``,
falling. Nothing else was transformed.
"""

from param_entry import entry as _p


_CODE = ("Cao et al. 2014 Text S1 (6-state stochastic model code), the "
         "IP3R-1 fit of Siekmann et al. 2012 at 0.1 mM ATP as extended by "
         "Cao et al. 2013; read 2026-09-24")


def _rate(key, a, b, value, what):
    return _p(f"pd.{key}", f"Rate q{key[1:]} ({a} -> {b})", value, "1/s",
              "empirical", "parkdrive", "cao2014", f"Constant transition "
              f"rate {what}.", _CODE, 0.0, 1e6)


def _hill(stem, label, base, amp, k, h, unit, falling=False):
    """The four registered numbers of one IP3 dependence (``base`` omitted
    when the printed function has none)."""
    shape = "falling" if falling else "rising"
    out = [] if base is None else [
        _p(f"pd.{stem}_base", f"{label}: IP3-independent part", base, unit,
           "empirical", "parkdrive", "cao2014", f"Constant term of {label} "
           f"(base + amp x {shape} Hill in IP3).", _CODE, 0.0, 1e4)]
    return out + [
        _p(f"pd.{stem}_amp", f"{label}: IP3-dependent amplitude", amp, unit,
           "empirical", "parkdrive", "cao2014", f"Amplitude of the {shape} "
           f"Hill term of {label}.", _CODE, 0.0, 1e4),
        _p(f"pd.{stem}_k", f"{label}: IP3 half-constant", k, "uM",
           "empirical", "parkdrive", "cao2014", f"IP3 at half of the Hill "
           f"term of {label}.", _CODE, 1e-4, 1e3),
        _p(f"pd.{stem}_h", f"{label}: IP3 Hill coefficient", h, "",
           "empirical", "parkdrive", "cao2014", f"Hill coefficient in IP3 of "
           f"{label}.", _CODE, 0.1, 10.0)]


def _const(key, name, value, unit, what, lo=0.0, hi=1e4):
    return _p(f"pd.{key}", name, value, unit, "empirical", "parkdrive",
              "cao2014", what, _CODE, lo, hi)


PD = [
    # States: drive C1 C2 C3 O6, park C4 O5; the modes meet at C2 <-> C4.
    _rate("q12", "C1", "C2", 1240.0, "within drive mode"),
    _rate("q21", "C2", "C1", 88.0, "within drive mode"),
    _rate("q23", "C2", "C3", 3.0, "within drive mode"),
    _rate("q32", "C3", "C2", 69.0, "within drive mode"),
    _rate("q26", "C2", "O6", 10500.0, "within drive mode (opening)"),
    _rate("q62", "O6", "C2", 4010.0, "within drive mode (closing)"),
    _rate("q45", "C4", "O5", 11.0, "within park mode (opening)"),
    _rate("q54", "O5", "C4", 3330.0, "within park mode (closing)"),
    # q42 = a42 + V42 m42 h42 (park -> drive)
    *_hill("v42", "V42", None, 110.0, 0.1, 2.0, "1/s"),
    *_hill("k42", "k42 (m42 half-activation)", 0.49, 0.543, 4.0, 3.0, "uM"),
    _const("n42", "Hill coefficient of m42 in Ca2+", 3.0, "",
           "Ca2+ Hill coefficient of the activating gate m42.", 0.1, 10.0),
    *_hill("kn42", "kn42 (h42 half-inhibition)", 0.41, 25.0, 6.5, 3.0, "uM"),
    _const("nn42", "Hill coefficient of h42 in Ca2+", 3.0, "",
           "Ca2+ Hill coefficient of the inhibitory gate h42.", 0.1, 10.0),
    *_hill("a42", "a42", None, 1.8, 0.58, 2.0, "1/s"),
    # q24 = a24 + V24 (1 - m24 h24) (drive -> park)
    *_hill("v24", "V24", 62.0, 220.0, 2.0, 2.0, "1/s", falling=True),
    _const("k24", "m24 half-activation", 0.35, "uM",
           "Ca2+ at half-activation of m24.", 1e-3, 1e3),
    _const("n24", "Hill coefficient of m24 in Ca2+", 3.0, "",
           "Ca2+ Hill coefficient of m24.", 0.1, 10.0),
    _const("kn24", "h24 half-inhibition", 80.0, "uM",
           "Ca2+ at half-inhibition of h24.", 1e-3, 1e4),
    _const("nn24", "Hill coefficient of h24 in Ca2+", 2.0, "",
           "Ca2+ Hill coefficient of h24.", 0.1, 10.0),
    *_hill("a24", "a24", 1.0, 20.0, 0.5, 2.0, "1/s", falling=True),
    # Relaxation rates of the gating variables (dG/dt = lambda (G_inf - G)).
    _const("lam_m24", "Relaxation rate of m24", 100.0, "1/s",
           "Rate at which m24 approaches its equilibrium."),
    _const("lam_h24", "Relaxation rate of h24", 40.0, "1/s",
           "Rate at which h24 approaches its equilibrium."),
    _const("lam_m42", "Relaxation rate of m42", 100.0, "1/s",
           "Rate at which m42 approaches its equilibrium."),
    _p("pd.lam_h42_closed", "Recovery rate of h42 (channel closed)", 0.5,
       "1/s", "empirical", "parkdrive", "cao2014", "Rate of h42 while the "
       "channel is closed: the slow recovery from Ca2+ inhibition that Cao "
       "et al. 2013 found necessary for realistic inter-puff intervals.",
       _CODE + ". Cao 2013 treats it as a free parameter (a_h42, 1 s^-1 in "
       "its Fig. 2); 0.5 is the value in the 2014 code", 0.01, 100.0),
    _p("pd.lam_h42_open", "Inhibition rate of h42 (channel open)", 20.0,
       "1/s", "empirical", "parkdrive", "cao2014", "Rate of h42 while the "
       "channel is open: fast inhibition by the Ca2+ at its own mouth.",
       _CODE + ". Cao 2013 Eq. 10 writes it as a steep Ca2+ switch at 20 uM "
       "with V_h42 = 100 s^-1; the 2014 code switches on the open state with "
       "20 s^-1, which is used here", 0.1, 1e3),
    _p("pd.ca_mouth", "Ca2+ at an open channel's mouth", 120.0, "uM",
       "empirical", "parkdrive", "cao2014", "Added to the cluster Ca2+ seen "
       "by an open receptor (its own nanodomain), after Rudiger's two-"
       "concentration puff model as used by Cao et al. 2013 (c_m = c + c_h).",
       _CODE + ": 120 uM x (store Ca2+ / 100 uM); the store is held full "
       "here, so 120 uM", 0.0, 1e3),

    # --------------------------------------------- puff comparison (method)
    _p("puff.pd_dt", "Park/drive cluster time step", 1e-4, "s", "method",
       "puff", "cao2013", "Step of the park/drive cluster simulation. The "
       "constant-rate transitions inside each mode are taken exactly over a "
       "step (matrix exponential); the Ca2+-dependent mode switches and the "
       "gating variables are updated once per step.",
       "Cao et al. 2013: 'We choose a maximum time step size of 10^-4 s'",
       1e-6, 1e-3),
    _p("puff.pd_ca_per_open", "Cluster Ca2+ per open channel (park/drive)",
       0.1, "uM", "method", "puff", "method_choice", "Mean-field cluster "
       "Ca2+ increment per open channel for the park/drive cluster; the "
       "same coupling scheme as puff.ca_per_open.",
       "Chosen by the rule used for puff.ca_per_open: the coupling (0-0.5 "
       "uM at 0.1-0.5 uM IP3) that most raises the Fano factor; it peaks at "
       "0.1-0.2 (Fano 2.7-3.6). Cao et al. 2014's own microdomain gives "
       "about 0.11 uM per open channel (k_IPR x store Ca2+ / k_diff)",
       0.0, 20.0),
    _p("puff.record_dt", "Observation bin", 1e-3, "s", "method", "puff",
       "method_choice", "Both cluster simulators record, per bin, the "
       "number open at the bin's start (for the Fano factor) and the most "
       "open at once within it (for events).",
       "Of the order of the frame time of fast puff imaging (Smith & Parker "
       "2009 record at ~2 ms); it resolves the fast park-mode flickers "
       "(0.3 ms) as brief one-channel events rather than missing them", 1e-5,
       0.1),
    _p("puff.large_fraction", "Large-event fraction of the cluster", 0.5, "",
       "method", "puff", "method_choice", "An event whose peak reaches this "
       "fraction of the cluster counts as recruiting the cluster (a puff); "
       "the recruitment share is the fraction of multi-channel events that "
       "do.",
       "Half the cluster: independent coincidences essentially never reach "
       "it (DYK uncoupled peaks at 6 of 20), while an all-or-none release "
       "should", 0.05, 1.0),
    _p("puff.scan_min", "Smallest non-zero coupling in the scan", 0.05, "uM",
       "method", "puff", "method_choice", "The coupling scan runs 0 and then "
       "geometrically from this to puff.scan_max.",
       "The step of the scan that chose puff.ca_per_open", 1e-3, 10.0),
    _p("puff.scan_max", "Largest coupling in the scan", 2.0, "uM", "method",
       "puff", "method_choice", "Top of the coupling scan.",
       "The top of the scan that chose puff.ca_per_open", 0.01, 50.0),
    _p("puff.scan_points", "Non-zero couplings in the scan", 7.0, "",
       "method", "puff", "method_choice", "Number of geometric steps from "
       "puff.scan_min to puff.scan_max.",
       "About two points per factor of two over 0.05-2 uM", 2.0, 50.0),
]
