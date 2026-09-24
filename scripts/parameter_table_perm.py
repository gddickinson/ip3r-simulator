"""Unitary conductance: the drift-diffusion transport constants, the fixed-
charge map and the published conductances it is compared with. Imported into
``parameter_table.P``.

The measured values were read from the papers' own text on 2026-09-24
(PMC2995152, PMC2217211). Both recordings are nuclear patch clamp of rat
InsP3R-3 in symmetric 140 mM KCl with no Mg2+, at room temperature.
"""

from param_entry import entry as _p

_VAIS = ("Vais et al. 2010 (PMC2995152), Results: 'symmetric 140-mM KCl "
         "solutions ([K+]f = [Cl-]f ~ 104 mM) in the absence of Mg2+'; "
         "'current traces were acquired at room temperature'")

PERM = [
    # ------------------------------------------------ transport constants
    _p("permeation.temperature", "Temperature", 296.0, "K", "convention",
       "permeation", "vais2010", "Absolute temperature in the thermal "
       "voltage and the Nernst-Einstein conductivity.",
       _VAIS + "; 23 C taken for room temperature", 270.0, 330.0),
    _p("permeation.bath_concentration", "Bath KCl concentration", 0.14, "M",
       "convention", "permeation", "vais2010", "Symmetric bath salt. The "
       "nominal concentration, not the ~104 mM free (activity) value the "
       "paper also gives; the model has no activity coefficients.",
       _VAIS, 0.001, 3.0),
    _p("permeation.test_voltage", "Test voltage", 0.02, "V", "method",
       "permeation", "method_choice", "Voltage at which the chord "
       "conductance is evaluated.", "Inside Vais 2010's -25 to +15 mV ramp, "
       "where their I-V is linear, so chord equals slope conductance",
       0.001, 0.2),
    _p("permeation.diffusion_potassium", "K+ bulk diffusivity", 1.96e-9,
       "m^2/s", "physical", "permeation", "convention", "Bulk diffusion "
       "coefficient of K+.", "K+ in water at 25 C; the standard tabulated "
       "value (as PIEZO1)", 1e-11, 1e-8),
    _p("permeation.diffusion_chloride", "Cl- bulk diffusivity", 2.03e-9,
       "m^2/s", "physical", "permeation", "convention", "Bulk diffusion "
       "coefficient of Cl-.", "Cl- in water at 25 C; the standard tabulated "
       "value (as PIEZO1)", 1e-11, 1e-8),
    _p("permeation.radius_potassium", "K+ crystal radius", 1.38, "A",
       "physical", "permeation", "convention", "Subtracted from the free "
       "radius to give the area a K+ centre can occupy.", "Shannon radius. "
       "Crystal rather than hydrated, because an ion can shed water at a "
       "constriction - itself a modelling choice (as PIEZO1)", 0.3, 4.0),
    _p("permeation.radius_chloride", "Cl- crystal radius", 1.81, "A",
       "physical", "permeation", "convention", "Subtracted from the free "
       "radius for Cl-.", "Shannon radius", 0.3, 4.0),
    _p("permeation.diffusion_scale", "In-pore diffusivity fraction", 0.5, "",
       "method", "permeation", "unverified", "In-pore diffusivity as a "
       "fraction of bulk, for both species.", "Not measured for any IP3R. "
       "PIEZO1's value, kept so the two projects are comparable; with the "
       "ion radius it is what the answer is most sensitive to, so it is "
       "swept (permeation.sweep_*), never tuned", 0.01, 1.0),
    _p("permeation.permittivity_pore", "In-pore relative permittivity", 40.0,
       "", "method", "permeation", "unverified", "Enters only the reported "
       "Debye length.", "Nanopore water is less polarisable than bulk (80); "
       "40 is a common compromise (as PIEZO1)", 2.0, 80.0),
    _p("permeation.sweep_scale_low", "Sweep: lowest diffusivity fraction",
       0.25, "", "method", "permeation", "method_choice", "Low end of the "
       "in-pore diffusivity sweep.", "PIEZO1's plausible range, 0.25-1.0 of "
       "bulk", 0.01, 1.0),
    _p("permeation.sweep_scale_high", "Sweep: highest diffusivity fraction",
       1.0, "", "method", "permeation", "method_choice", "High end of the "
       "in-pore diffusivity sweep (bulk water).", "PIEZO1's plausible range",
       0.01, 1.0),
    _p("permeation.sweep_radius_low", "Sweep: smallest ion radius", 1.0, "A",
       "method", "permeation", "method_choice", "Low end of the ion-radius "
       "sweep.", "PIEZO1's plausible range, 1.0-2.0 A", 0.3, 4.0),
    _p("permeation.sweep_radius_high", "Sweep: largest ion radius", 2.0, "A",
       "method", "permeation", "method_choice", "High end of the ion-radius "
       "sweep.", "PIEZO1's plausible range", 0.3, 4.0),
    # ------------------------------------------------ measured conductances
    _p("permeation.published_itpr3_dt40", "ITPR3 conductance (DT40)", 545.0,
       "pS", "empirical", "permeation", "vais2010", "Rat InsP3R-3 in the "
       "outer nuclear membrane of DT40-KO cells, symmetric 140 mM KCl.",
       "Vais 2010: 'was 545 +/- 7 pS (n = 16)'", 1.0, 2000.0),
    _p("permeation.published_itpr3_oocyte", "ITPR3 conductance (oocyte)",
       358.0, "pS", "empirical", "permeation", "mak2000", "Rat InsP3R-3 in "
       "the Xenopus oocyte outer nuclear membrane, symmetric KCl, 0 Mg2+.",
       "Mak 2000 (PMC2217211): 'the InsP3R-3 I-V relation became linear "
       "with a conductance of 358 +/- 8 pS'. Vais 2010 cites this paper as "
       "370 +/- 8 pS; the primary text is used", 1.0, 2000.0),
    # ------------------------------------------------ fixed charge
    _p("pore_charge.lining_margin", "Charge lining margin", 3.0, "A",
       "method", "pore_charge", "method_choice", "A charge centre lines the "
       "pore when it is at most this far outside the profile's atom-centre "
       "radius at its own height.", "About one water diameter (2.8 A): a "
       "charge one hydration layer behind the lumen wall still acts on the "
       "ions in it. Swept in SESSION_LOG, not tuned", 0.0, 15.0),
    _p("pore_charge.stub_reach", "Stubbed side-chain reach", 7.3, "A",
       "convention", "pore_charge", "convention", "A stubbed ionisable "
       "residue is counted as unplaced when its C-alpha is within this plus "
       "the lining margin of the lumen.", "C-alpha to guanidinium centre of "
       "an extended Arg, the longest ionisable side chain (PIEZO1's "
       "pore_charge.reach_arg)", 1.0, 12.0),
    _p("pore_charge.smoothing", "Axial charge smoothing", 3.0, "A", "method",
       "pore_charge", "method_choice", "Gaussian width each charge is spread "
       "over along the axis.", "About a side chain's positional uncertainty "
       "at 3 A resolution; the total charge is conserved whatever it is "
       "(as PIEZO1)", 0.5, 15.0),
    _p("pore_charge.salt_bridge_cutoff", "Salt-bridge N-O distance", 4.0,
       "A", "convention", "pore_charge", "barlow1983", "An acidic and a "
       "basic side chain are an ion pair (and cancel) when a charged oxygen "
       "and a charged nitrogen are at most this far apart.", "Barlow & "
       "Thornton 1983, abstract: ion pairs 'less than or equal to 4 A "
       "between charged groups', from the like- vs opposite-charge distance "
       "distributions in 38 proteins. On 8TKF the answer depends on it: "
       "K2482-D2400 and D2518-R2524' sit at 4.2-4.4 A, so it is swept in "
       "SESSION_LOG, not tuned", 2.0, 8.0),
    _p("pore_charge.max_concentration", "Counterion packing ceiling", 10.0,
       "M", "method", "pore_charge", "method_choice", "In-pore concentration "
       "above which a result is flagged; nothing is clipped to it.",
       "Close-packed hydrated K+ is about 8 M, so 10 M is where a continuum "
       "of point ions has certainly stopped describing a solution (as "
       "PIEZO1)", 0.1, 100.0),
]
