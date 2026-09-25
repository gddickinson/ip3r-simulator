"""The puff microdomain (Round 7.2): Cao et al. 2014's Ca2+ pools around a
park/drive cluster, with fluo-4 bound in the microdomain as in Cao et al.
2013's Eq. 12. Imported into ``parameter_table.P``.

The pool constants were read from Cao et al. 2014 Text S1 (the stochastic
6-state script) on 2026-09-25. Cao 2013's own point-domain constants are in
its Table S2, which is behind PMC's browser challenge. Fluo-4 comes from
Shuai, Rose & Parker 2006, Table 1. That table gives K_d = 2 uM, the value
Cao 2013's text states for fluo-4.
"""

from param_entry import entry as _p

_CODE = ("Cao et al. 2014 Text S1 (6-state stochastic model code), read "
         "2026-09-25")
_FLUO = "Shuai, Rose & Parker 2006, Table 1 (fluo-4 dextran)"


def _pool(key, name, value, unit, what, note=_CODE, lo=0.0, hi=1e4):
    return _p(f"domain.{key}", name, value, unit, "empirical", "domain",
              "cao2014", what, note, lo, hi)


def _method(key, name, value, unit, what, note, lo, hi, citation="method_choice"):
    return _p(f"domain.{key}", name, value, unit, "method", "domain",
              citation, what, note, lo, hi)


DOMAIN = [
    # ------------------------------------------------------------- pools
    _pool("gamma1", "Cytosol : microdomain volume ratio", 100.0, "",
          "Fluxes into and out of the microdomain are multiplied by this, "
          "so the microdomain relaxes 100x faster than the cytosol.",
          lo=1.0),
    _pool("gamma2", "Cytosol : ER volume ratio", 10.0, "",
          "Store Ca2+ is cs = gamma2 (ct - cb/gamma1 - c).", lo=1.0),
    _pool("k_ipr", "Release rate constant per open receptor", 0.0025, "1/s",
          "Flux through one open receptor, k_ipr (cs - cb), in cytosolic "
          "units.", _CODE + ": kipr = 0.05/Num_IPR with Num_IPR = 20. "
          "Registered per receptor, the constant single-channel flux Cao "
          "2013 assumes, so a larger cluster releases more", 0.0, 10.0),
    _pool("k_diff", "Microdomain -> cytosol exchange", 10.0, "1/s",
          "Jdiff = k_diff (cb - c); with gamma1 the microdomain's own "
          "relaxation is gamma1 x k_diff = 1000/s.", lo=0.01),
    _pool("k_leak", "ER leak", 0.0032, "1/s", "Jleak = k_leak (cs - c)."),
    _pool("v_serca", "SERCA maximum", 10.0, "uM/s",
          "Jserca = V c^n / (c^n + K^n)."),
    _pool("k_serca", "SERCA half-activation", 0.26, "uM", "SERCA K."),
    _pool("n_serca", "SERCA Hill coefficient", 1.75, "", "SERCA n.",
          lo=0.5, hi=6.0),
    _pool("v_pm", "Plasma-membrane pump maximum", 0.8, "uM/s",
          "Jpm = V c^2 / (K^2 + c^2)."),
    _pool("k_pm", "Plasma-membrane pump half-activation", 0.5, "uM",
          "PMCA K."),
    _pool("j_leak_in", "Plasma-membrane leak in", 0.03115, "uM/s",
          "Constant Ca2+ entry."),
    _pool("v_rocc", "Receptor-operated entry per uM IP3", 0.2, "1/s",
          "Jrocc = v_rocc x p (uM/s with p in uM)."),
    _pool("v_socc", "Store-operated entry maximum", 1.6, "uM/s",
          "Jsocc = V K^4 / (cs^4 + K^4): entry rises as the store empties."),
    _pool("k_socc", "Store-operated entry half-point", 100.0, "uM",
          "SOCC K, in store Ca2+."),
    _pool("mouth_per_store", "Mouth Ca2+ per store Ca2+", 1.2, "",
          "An open receptor sees mouth_per_store x cs at its own mouth "
          "(cm = 120 (cs/100) in the code).",
          _CODE + ". At the model's resting store (about 490 uM) this is "
          "about 590 uM. pd.ca_mouth = 120 uM is the same rule at a store "
          "of 100 uM", 0.0, 10.0),
    _pool("ct0", "Starting total Ca2+", 45.0, "uM",
          "The code's initial total (cytosolic units). With c = cb = "
          "puff.ca_rest it fixes the starting store, gamma2 (ct0 - c0/gamma1 "
          "- c0) = 449 uM, which the clamped store keeps.",
          _CODE + ": ct0 = 45, 'total calcium concentration'. It is not the "
          "model's steady state (that store is 630-1340 uM at 0.05-0.5 uM "
          "IP3, reached over minutes)", 0.0, 1e3),
    # ------------------------------------------------------------- fluo-4
    _p("domain.fluo_total", "Fluo-4 total", 25.0, "uM", "empirical",
       "domain", "shuai2006", "Total indicator in the microdomain "
       "(Cao 2013 Eq. 12's B_fluo4). Set to 0 for the dye-free control.",
       _FLUO + "; the oocyte loading of Rose et al.", 0.0, 1e3),
    _p("domain.fluo_kon", "Fluo-4 on-rate", 150.0, "1/(uM s)", "physical",
       "domain", "shuai2006", "Ca2+ binding to fluo-4.", _FLUO, 1.0, 1e4),
    _p("domain.fluo_koff", "Fluo-4 off-rate", 300.0, "1/s", "physical",
       "domain", "shuai2006", "Ca2+ unbinding; K_d = koff/kon = 2 uM.",
       _FLUO + "; Cao 2013's text: 'the dissociation constant of fluo-4 "
       "of 2 uM'", 1.0, 1e5),
    # ------------------------------------------------------------- method
    _method("puff_threshold", "Puff threshold in mean blips", 1.875, "",
            "A fluorescence event is a puff when its peak dF/F0 exceeds "
            "this multiple of the mean blip (one-channel event) amplitude.",
            "Cao 2013 set dF/F0 > 3 with a mean blip of 1.6; their absolute "
            "3 depends on their Table S2, so the ratio 3/1.6 is carried",
            1.0, 20.0, citation="cao2013"),
    _method("detect_fraction", "Event detection level", 0.5, "",
            "A fluorescence event is a run of bins whose dF/F0 exceeds this "
            "fraction of the steady dF/F0 one open receptor gives (the "
            "blip). Its peak open count says blip (one) or more.",
            "Half a blip: a park-mode flicker (0.3 ms) stays below it, and "
            "a drive-mode opening (tens of ms) crosses it", 0.01, 10.0),
    _method("ipi_duration", "Simulated time per IPI record", 300.0, "s",
            "Length of each run behind an inter-puff-interval distribution.",
            "Cao 2013's inter-puff intervals are 2-10 s, so 300 s gives "
            "30-150 intervals per run", 10.0, 1e5),
    _method("ah42_min", "Slowest h42 recovery in the scan", 0.1, "1/s",
            "The a_h42 scan (pd.lam_h42_closed) runs geometrically from "
            "this to domain.ah42_max.", "Cao 2013 Fig. 4's range",
            1e-3, 100.0, citation="cao2013"),
    _method("ah42_max", "Fastest h42 recovery in the scan", 5.0, "1/s",
            "Top of the a_h42 scan.", "Cao 2013 Fig. 4's range", 1e-3, 100.0,
            citation="cao2013"),
    _method("n_min", "Smallest cluster in the N scan", 3.0, "",
            "The amplitude-vs-N scan runs from this to domain.n_max.",
            "Cao 2013 Figs. 6 and 8 (3-25 receptors)", 1.0, 100.0,
            citation="cao2013"),
    _method("n_max", "Largest cluster in the N scan", 25.0, "",
            "Top of the amplitude-vs-N scan.", "Cao 2013 Figs. 6 and 8",
            1.0, 200.0, citation="cao2013"),
    _method("store_scale_max", "Largest release-rate multiple", 20.0, "",
            "The store scan multiplies domain.k_ipr by 1 up to this, which "
            "raises the coupling per open receptor from ~0.1 uM toward the "
            "~0.5-2 uM where the mean-field cluster stays open.",
            "Reaches the mean-field couplings of the Round 4 question",
            1.0, 1e3),
    _method("ah42_points", "Points in the a_h42 scan", 5.0, "",
            "Geometric steps from domain.ah42_min to domain.ah42_max.",
            "About two per factor of five over Cao 2013's range", 2.0, 50.0),
    _method("n_points", "Points in the N scan", 7.0, "",
            "Geometric steps (rounded, distinct) from domain.n_min to "
            "domain.n_max.", "Enough to see a bend near Cao's N = 12",
            2.0, 50.0),
    _method("store_points", "Points in the store scan", 5.0, "",
            "Geometric steps of the release-rate multiple from 1 to "
            "domain.store_scale_max.", "About two per factor of five",
            2.0, 50.0),
    _method("store_duration", "Simulated time per store-scan run", 60.0,
            "s", "Length of each clamped / free / mean-field run in the "
            "store scan.", "Long enough for the store to fall by half at "
            "the sustained open fraction of the Round 4 question", 1.0, 1e4),
    _method("scan_seed", "Seed of every scan run", 0.0, "",
            "One seed for every point, so points differ only in the "
            "scanned constant.", "Reproducibility", 0.0, 1e9),
    _method("blip_scale", "Release multiple matching Cao's mean blip", 2.5,
            "", "The k_ipr multiple the N scan reads at by default: there the "
            "measured mean blip (one-channel event) is dF/F0 1.33-1.66, "
            "against Cao 2013's 1.6. At 1x (Cao 2014's release) it is "
            "0.58-0.72.", "Found by N scans at 1, 2.5 and 4x on 2026-09-25 "
            "(Round 7.2); Cao 2013's own release is in its Table S2, which "
            "was not reachable", 0.01, 100.0, citation="cao2013"),
]
