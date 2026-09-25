# ROADMAP_RYR — Round 6, the ryanodine receptors

Split from `ROADMAP.md` (which passed 490 lines) the way `docs/SCIENCE_RYR.md`
was split from `docs/SCIENCE.md`. Same conventions: `[ ]` planned, `[x]` done
with what it measured.


- [x] 6.1 RyR1 structures through the whole structural pipeline. There is
  a curated resource (`scripts/curate_ryr.py`, `make ryr`) with a sequence,
  Pfam domains and six deposits selected by five stated rules plus a
  same-paper morph partner. `RYR1` is a numbering, not a publication
  paralog. The measurements: all in P11716 numbering; the shut states gate
  at I4937; only 9HEO opens (5.05 Å). 9R8O → 9HEO morphs, with ANM overlap
  0.17 (null 0.03). 9HEO gives 136 pS neutral and 180 pS charged against
  801 pS measured, so the continuum is short on both receptors. The charge
  mutants (Xu 2006): direction right in 4/4 lining residues, E4955Q null
  right; D4899Q 0.20× measured vs 0.90× modelled. D4899 is bridged, so
  pairing predicts 1.00×, refuting the Round 4 paired reading. The GUI
  follows the loaded deposit's family: state panel, conductance, a
  mutant exhibit, and a primed → open preset.
  Emergent:
  - [ ] The filter charge is 4× too weak in the continuum. Try a filter-
    local treatment (charge not spread over 3 Å; a radial rather than
    cross-section-averaged Donnan) and hold it to the five mutants, not to
    the wild-type number.
  - [ ] RyR2 (cardiac), whose closed/open deposits exist, by the same rules.
- [x] 6.2 RyR1 gating and sparks (`physics/ryr_gating.py`,
  `physics/sparks.py`; Gating tab third model, Puffs tab third receptor;
  `ryr-gating`, `sparks` CLI). The sources were found by a readable-
  constants search: Murayama 2015's measured bell and Stern 1997's unfitted
  two-gate scheme. Stern's printed k_i is a sign typo, corrected from the
  text's 10 µM Kd; the printed value would remove inhibition (tested). The
  scheme's half-activation matches (3.9 vs 4.4 µM), but its inactivation
  is 6.7× too sensitive (48 vs 320 µM). Coupling is derived from Stern's
  current, diffusivity and spacing (8.25 µM, not tuned). Sparks: uncoupled
  only blips; coupled 1.5/s reaching 25–30 of 30, Fano 4; switched on at
  0.05–0.1× the derived coupling. The step was measured (2.5e-5 s). Sparks
  last ~120 ms against 6.3 ms measured (frog): the mean-field cluster sits
  at a self-sustaining point (~5.5 open) until it closes by chance.
  Emergent:
  - [x] A spatial Ca²⁺ field for the cluster: done as Round 6.3.
  - [x] Stern's inactivation against Murayama's bell: done as Round 6.4.
    Fitted, sparks do not terminate.
  - [ ] Murayama's S1 Table has 11 MH/CCD mutants (Amax, KA, KI at 25 and
    37 °C). A mutant bell could be drawn beside the variant spheres.

- [x] 6.3 Sparks in the junctional cleft (`physics/cleft.py`,
  `physics/sparks_cleft.py`; Puffs tab fourth receptor; `sparks --cleft`).
  Stern 1997's own geometry (60 × 15 nm cleft, V/C chessboard, 30 nm
  source disc, Eq. 13 leaky edges) is solved by finite volumes. The solve
  conserves Ca²⁺ exactly and is grid-converged. Per pA it gives 171/65/56 µM
  at the source, across the row and 30 nm along it, against ~185/75/73 read
  from their Fig. 9, so it falls off somewhat faster. Each channel sees its
  own Ca²⁺: 51 µM from itself and 11.4 µM from a diagonal neighbour. With
  all the others open that is 36 µM, against 239 µM in the mean-field
  cluster. The array is simulated exactly by Gillespie's method, with no
  step. Sparks last a median 21 ms (IQR 17–26) against ~130 ms mean-field,
  and they end with 16 of 30 channels inactivated (0 at the start): local
  inactivation, as Stern described. The duration holds at 16–22 ms over
  0.75–2× coupling and 15–45 nm source (own release 40–74 µM). At 0.5×
  there are no sparks. It is still 3× the 6.3 ms frog release.
  `docs/SCIENCE.md` passed 500 lines, so the RyR1 sections moved to
  `docs/SCIENCE_RYR.md`.
  Emergent:
  - [ ] Draw the cleft: the field and each channel's state at a spark's
    peak, in the Puffs tab (`cleft.field` exists).
  - [ ] Mobile buffers (their fura-2 runs, 3 mM suppressing release) need
    time-dependent diffusion in the cleft; the steady field cannot hold them.
  - [ ] The solve decays faster than Fig. 9 (15 vs ~25 µM/pA at 60 nm).
    Their edge coefficient or source may differ from the reading here; a
    finer digitisation of Fig. 9 would say which.
  - [ ] IP3R puffs with a spatial field: the cluster has no cleft, so the
    field is open cytosol (point sources with buffers). The geometry of an
    IP3R cluster would have to be sourced first.

- [x] 6.4 What ends a cleft spark (`physics/spark_termination.py`,
  `ryr_gating.fit_to_bell`; `spark-termination` and `sparks --cleft --fit`
  in the CLI; a fifth Puffs receptor; the fitted bell dashed in Gating).
  Stern's Ka and Ki were fitted to Murayama's two flanks, with off rates
  moved, on rates kept, and flanks exact. That gives Ka 4.9 µM and Ki
  249 µM at 25 °C, and Ka 16.9 µM and Ki 358 µM at 37 °C (S1 Table, now
  registered), against 7.1 and 10. In the cleft (4 seeds × 10 s), Stern's
  sparks last 19 ms and all end. Fitted to 25 °C, the array never shuts
  after the first spark (4/4 unended, 70 % open). Fitted to 37 °C, no
  spark starts. Ki scan: 17–19 ms up to 10 µM, 89 ms at 64 µM, 1.4 s at
  118 µM, unended from 217 µM. Rate scan: at Ki 10 the shortest is
  13–15 ms, reached only as sparks stop starting (none at 30×). At the
  fitted constants no rate from 1/30× to 30× ends a spark. The 3× gap was
  not caused by the inactivation. Stern's too-sensitive inactivation is
  what let sparks end at all.
  Emergent:
  - [x] Mg²⁺: done as Round 6.5.
  - [x] Luminal Ca²⁺ depletion as a terminator: done as Round 6.7. It is
    not the terminator.
  - [x] The V channels: done as Round 6.6.
  - [x] A one-Ca²⁺ inactivation gate has Hill slope 1 against the
    bell's 1.5: the two-site gate of Round 6.8 fits it.

- [x] 6.5 Mg²⁺ and the cleft spark (`physics/spark_mg.py`; `mg`,
  `k_mg_a`, `mg_i` in `SternParams`; `spark-mg` CLI). No open RyR1 bell
  with Mg²⁺ could be read (Meissner 1997 JBC is behind a browser
  challenge). The two actions come from readable sources instead:
  competition at the activation site, K_Mg,A 54 µM (Laver 2004 Table I,
  with ATP), and equal Ca²⁺/Mg²⁺ inhibition at the I1 site (Laver 1997).
  Fibre free Mg²⁺ is 1 mM. Nothing is fitted. Under Mg²⁺ no spark starts
  by itself, so sparks are triggered (available channels opened at
  t = 0). Fitted to the 25 °C bell, the array never shuts without Mg²⁺.
  With 1 mM at the activation site alone it shuts every time, in 17.5 ms
  (54 µM reading) or 49.5 ms (Laver's selectivity on the fitted Ka,
  521 µM), with **no channel inactivated**. The mechanism is induction
  decay, not inactivation. Adding the I1 site holds 24 of 30 channels
  shut at rest; then 6 open and the spark lasts 4–6 ms. Scan: the array
  never shuts up to 25 µM Mg²⁺, and does so in 131 ms at 63 µM and 4 ms
  at 1 mM. Stern's published Ki (10 µM) cannot take Mg²⁺ at the I1 site:
  all 30 are inactivated at rest.
  Emergent:
  - [x] K_Mg,A in Murayama's condition (`physics/mg_competition.py`):
    Meissner 1997 Table IV, same assay, Mg²⁺ Ki 18 µM (+AMP). Na⁺ holds
    the same site in Murayama's 0.17 M NaCl (Scheme 2, Ki 24 mM), so
    1 mM Mg²⁺ shifts Ka 2.3×, not 11×: K_Mg,A ≈ 769 µM (345–1,600 over
    the error bars). The transfer is checked on the paper's own Table II
    (three salts within 1.6×). Activation site alone: triggered spark
    435 ms (32 ms at best) vs 6.3 ms measured; both sites 6.0 ms at every
    reading. Caveat: Murayama's Ka is 7× what the same law predicts for
    their salt.
  - [ ] Mg²⁺ binding kinetics: rapid equilibrium is assumed; Laver 2004
    say Mg²⁺ may not re-equilibrate within an opening.
  - [x] The V-channel trigger: done as Round 6.6. Neither row gives the
    release waveform.
  - [x] Mg²⁺ in the GUI. Puffs: free Mg²⁺ and the K_Mg,A reading on
    every RyR1 receptor, plus "Triggered sparks vs Mg²⁺" on the cleft
    receptors (never shuts ≤ 25 µM; 131 → 4 ms from 63 µM to 1 mM, with
    6 of 30 opened and 24 inactivated before the trigger at 1 mM). Gating:
    the fitted bell under 1 mM Mg²⁺ keeps 17 % of its Mg²⁺-free peak, with
    half-activation at 77 µM (54 µM reading) or 21 % and 13 µM (521 µM).
    The smoke test found that a silent cluster crashed the event-size
    plot (log axis with no data); fixed, and now a step.

- [x] 6.6 The V channels: the couplon under voltage clamp
  (`physics/allosteric_v.py`, `physics/couplon.py`, `physics/ec_release.py`;
  `ec` CLI). Rios 1993's ten-state allosteric model is used at fiber 827,
  the row that gives Stern's V plateaus (0.457/0.052/0.003). Stern's text
  says the rates were doubled. Their couplon figure (12) follows ×2, while
  their stand-alone check (Fig. 11) follows ×1, so it predates the factor.
  V trajectories are exact, then fed into the C Gillespie as scheduled
  events. The Monte Carlo equals the master equation. Calibration: Stern's
  constants reproduce his Fig. 12 (C peak 0.64/0.38/0.11 vs
  0.64/0.39/0.12; plateau 0.12/0.09/0.04 vs 0.105/0.078/0.039), with
  release stopping at repolarisation. Finding: fitted to Murayama's bell,
  the C array gives flux peak/plateau 1.0–1.1 (Stern 2.7) under every
  Mg²⁺ arrangement and reading. No Mg²⁺ or activation-site Mg²⁺: release
  continues after repolarisation (C Po 0.72, 0.39). Both sites: control
  is kept but the C channels barely open (0.12 at 0 mV, 0.014 at
  −30 mV). The time-dependent terminator is missing.
  Emergent:
  - [x] Luminal depletion: done as Round 6.7.
  - [ ] The couplon in the GUI: V and C open probability and flux per
    step, and the configuration selector (a Puffs sub-tab or its own).
  - [ ] Stern's other protocols as calibrations: the peak/plateau bell
    against voltage (Fig. 2 B), quantal two-pulse release (Figs. 4–6),
    couplon length (Fig. 13).
  - [ ] Rios 1993's Appendix model (eight charges per channel), which
    saturates Po before charge; it may change how many C channels a V
    channel recruits at low voltage.

- [x] 6.7 SR depletion as the terminator (`physics/lumen.py`;
  `ec --depletion [--pool-scan] [--large]`). The model is Stern's own
  Fig. 20: one well-mixed pool of 2 mM, 4.8 couplons/µm³, every unitary
  current proportional to the content, and uptake first order with
  τ 0.43 s (digitised from Fig. 20 C). The shared path is made consistent
  with the ensemble by damped fixed-point iteration. Calibration, Stern's
  constants on his 28-channel couplon: 0.54 mM left at 100 ms [0.655] at
  0 mV and 1.19 [1.34] at −30 mV; corrected plateau 1.11 pA [1.2]. It
  over-depletes slightly, the generous direction. Finding: at 2 mM the
  fitted scheme regains control (C after repolarisation 0.001, against
  0.72) and a corrected peak (2.9), but it releases 87–89 % of the store.
  Stern's constants release Rios's measured 50–60 % only at twice the pool
  (57 %). At that pool the fitted scheme releases 73–83 %, with a C Po
  peak/plateau of 1.5–2.7 against Stern's 4.6, which does not depend on
  the pool. **Depletion makes a peak only by emptying more store than a
  fibre loses.** This agrees with Launikonis 2006: skraps are ≤ 7.4 %.
  Emergent:
  - [x] Two-site Ca²⁺ inactivation: done as Round 6.8. It makes things
    worse.
  - [ ] Inactivation kinetics: Stern's I-gate rates are not measured by a
    steady-state bell. A sourced time course (e.g. RyR1 inactivation after
    a Ca²⁺ step in bilayers) would constrain them.
  - [ ] The 60-channel couplon through the depletion panel (not run; it
    doubles the channel density at a fixed couplon density).
  - [ ] A Ca²⁺-driven pump (Brum 1988 removal model, González & Ríos 1993
    Fig. 2 constants, in Rios 1993's Fig. 1 legend) in place of first-order
    uptake, with the global cytosolic Ca²⁺ that drives it.

- [x] 6.8 Two-site Ca²⁺ inactivation (`physics/ryr_two_site.py`;
  `ec --two-site`, `spark-mg --two-site`). Two Ca²⁺ bind in sequence (Stern's
  on rate for each) and inactivate only together, giving six states. The
  simulators now read each scheme's own tables (`dest`, `open_mask`,
  `inact_mask`, `trigger_map`), and Stern's rows reproduce Round 6.6
  unchanged. The scheme was fitted to three numbers of Murayama's 25 °C bell
  on one ruler (both half-peak points and the log-log slope at half
  inhibition, −0.84 against the one-site −0.56). Result: Ka 4.57, K1 512,
  K2 115 µM; gate Hill slope 1.63. Murayama *fixed* nI 1.5, so the slope is
  a convention of their fit. The finding: with the half point fixed, a
  steeper gate inactivates less below it. At 30/50 µM (cleft Ca²⁺) it
  inactivates 1.4/3.7 % against 11/17 % for the one-site gate. In the
  couplon, with no Mg²⁺, the C channels stay 0.885 open after
  repolarisation [one-site 0.72]; with activation-site Mg²⁺, 0.73 [0.39];
  with both sites, C Po is 0.09 [0.12] at 0 mV with no peak. At −50 mV
  there are only blips. Triggered sparks without Mg²⁺ never shut (0/20).
  **No steady-state inactivation consistent with the bell terminates the
  couplon.** Stern's works because its Ki (10 µM) lies inside the cleft's
  range.
  Emergent:
  - [x] Why the fibre's effective Ki would be ~30× below the bell's:
    answered in Round 6.9, and not by the bell's conditions. Crediting
    part of the descending limb to a use-dependent gate puts Ki at
    21-59 µM.
  - [x] A flux-driven (non-equilibrium) inactivation, which a binding
    bell cannot show, as the remaining mechanism, if a source exists: done
    as Round 6.9 (Laver & Lamb 1998).

- [x] 6.9 Use-dependent (flux-driven) inactivation (`physics/ryr_use.py`;
  `spark-termination --scan use [--bell] | --scan ratio`, `ec --use`).
  Laver & Lamb 1998's inactivation, whose rate is set by how much the
  channel is *open* and not by any ligand: a third gate entered only from
  the conducting state, eight states, a cycle with no detailed balance.
  **The two inactivations partly cancel**, because shutting the channel by
  the Ca²⁺ route protects it from the use gate; so the gate does not hide
  in a bell but widens it (1.86 -> 2.90 decades at the measured τ), and
  the Ca²⁺ gate must be refitted *beside* it or the bell is counted twice.
  Doing that: the fastest use gate the bell permits is τ ≈ 1 s, the fast
  end of the measured 1-3 s (an independent agreement), and the refitted
  Ki is 21-59 µM rather than 249 — within 2-6× of Stern's 10 µM, which is
  where the "30×" came from. The composite scheme is the first in Rounds
  6.3-6.9 that reproduces the bell *and* terminates cleft sparks (9-23 ms
  against a measured 6.3) *and* keeps couplon control (C after
  repolarisation 0.0002-0.0094 against 0.7265) *and* shows a release peak
  (flux peak/plateau 1.3-1.8 against Stern's 2.71, measured ~4.6). But it
  holds only for a recovery τ near 10 s, which **could not be sourced**:
  faster and Ki climbs back to 249 and nothing terminates, slower and no
  spark starts. Round 6.9 therefore yields a mechanism with a falsifiable
  requirement, not a settled answer.
  Emergent:
  - [x] Laver & Lamb 1998's recovery rate: read from the full text in
    Round 6.10. It is not reported; the round rests on a ratio instead.
  - [x] Per-channel heterogeneity: only half to two-thirds of RyRs
    inactivate, stably, so a cluster keeps channels that never do. Done as
    Round 6.11: it does make termination harder, through the bell.
  - [ ] A time-dependent cleft field, for Laver et al. 2013's induction
    decay. With neither inactivation gate this model's sparks never end
    (open fraction 0.942), because `sparks_cleft` holds the field at its
    steady state between events, and induction decay lives in the time
    course of the gradients.

- [x] 6.10 Laver & Lamb 1998 read in full (`ryr_use`: `use_rate_at`,
  `residual_bound`, `ratio_values`; `spark-termination --scan ratio`). The
  paper reports **no recovery rate at a fixed potential** (recovery only on
  voltage reversal, in ms). Its 1-3 s is at +40 mV; Fig. 4 (log10, checked
  against the abstract and z delta) puts the entry rate at 0 mV at
  0.056 s^-1, now `ryr.k_use_on`. Because both use-gate rates are seconds,
  the bell, the refitted Ca2+ gate and termination depend only on
  **rho = k_use-/k_use** (tested), now `ryr.use_recovery_ratio`. Sparks of
  the fibre's order need rho 0.1-0.4 (none fire below 0.1; >0.6 they stop
  ending); the couplon at rho 0.2 keeps control at every speed. Fig. 8's
  residuals bound rho <= 0.03-0.52 at +40 mV only, a span containing the
  band; at 0 mV rho is unmeasured. **Corrected:** Round 6.9's "no Ca2+
  gate fits beside a gate faster than tau ~0.3 s, so the bell permits
  tau >= 1 s, independently matching 1-3 s" was a stalled single-start
  solver (a fit exists down to rho 0.002); `fit_with_use` is now
  multi-start and the claim is withdrawn.
  Emergent:
  - [ ] rho at 0 mV: the steady residual activity of RyR1 held open at
    0 mV (the model needs ~10-30 %). Sitsapesan 1995 and Laver & Curtis
    1996 read in full (Round 6.11): neither records at 0 mV. Needs a
    bilayer at 0 mV or a longer +40 mV record than Fig. 8's 5 s.

- [x] 6.11 A mixed cluster (`physics/ryr_mixed.py`; `spark-termination
  --scan fraction`). A registered 0.8 of channels carry the use gate
  (Laver & Lamb 1998: 12 of 15 skeletal RyRs; Laver & Curtis 1996 70 %;
  Ma 1995 50-70 %), drawn per run; one Ca2+ gate shared, **fitted to the
  population's mean bell** because [3H]ryanodine binding is a population
  measurement. At rho 0.2: f 1 -> Ki 40 uM, sparks 14 ms; **f 0.8 -> Ki
  104 uM, sparks 228 ms**; f <= 0.6 never end; f >= 0.9 needed for sparks
  under ~30 ms. The couplon at 0.8 loses control (C after repolarisation
  0.33 vs 0.0018). **Taken apart, the non-inactivating channels alone
  cost nothing (14 ms); the weaker Ca2+ gate the population bell forces
  costs it all (151 ms).** The use gate reconciles bell and termination
  only if nearly every channel carries it, against a measured 0.5-0.8.
  Laver & Curtis 1996 is all +40 mV. Sitsapesan 1995, read in full after
  the round, is -40/+40 mV only; its -40 mV record (no decline in 5 s at
  Po 0.86) admits the registered entry rate (3.7 % predicted) and rejects
  Round 6.9's 0.5 s^-1 (29 %, beyond 2 SD; `decline_bound`), bounding the
  0 mV rate at about 0.44 s^-1. Rho at 0 mV remains unmeasured.
  Emergent:
  - [x] Separate Ca2+ gates for the two populations: Copello 1997 read
    in Round 6.12. It reports a second *population* (35 % low-activity,
    Po <= 0.1), not a second gate; putting it in the bell shortens sparks
    from 228 to 48-157 ms, not to 6.3. Does not rescue the use gate.
  - [ ] Whether a triad's inactivating fraction is the bilayer's (Hain
    1994/95: about half of RyRs phosphorylated; Laver & Lamb suggest
    phosphorylation or FKBP12 as the switch).

- [x] 6.12 Copello et al. 1997's low-activity channels (`ryr_mixed.
  low_activity`, `LA_READINGS`; `spark-termination --scan low-activity`).
  About 35 % of skeletal RyRs (7 of 20) have Po <= 0.1 at every Ca2+,
  activated at 70-150 uM and inhibited at 100-300 uM. Placed in the
  population bell (not in the cleft, where they could only lengthen a
  spark), they force the high-activity channels' Ca2+ gate down from
  Ki 104 to 65-106 uM and sparks from 228 ms to **48-157 ms** -- the right
  direction, an order of magnitude short of the measured 6.3 ms.
  **Round 6.11's shared-gate assumption is addressed and survives**: the
  documented spread is real and too small. Copello's heterogeneity is in
  activation and overall activity; its skeletal IC50 spread (0.16-1.1 mM)
  straddles the bell's Ki and never approaches Stern's 10 uM.
  Emergent:
  - [ ] Whether a triad's inactivating fraction is the bilayer's (carried
    over from 6.11; Hain 1994/95 phosphorylation, FKBP12).
