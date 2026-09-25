# ROADMAP — what is not done

Work proceeds in rounds of one focused task each. Every round: implement,
test (`make test`, `make screenshots` if the UI changed), update the docs,
commit, push. `[ ]` planned, `[x]` done (with what it measured). The
completed Round 1 is recorded in `SESSION_LOG.md`.

**Next:** Round 7.11 (science, to be chosen from the open IP3R items below; the wall charge in 3-D is the natural follow-on to 7.10). **Priorities changed 2026-09-25 (user):** the RyR work (Round 6,
`ROADMAP_RYR.md`) is parked after 6.13, and its next item (6.14, the C/V
flux ratio) waits there. Rounds 7.x alternate IP3R science with GUI
upgrades, drawing on the open IP3R items of Rounds 2-5.

## Destination

An application in which every result of the `ip3r_genes` series can be *seen*
on the receptor and *re-derived* by an independent route, beside physics
models that make the channel's function visible: IP3 and Ca²⁺ in, a gate
that opens, Ca²⁺ out.

## Round 1 — the port and the first 26 checks  [x]

- [x] Engine ported from PIEZO1 (reader, renderer, camera, GL widget,
  parameter registry). Resources imported from ip3r_genes with source hashes.
- [x] Structure measurement: C4 axis by superposition, pore profile, filter
  and gate, IP3 contacts, numbering check (stub-aware). 6DQN reproduces S0 to
  ≤ 0.03 Å.
- [x] Physics: C4-labelled ANM; DYK/Li–Rinzel gating, cell oscillations
  (window 0.36–0.63 µM), stochastic puffs (Fano 1.43 vs 0.96).
- [x] 26 findings checks across S0 and papers 1–6, each calibrated by a
  planted flip; 24 confirmed, 2 genuine discrepancies (see README).
- [x] GUI: structure, channel (with ITPR3 state comparison), modes, dynamics,
  findings, variants; scripted smoke test.

## Round 2 — gating transition on the structure  [x]

- [x] Morph between two states of the same paralog (8TKG resting → 8TKF
  activated), residue-matched (2,194 residues × 4 subunits), restrained
  (worst Cα–Cα error 0.00 Å vs 2.7 Å linear), labelled as an interpolation.
  Drawn end = deposited 8TKF as a shape to < 0.001 Å (tested).
- [x] Overlap of each ANM mode with the observed displacement. Lowest
  collective A mode of 8TKG: 0.415 (null 0.021); 20 modes 0.638 (null
  0.048); robust over strides 1–4. From the 8TKF end: lowest A mode 0.17.
  Needed a collectivity guard: at stride 4 the naive "lowest A mode" was a
  weakly attached fragment (κ 0.01).
- [x] Per-residue displacement painted on a fixed 0–25 Å scale: RIH_N 22 Å,
  MIR 19 Å, pore domain 2.5 Å, filter 0.8 Å (pore fit).

Emergent (not scheduled):
- [x] Where do the stride-3 local modes come from (`physics/network_checks.py`)?
  All five (#11–15 of 8TKG) are one flap: residue 86, the first after the
  unresolved loop 77–85, on 5 springs (median 17). The same flap takes modes
  2–5 at stride 4; strides 1–2 have none. RyR1 9R8O has its own (residue
  896, even at stride 1). Sequence-neighbour springs (1–100 × γ) change
  nothing, so the network is not bridged; the modes are named in the
  report and the Modes tab. RMSIP20 vs stride 1: 0.97 / 0.89 / 0.77 at
  strides 2 / 3 / 4. **Found on the way:** the lowest-A-mode overlap is a
  cutoff choice. At stride 1 it falls 0.48 → 0.24 over 12–21.6 Å while the
  three lowest A modes together hold at 0.66–0.68 (`--cutoff-scan`).
  RyR1's lowest A (0.12) is flat in the cutoff.
- [x] Default stride 2 (user's decision, 2026-09-24): no local modes on 8TKG,
  RMSIP20 0.97 vs stride 1, 1.2 s. Collective A together 0.669; lowest A
  0.390, best single #10 0.498. RyR1 keeps one local set (residue 1988).
- [x] Plot the A-subspace overlap as the Transition tab's headline (Round
  7.3: `TransitionOverlap.subspace`, 0.67 against 0.039).
- [x] The gate radius along the morph (`structure/morph_pore.py`): every
  heavy atom both deposits resolve interpolated (Cα offset, 70,952 atoms,
  none unmatched), axis re-found per frame. Endpoint frames measure as the
  deposits (8TKG 2.73 Å F2513, 8TKF 5.85 Å N2510; 9R8O 3.32, 9HEO 5.05);
  rigid side chains end 0.69 Å short on ITPR3, 0.19 Å long on RyR1. Gate
  opens monotonically, half-way at 0.42 (ITPR3) / 0.44 (RyR1); in ITPR3
  the constriction hands over from I2517 to N2510 at t ≈ 0.7. Lining
  side-chain chord ≤ 0.23 Å. The viewer now draws the interpolated atoms.

## Round 3 — more of the publication, on the structure

Emergent (not scheduled):
- [ ] The JSD falls most at 11.5 Å, but the FEL purifying fraction falls
  after the second shell (S22 §8). Is that the instruments or the residues?
  Re-derive the §8 shell FEL table and compare residue by residue.
- [x] "Show on structure" for the shells frames the whole tetramer, so the
  pocket is small on screen. Add a camera preset centred on one IP3 site.
  (Round 7.1: `SceneController.site_view`, Ctrl+3.)
- [ ] Paper 2's model-robustness claim ("nine claim clades held, none
  weakened" under the extra NNI round) can be re-derived the same way from
  `itpr_ml_bnni.contree` and `claim_members.tsv`, with the two trees shown
  side by side.

- [x] Paper 6 module contrast, re-derived from the deep alignments with this
  project's own column map, modules and statistics. Every published number
  reproduces: 260/246/262 tips, core − pore −0.0239/−0.0005/−0.0199,
  ITPR2's signed-rank p 0.134144 to all six figures. With the luminal loop
  counted as pore, the answer reverses in all three paralogs (+0.055/+0.044/+0.056).
  3 checks (29 total: 27 confirmed, 2 discrepancies, unchanged). Modules
  drawn as Cα traces from "Show on structure".
- [x] Ligand shells. The pocket is recomputed from the six IP3-bound
  depositions: 125 residues, 12/14/40/59 per shell, every shell identical,
  medians to 2×10⁻⁵ Å. It is carried to ITPR1/2 by an own Gotoh/BLOSUM62
  alignment, and all 250 positions agree with S22's MAFFT transfer. Every
  shell_constraint and shell_trend field reproduces. "No step at 4.5 Å"
  holds: the contact shell does not beat the second (p 0.70/0.31/0.39), and
  the largest drop between shells is at 11.5 Å in all three paralogs.
  4 checks (33 total: 31 confirmed, 2 discrepancies, unchanged). Residues are
  painted by shell; the exhibit plots JSD against distance.
- [x] Paper 2 tree viewer. `rooted.nwk` is drawn with this project's reader
  (Tree tab, and as the P2 exhibits). Paralog clades and the RyR outgroup are
  boxed, cyclostome tips marked, and nodes clearing 80/95 dotted. There is a
  one-click vertebrate zoom. Three new rederived checks, all confirmed:
  clades 19/13/19 at 100/100 with 4/2/1 unnamed tips and only the six
  cyclostomes outside; two cyclostome-only clades, each holding both
  species (4 tips first among vertebrates at 99.5/100, 2 joining
  ITPR2+ITPR3 at 83.1/77); 91 of 131 bipartitions clear both bars. The
  naive node count was 92/132 because the root's twin edge was counted
  twice. 36 checks: 34 confirmed, 2 discrepancies (unchanged).
- [x] Paper 3/4 genome × paralog grid (Genomes tab: 309 genomes × ITPR1–3
  + RyR). It has four layers (search grade, misses, S15a state, recovery
  channel) and an N50 strip on a fixed scale. It sorts by N50, class or
  name, and the contiguity bar is drawn. Three new rederived checks, all
  confirmed. Every contiguity test reproduces from the per-cell table with
  this project's statistics: the miss rate is 140/923 against 42/309
  (Fisher p 0.578), and the odds of finding the gene rise 8.10× and 19.99×
  per tenfold N50. The Mann–Whitney p is one-sided (the first checker
  doubled it). Above the registered bar there are 189 genomes and 5/563
  misses, and the split matches `spans_gene` in every cell. The recovery
  channel rebuilt from counts matches the label in 1,236/1,236 cells, and
  the reasons come to 289/186/254/15 with 179 reachable. 39 checks: 37
  confirmed, 2 discrepancies (unchanged).
Emergent (not scheduled):
- [x] The grid's miss layer shows *which* genomes fail, not only how many.
  Paper 3 says the ITPR3 lesion excess is a bird result below the bar. That
  could be drawn as another layer from `lesion_by_class.tsv` /
  `integrity_pairs.tsv`, so the reader can see it sits in the same rows.

- [x] Paper 1: presence/absence across eukaryotic clades (Range tab: one bar
  per clade, fraction of 6,928 proteomes with a call, genome-level absences
  marked). Seven rederived checks, all confirmed (`P1.absences` upgraded
  from reading the summary to rebuilding it). The numbers: 662/6,928
  proteomes and 45/135 clades. Presence from the 2,012 assignment records
  agrees with the presence table in 6,854/6,854 taxa. All 13 named phylum
  counts match (Streptophyta 0/384 vs Chlorophyta 15/48; Dikarya
  0/1,353). All 24 cells of the relaxed-sensitivity table match (PF08709
  0/26/0/16, MIR 633/4,376). S23's G3 rule, rebuilt from its sentence,
  gives the same 35 absence targets, and all 35 genome rows rebuild from
  the per-genome ledgers. Copies: 117/43 of 194 genomes, with Macrostomum
  18, Stentor 13, Dysidea 8 and Cymbomonas 3. Chase: 47 real genes, 52
  fragments. 46 checks: 44 confirmed, 2 discrepancies (unchanged).
Emergent (not scheduled):
- [x] The Range tab stops at clade level. S23's 194 genomes (copy number,
  control verdict) could be drawn per genome inside a clicked clade, the way
  the Genomes tab draws Paper 3.
- [ ] Paper 1's profile calibration (11,875 agree, 1 disagree) is not yet
  re-derived; `hmm_sweep/` holds the per-record inputs. (The family-call
  benchmark was re-derived in Round 7.7.)

## Round 4 — better physics

- [x] Mak et al. 1998 Hill-type gating model (`physics/gating_mak.py`,
  Gating panel model choice, `gating --model mak`). All seven constants
  were read from the paper's text before registration. Both models are
  measured with one ruler (`physics/bell.py`). From 33 nM to 10 µM IP3,
  Mak moves half-inhibition 6.22× and half-activation 1.016×; DYK moves them
  2.76× and 2.01×. A planted K_act dependence is caught. The plateau is
  0.77–0.81 over 1–20 µM Ca²⁺ (paper: "≈0.8").
  Emergent: the model has no kinetics. A kinetic scheme fitted to the same
  data (e.g. Siekmann's park/drive, the next item) would let the Mak
  receptor drive the cell and puff models.
- [x] Park/drive puff model (Siekmann 2012 + Cao 2013 gating variables;
  41 constants read from the authors' code, Cao 2014 Text S1), in the same
  mean-field cluster as DYK and read with one ruler (`puff_compare`). Over
  30 s at 0.2 µM IP3, the DYK cluster reaches half its channels (10/20) in
  ≤ 1 event at every coupling from 0 to 2 µM. Park/drive does so 10–17 times
  at 0.09–0.32 µM (Fano 2.8 vs DYK ≤ 1.32), with a valley in the event-size
  distribution. Seeds 1–3 at 0.1 µM: 12–14 against 0–1. Resting open
  fraction is 1.7 % against 5.2 %. The split step reproduces the clamped
  stationary P_open to 0.2 %, and a 5 ms step is caught (10.5 %).
  Emergent:
  - [x] Park/drive's stationary bell is in the Gating panel beside DYK
    and Mak (Round 7.3; it passes Mak's flank test, 1.09× / 42.8×).
  - [x] The cluster Ca²⁺ is instantaneous mean-field. Round 7.2 adds Cao's
    microdomain with fluo-4 (`physics/microdomain.py`): the IPI trend with
    a_h42 and the dF/F0 bend at N ≈ 12 are reproduced.
  - [x] Above ~0.5 µM coupling, park/drive sits at a sustained 9 % open.
    Round 7.2: the model, not the store. The state survives the
    microdomain's kinetics in a fixed bath, and a free store raises
    activity.
- [x] Unitary current from the pore profile (1-D drift–diffusion ported from
  PIEZO1, and wall charge from each deposit's own side-chain atoms), across
  the ITPR3 state panel. Six of seven states are sterically shut (r_free
  ≤ 1.03 Å). Activated 8TKF gives 65 pS neutral (series check 64 pS;
  25–150 pS over diffusivity 0.25–1× and ion radius 1–2 Å) and 33 pS
  charged. The measured values are 358 ± 8 pS (Mak 2000; Vais 2010 cites it
  as 370) and 545 ± 7 pS (Vais 2010), so the model is 2.4× short at best.
  The lining rings of alternating sign act as junctions in series and
  *lower* conductance. Acidic rings alone give 174 pS. The charged number
  is fragile (7–151 pS over margin/smoothing; peak partition 18.6 M, above
  the packing ceiling).
  Emergent:
  - [x] Salt bridges cancelled before the wall charge counts
    (`physics/salt_bridges.py`; Barlow & Thornton's ≤ 4 Å N–O, one-to-one,
    closest first, partners searched in the whole deposit). A fourth
    reading, *paired*, appears in the CLI and the Channel panel. On 8TKF the
    four D2478–R2471′ bridges (2.45–2.58 Å) drop out: −8 → −4 e, and
    33 → 23 pS (8–54 pS swept). Over cutoffs of 3–6 Å the reading is
    23–38 pS, because D2518–R2524′ (4.2–4.3 Å) and K2482–D2400 (4.3–4.4 Å)
    sit just past the line. The 18.6 M peak is on the D2518/D2522 rings and
    does not move. Pairing does not rescue the charged model.
    **Later the same day (Round 6.1): RyR1 refutes the paired reading.**
    D4899 is D2478's homologue, bridged the same way, and D4899Q cuts
    RyR1's conductance to 0.20×, where pairing predicts 1.00×. *Charged*
    is the better reading; *paired* stays only as a bound.
    Emergent:
    - [x] A pKa estimate for the eight D2518/D2522 carboxylates in a 4.4 Å
      lumen (Round 7.4): charged at pH 7.3 by both routes (network pKa
      4.1, PROPKA 5.3). The ring is 8.8 Å across and R2524′ sits beside it.
  - [x] Why 2.4× short even uncharged? Candidates: 8TKF is not maximally
    open (a subconductance state?), the continuum fails at 3 Å, or the
    cytosolic exit is not where the profile window ends. Compare RyR1 open
    deposits, whose ~750 pS a structure-based model should also meet.
    Round 7.6: mostly the inscribed circle. The real lumen in 3-D conducts
    1.3–2.0× more; not a substate (7T3T agrees); not the window.
  - [x] Selectivity and the Ca²⁺ current (`physics/selectivity.py`,
    `python -m ip3r selectivity`). Vais 2010's three lum-out protocols run
    through the solver in their own solutions and are read with their GHK
    Eq. 1. A ratio does not depend on the unmeasured diffusivity, so it
    tests the wall. Needed: each mouth in Donnan equilibrium with its own
    bath (TMS; symmetric results unchanged). Calibrated against Planck, TMS
    at four charges, the Nernst limit, an excluded NMDG⁺, and Vais's own
    arithmetic. On 8TKF: P_Ca:P_K 0.17 / 0.00 / −0.07 / 0.69
    (neutral / charged / paired / acidic only) against 15.2, and i_Ca
    ≤ 0.046 pA/mM against 0.30. The lining bases (K2529, K2482) are Ca²⁺
    barriers under local Donnan. Even acidic-only, the charge is in rings,
    and the calibration shows a ring cannot give Ca²⁺ selectivity (< 0.6 at
    any charge) where a charged tract can (60 at −30 M). Every charged
    reading is anion-tight (P_Cl:P_K ≤ 0.05 against 0.27). The model obeys
    GHK to 4 %, so it cannot show Vais's 8–10× i_Ca shortfall below GHK.
    Emergent:
    - [x] The same question on RyR1's open deposit, the control (Round
      7.4, Xu 2006's own protocol): 0.46 against 7.0, and the mutants'
      order is missed.
    - [x] Screening in the wide vestibule (`physics/radial_pb.py`,
      `selectivity --closure radial`). Cylindrical Poisson–Boltzmann across
      every slice, with Gauss's law at the wall, replaces local Donnan. Each
      species sees its own cross-section partition. Calibrated against the
      Donnan limit, the Debye–Hückel Bessel closed form (0.1 %), exact
      discrete Gauss and grid convergence. At K2529 (R 9.9 Å, R/λ_D 1.7)
      the Ca²⁺ barrier falls from +155 to +111 mV, a factor of ~6 in
      partition. Charged P_Ca:P_K goes from 0.00 to 0.04 (0.01–0.16 over
      ε 80–10), and g from 33 to 46 pS. Acidic-only stays at 0.69, since
      those rings sit in slices narrower than λ_D. **The closure is not the
      Ca²⁺ barrier.**
      Emergent:
      - [x] The lysines' own protonation (Round 7.4): charged at any
        permittivity (network pKa 12.9 and 10.9).

## Round 5 — usability

- [x] Parameter editing in the GUI (Help → Parameters, `Ctrl+Shift+P`) with
  an amber full-width "modified" banner, driven by registry change
  listeners. Import/export use the `IP3R_PARAMETERS` format; nothing is
  persisted between sessions. Listeners also clear the two
  parameter-dependent memo caches the checks use (`_summary`, the shell
  pocket). Without that, a measurement made under an edit was served to a
  check after reset (the test fails with the subscription removed). The
  smoke test edits, sees the banner, sees a check refuse, and resets.
  Emergent:
  - [x] Panels that read a parameter when they are built (spin-box
    defaults) do not follow an edit. Audit them, and either re-read on
    change or document it. Done in Round 7.5.
- [x] Session save/restore (File → Save/Open session, `Ctrl+Shift+S`/`Ctrl+O`,
  `--session`). A session holds the view: deposit, style, colour, layer,
  subunits, sites, pore, camera, tab, and the transition spec (end, fit,
  method, frame). It holds no coordinates or results. It does record the
  parameter overrides it was saved under, and restoring asks before applying
  a different set (applying re-measures). The smoke test saves an 8TKG→8TKF
  view at frame 5 under an edit, moves to 6DQN with defaults, and restores.
  All 17 view fields come back, the camera to 1e-9. It fails when the camera
  restore is removed (checked). It caught one bug on the first look: the
  deposition list stayed on the previous deposit.
  Emergent:
  - [x] Mode animation and Dynamics-panel settings (IP3, coupling, model) are
    not in a session. Add them if a saved view turns out to need them.
    Done in Round 7.5.
- [x] Variants as spheres on every visible subunit, in class colours
  (Variants tab, "Draw on structure"; refused on a deposit in another
  numbering). The spheres follow a morph frame. With a layer chosen, each VUS
  is placed against its gene's labelled medians (Paper 5 §8) and takes its
  stratum's colour. New rederived check `P5.vus_stratification`: all 12
  gene × layer rows of `vus_stratification.tsv`, every count, median and
  fraction, rebuilt from `variants.tsv` and the per-residue tables. On the
  family layer, 10/6/11 % of VUS reach the P/LP median and 16/36/52 % sit at
  or below the B/LB median. The table includes the curated UniProt P/LP
  records, which the AUC test excludes; run ClinVar-only it is identical. The
  viewer's resource route reproduces the table too (tested). 46 checks: 44
  confirmed, 2 discrepancies. The earlier "46 checks" was a miscount of 45.
  Emergent:
  - [x] The variants view (class, layer, drawn) is not in a session
    (Round 7.5).
  - [x] Per-paralog: ITPR2's P/LP median is one position's score. Draw the
    thresholds' uncertainty (bootstrap the medians) so a stratum near a
    median is shown as such.
- [x] AlphaFold models for the unresolved stretches, seams shown
  (Representation → Completeness; `python -m ip3r graft`). Each stretch is
  placed by a local Cα fit on sequence-matched, unstubbed anchors either side
  (a terminus on one side, only when asked). It is drawn as its own
  pLDDT-coloured layer, and each seam is a bond, red when broken. The fill
  follows morph frames. No measurement runs on it. The model is chosen by
  identity by number: the 7 ITPR3 deposits match `AF-Q14573-F1` at
  99.1–100 %. AlphaFold DB has only isoforms for ITPR1/2, so 7LHF (22 %) and
  9YKK (18 %) are refused. On 8TKG: 56 stretches, 1,592 residues, mean pLDDT
  38; 12 of 224 seams broken; 108 residues clash. Calibration: 31 stretches
  that another deposit misses but 8TKG/8TKF resolve were hidden and filled.
  Median 1.46 Å, against 5.59 Å for a straight line and 10.5 Å for a global
  fit. The fill beats the line 28/31, and pLDDT predicts the error (ρ −0.79).
  The join tolerance was moved from 4.5 to 5.5 Å because true seams reached
  5.42 Å. The smoke test caught a fill drawn at the deposit instead of the
  morph frame on its first run.
  Emergent:
  - [x] Long stretches (Round 7.9): right at any length where AlphaFold is
    confident; the one very-low stretch with a truth is 58 Å off.
  - [x] 7LHF filled from rat isoform 8 through an alignment (Round 7.9),
    SI and SII left unfilled and named.
  - [ ] Fills on neighbouring subunits are not checked against each other,
    only against the deposit.
  - [ ] The Completeness choice is in a session; the transition spec does not
    say whether the fill was drawn while the frame was saved (it is simply
    rebuilt).

## Round 7: IP3R first, and the GUI (planned 2026-09-25)

The user reprioritised after Round 6.13: IP3R science and GUI upgrades
come before any more RyR work. GUI rounds alternate with science rounds,
so each new result lands somewhere it can be seen. (Rn) marks an open
item carried over from round n; it stays listed there too.

- [x] **7.1 GUI: viewport and navigation.** The fit itself was right. The
  viewport was the problem: the Analysis dock's minimum was 851 px (one
  row of Genomes controls), so the view got 209-319 px of a 1512 px
  window. Wide rows now split, wide combos may shrink, and the docks
  start at 380/580 px. The viewport gets 36 % of the window and the
  molecule reaches 0.94 of it. The fit is redone on resize and on a
  subunit change until the user moves the camera. Ctrl+3 and "Show" for
  the contact and shell checks centre on one IP3 site: side-on from
  outside, the 15 Å pocket filling the view, the front clipped at the
  pocket, and the clip saved in sessions. Deposits without IP3 fall back
  to the whole structure and say so. The deposition list is a tree by
  family with RyR1 collapsed. The Puffs rows for Mg²⁺, K_Mg,A and
  triggered sparks are hidden for IP3R, and IP3 is hidden for RyR1. All
  of this is held by the smoke test (`scripts/screenshot_view.py`).
- [x] **7.2 Science: the IP3R puff in a microdomain.** Cao 2014's pools
  (every constant from their code) with fluo-4 in the microdomain (Cao
  2013 Eq. 12, Shuai 2006's constants), driving the park/drive receptors
  (now `puffs_pd.ReceptorCluster`, shared and bit-identical). Puffs are
  read from F/F0 and IPIs fitted by Thurley's Eq. 14 by maximum likelihood.
  *IPIs* (N 10, 0.1 µM, 900 s): as a_h42 goes 0.1 → 5 /s the rate rises
  5.5× and the CV goes 0.79 → 0.93 (Cao: 0.65-0.95). At the slow end λ 0.16
  and ξ 0.61 land on Cao's values, but the refractory fit is not
  significant on 114 intervals (LR 2.2). *Amplitude vs N*: at Cao 2014's
  release, puffs stay under K_d and nothing bends. At 2.5× (mean blip
  dF/F0 1.33-1.66, Cao's 1.6) dF/F0 bends at N ≈ 12 (0.5 → 0.17 per
  receptor) while Ca2+ bends less. *Sustained 9 %*: five conditions from
  the mean field to a free store. The state survives kinetics and dye in a
  fixed bath (6.5-8 %). A free store *raises* activity against a held one,
  and only the whole cytosol filling (27-60 µM) lowers it. So: the model,
  not the store. Incidental: `pd.ca_mouth` = 120 µM is the code's rule at a
  store of 100 µM, while the code's store is 449 µM (539 µM mouth). Deepening
  the mouth moves the open fraction by < 1 point.
  Emergent:
  - [ ] Cao 2013's Table S2 (its release and J_decrease) is still unread.
    It would replace `domain.blip_scale`, which is matched, not read.
  - [ ] Cao 2013's Eq. 10 λ_h42 (a Ca2+ switch at 20 µM, V = 100 /s) against
    the 2014 code's open/closed switch used here.
  - [ ] Seeds: the a_h42 scan runs one seed per point (ξ at 0.27 /s is out
    of line with its neighbours).
- [x] **7.3 GUI: the gating models side by side.** The Gating panel
  (`ui/gating_panel.py`, split out of `dynamics_panel`) adds park/drive's
  stationary bell (with the parked fraction) and a "three models side by
  side" view at Mak's two IP3 levels. `compare_flanks` now measures all
  three on one ruler, 33 nM → 10 µM: park/drive moves half-activation
  1.09× and half-inhibition 42.8× (Mak 1.016× / 6.2×, DYK 2.0× / 2.8×), so
  it passes Mak's test although its rates were not fitted to it. The
  Transition tab's headline is now the A subspace
  (`TransitionOverlap.subspace`: the collective A modes added lowest first,
  beside the same-symmetry null and the √(A share) ceiling; 8TKG → 8TKF
  0.67, random 0.039). A selector below it picks the gate path, every
  mode, or the element means: four plots at half the dock's width did not
  fit. For park/drive, the Puffs panel gets a microdomain box
  (`ui/puffs_domain_view.py`, `domain.gui_duration`): F/F0 with the puffs
  marked, the number open, and IPIs against Thurley's Eq. 14. At the
  panel's defaults (IP3 0.2 µM, N 20, 30 s): 44 puffs, CV 0.68, and a
  significant refractory period (LR 8.4).
- [x] **7.4 Science: protonation in the IP3R pore.** Two independent
  routes: a Tanford–Kirkwood network (`physics/pka.py`: Thurlkill/Fitch
  model pKas, Mehler–Solmajer ε(r), Monte Carlo with a ring heat-bath move,
  held to exact enumeration) and PROPKA 3 (`physics/pka_propka.py`). Both
  keep every lining group of 8TKF charged at pH 7.3. K2482/K2529 stay
  charged even at a uniform ε of 4, where the acid rings fall to −0.5 to
  −0.8. P_Ca:P_K is ≤ 0.00 under every reading. A bound that needs no pKa
  (`protonation --corners`: 64 ring on/off states, plus interior samples)
  tops out at 0.69 with both lysines neutral, and at 0.05 with them
  charged, against 15.2. **The control fails the same way**: RyR1 (9HEO)
  under Xu 2006's protocol gives 0.46 (0.07–0.46 over readings) against
  7.0, and the model gets the mutants' order wrong (D4899Q ×0.74 against
  ×0.14 measured; E4900N ×0.24 against ×0.64). So the continuum model is
  missing the selectivity physics, and protonation is not the IP3R
  wall's problem.
  Emergent:
  - [ ] Charge–space competition: finite ion size (a local excess chemical
    potential, e.g. Bikerman or a mean-spherical-approximation term) in
    the drift-diffusion pore. Test it first on RyR1, where Xu's six
    P_Ca:P_K values and conductances are the calibration (Gillespie's
    PNP-DFT reproduces them), then on 8TKF against 15.2.
  - [ ] PROPKA buries RyR1's E4900 (pKa 8.0), which Xu's E4900N shows is
    charged. Is it the deposit's rotamer or PROPKA's desolvation?
- [x] **7.5 GUI: sessions and live parameters.** Sessions now hold the
  Dynamics controls (gating model, oscillation, puffs, microdomain, the
  sub-tab), the animating mode and its amplitude, and the Variants view
  (`Session.dynamics/modes/variants`, flat name → scalar, no format bump:
  older files open). Controls only: no simulation is re-run. A mode is
  recomputed and restarted after any transition frame, because a frame
  stops it. A value that cannot be set (unknown key, out of range) is
  noted in the status line, not forced. The audit found one real fault:
  the Puffs cluster size and coupling were copied into spin boxes when the
  receptor was chosen, so editing `puff.n_channels` changed nothing the GUI
  simulated. Now such controls (`ui/view_state.Seeded`: those two, the
  microdomain run length) follow an edit unless the user typed their own.
  The Gating plot and the displacement/shell colours with their legend are
  redrawn. Results the user ran keep their values, as the banner says.
  The Range panel's minimum is a display filter, not
  `range.absence_min_proteomes`. The smoke test sets and then disturbs
  every new field, restores, and compares the whole session. It animates a
  mode, restores it, edits the seeding parameters (a restored 90 s
  duration must not follow; the coupling must), and resets.
  Found on the way: the Transition tab's Stop, pressed during a mode
  animation, froze the atoms mid-swing. It now puts the deposit back.
- [x] **7.6 Science: the conductance shortfall.** Each candidate was tested
  by one measurement (`python -m ip3r shortfall [--scan]`), with RyR1's 9HEO
  as the control. The new instrument is a 3-D ohmic solve of the voxelised
  ion-accessible lumen (`structure/pore_volume.py`, `physics/ohmic3d.py`).
  It uses the same electrolyte as the 1-D model and is calibrated on
  Hall's cylinder (0.97 at 0.5 Å), two pores, a blind hole, a sideways
  exit and the membrane seal. The results:
  - **Substate: no.** Schmitz 2022's independent active 7T3T (now a
    registered `open_control`, kept out of S11's state panel) reads 85 pS
    against 8TKF's 106. Both are 3.4–4.2× short at the registered
    diffusivity, and so is RyR1.
  - **Exit window: no.** 6–24 Å moves the 1-D reading by 2 %, access is
    2 % of R, and the 3-D solve (no window, bath on the box sides) moves
    < 0.5 % with the box.
  - **Continuum geometry: yes.** The inscribed circle leaves out the
    lumen's corners. The real shape conducts 1.3–2.0× more (8TKF 65 →
    106 pS, 9HEO 136 → 278).
  - At bulk diffusivity and a 1 Å K+ exclusion, RyR1 gives 787 pS against
    801 measured, and ITPR3 gives 261–278 pS against 358–545 (1.3–2.1×).
  Emergent:
  - [ ] The wall charge in 3-D. In 1-D, rings of alternating sign act as
    junctions in series and lower g. Do they still, with the charges in
    the lumen's corners? A Donnan-partitioned conductivity per voxel
    first, then (if it matters) 3-D PNP. RyR1's mutants are the
    calibration.
  - [x] Show the lumen and its potential in the viewer (Round 7.10).
- [x] **7.7 Publication: re-derive what is still read.** Five new checks,
  all confirmed (49: 47 confirmed, 2 discrepancies, unchanged).
  - **S22 §8** (`P6.shell_rates`, `P6.module_rates`): S17's FEL sites
    joined to our own pocket and modules. The literature core (Bosanac
    2002, ITPR1 224–604) is now a registered module, carried by our
    alignment onto S17's transfer exactly, so `P6.module_map` compares all
    twelve spans. All 21 rows reproduce, including the BH q-values, which
    need every row of their family.
  - **`--bnni`** (`P2.bnni_robustness`): the 13 claim sets are rebuilt from
    census prefixes. A core is the largest *pure* clade: "every tip naming
    it" gave 15/11/18, not 13/10/16. Both trees are rooted with a new
    `newick.reroot`. 9 of 10 clades held; the ITPR1 core, never strong,
    went 47.8/95 → 47.5/73. The Tree tab's "Beside --bnni" shows the two
    trees.
  - **The benchmark** (`P1.bait_margin`, recomputed; `P1.benchmark_counts`):
    every control aligned pairwise to the six human baits from the
    committed sequences. Full identity: 31/31 calls, gap 0.582 (0.601),
    iplA +0.048 (+0.065) inside ±0.10. Covered identity is not
    re-measurable pairwise: unrelated long sequences floor at 0.25–0.28,
    and there iplA reads −0.002. The fly Itpr miss (0.342 vs the 0.35 bar)
    reads 0.361 pairwise, so the 0.008 is inside aligner noise. Scores
    rebuilt with registered points and our own caps: 24/25, 31/31, 6/6.
- [x] **7.8 GUI: publication views.** 52 checks: 50 confirmed, 2
  discrepancies (unchanged). Round 7.7's "49" was a miscount of 51.
  - **Lesions** (`P3.lesion_strata`, new, rederived). S15b §8's class
    strata are rebuilt from `integrity_loci.tsv` alone: all 744 matched
    pairs, all 38 strata with p and q, and all five bar splits reproduce.
    ITPR3 Aves is 25:2 (q 4.52e-05) against Actinopteri 7:6; below the bar
    19:2, above 6:0. The Genomes tab's new lesion layer draws those pairs,
    and "Show" opens it on the birds in N50 order.
  - **Range per genome.** View → "Genomes (S23)", or double-click a clade.
    Each of the 194 genomes is filed under its S20 clade by taxid (109),
    else phylum or class; 13 sit in phyla S20 never swept and are shown as
    such. Each is drawn with its control verdict, copy-ledger status, its
    contiguity against its own bar, and its gene models on a fixed 0–20
    scale.
  - **VUS thresholds.** An exact order-statistic interval, not a
    bootstrap: a percentile bootstrap of ITPR2's single P/LP position has
    zero width. At 95 %, the P/LP medians of ITPR2 (1 position) and ITPR3
    (5) are unbounded. No ITPR2 or ITPR3 VUS is firmly pathogenic-like on
    any layer, while ITPR1 keeps 136 of its 185 on deep. The B/LB medians
    are bounded for all three genes.
- [x] **7.9 Fills.** Long stretches calibrated; 7LHF filled through an
  alignment.
  - **7LHF by alignment.** The deposit's whole construct
    (`_pdbx_poly_seq_scheme`, unbuilt residues included) is aligned to each
    model. Rat isoform 8 gives 1.000 identity over 2,681 pairs. Only SI
    (318–332) and SII (1693–1732) are unpaired, and stretches touching them
    are skipped by name. The result is 40 stretches and 1,264 residues, each
    the construct's own residue. The bar `graft.align_min_identity` (0.99)
    admits only the same protein: human ITPR1 isoform 4 gives 0.988, and
    9YKK's best is 0.713, so it is still refused.
  - **Length does not break a fill.** Hidden resolved windows of 10–60
    residues fill to 0.5–1.5 Å (line 4–16 Å) in 8TKG, 8TKH and 7LHF. The
    seams of 802 right fills stay ≤ 5.22 Å, so the 5.5 Å tolerance holds.
    But those residues are ordered (pLDDT ~82).
  - **Low confidence does.** No deposit resolves a pLDDT < 50 run with
    anchors of its own. Islands do the job instead: a resolved run between
    two gaps is hidden and the whole 27–137-residue span filled. Over 28
    islands, the 26 at pLDDT ≥ 70 land at a median 1.26 Å (line 17.8), all
    better than the line. The only island below 50, 8TKH 926–943 (pLDDT
    31), is 57.9 Å off, against 24.4 Å for the line. The fill summary now
    counts the residues below 50 (1,308 of 8TKG's 1,592) and warns that
    they are not positions.
  Emergent:
  - [ ] A Completeness choice that leaves out the very-low residues, or
    draws them as a tether rather than a chain. The one test says they are
    not positions, but it is one stretch.
  - [ ] Two fills of the same span in different deposits (e.g. 8TKG vs
    8TKF) as a consistency measure where no truth exists.

- [x] **7.10 GUI: the lumen and where the voltage falls.** Channel panel
  → "Draw the lumen": Round 7.6's voxelised lumen, solved on a worker, is
  drawn as a surface (`render/lumen_mesh.py`) coloured by the 3-D Laplace
  potential on a fixed 0–1 ramp. It is hidden on any morph or mode frame,
  because it was solved on the deposit. The panel plots the lumen's area and
  the normalised φ along S0's window, 3-D against the 1-D model's inscribed
  circle (`physics/lumen_field.py`, `python -m ip3r lumen`). Calibrated on
  a cylinder (a linear drop, Hall's access share outside), a neck (the 1-D
  share by hand on both routes) and a blocked pore. **Found:** the corners
  that give the 3-D shape its 1.3–2.0× conductance do not move the field.
  In 8TKF, 7T3T and 9HEO, 99 % of the voltage falls in the window, and the
  half-drop point is within 1.3 Å of the 1-D model's. The filter's ±3 Å
  holds 32/32/24 % in 3-D against 35/35/25 % in 1-D. So the 1-D field
  profile is sound where its magnitude is not.
  Emergent:
  - [ ] Colour the lumen by the charged reading's potential once the wall
    charge is solved in 3-D (Round 7.6's open item); the neutral φ is only
    geometry.

## Round 6 — ryanodine receptors

In `ROADMAP_RYR.md`: RyR1 structures, gating, sparks in the cleft, Mg²⁺,
the couplon under voltage clamp, SR depletion, two-site inactivation.

## Deliberately not doing

- All-atom MD in the interactive loop.
- Editing `ip3r_genes`. A discrepancy found here is reported to the user
  and recorded in `SESSION_LOG.md`; fixing it is that project's work.
