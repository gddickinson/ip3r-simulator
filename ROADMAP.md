# ROADMAP — what is not done

Work proceeds in rounds of one focused task each. Every round: implement,
test (`make test`, `make screenshots` if the UI changed), update the docs,
commit, push. `[ ]` planned, `[x]` done (with what it measured). The
completed Round 1 is recorded in `SESSION_LOG.md`.

**Next: Round 5, item 1** (parameter editing in the GUI), unless an emergent
Round 4 item is preferred.

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
- [ ] Where do the stride-3 local modes (#11–15 of 8TKG, κ ≤ 0.11) come
  from — which sites are weakly attached, and should the network bridge
  them (e.g. a sequence-neighbour spring) rather than just flag them?
- [ ] The gate radius along the morph is not measured: side chains ride
  their Cα rigidly, so a profile of intermediate frames would be wrong at
  exactly the gate. Needs side-chain interpolation (or a rotamer-free
  backbone-only pore measure) first.

## Round 3 — more of the publication, on the structure

Emergent (not scheduled):
- [ ] The JSD falls most at 11.5 Å, but the FEL purifying fraction falls
  after the second shell (S22 §8). Is that the instruments or the residues?
  Re-derive the §8 shell FEL table and compare residue by residue.
- [ ] "Show on structure" for the shells frames the whole tetramer, so the
  pocket is small on screen. Add a camera preset centred on one IP3 site.
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
- [ ] The grid's miss layer shows *which* genomes fail, not only how many.
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
- [ ] The Range tab stops at clade level. S23's 194 genomes (copy number,
  control verdict) could be drawn per genome inside a clicked clade, the way
  the Genomes tab draws Paper 3.
- [ ] Paper 1's family-call benchmark (24/25 recall, 31/31 specificity, the
  0.10 labelled-bait margin) and the profile calibration (11,875 agree,
  1 disagree) are not yet re-derived; `benchmark_controls/` and
  `hmm_sweep/` hold the per-record inputs.

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
  - [ ] Park/drive's stationary bell is not in the Gating panel beside DYK
    and Mak (`park_drive.bell_at` exists; at 0.2 µM IP3 it peaks at 0.32
    near 0.70 µM Ca²⁺).
  - [ ] The cluster Ca²⁺ is instantaneous mean-field. Cao integrate a
    microdomain ODE with fluo-4, which is what their inter-puff-interval and
    puff-shape results rest on. Adding it would let the app reproduce their
    IPI distribution and the amplitude-vs-N saturation.
  - [ ] Above ~0.5 µM coupling, park/drive sits at a sustained 9 % open
    rather than puffing. Is that the missing store depletion, or the model?
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
  - [ ] Salt bridges: filter D2478 is 2.5 Å from R2471 of the next subunit,
    and K2482 is 4.4 Å from D2400. Neutralise ion pairs, or use a pKa
    estimate, before trusting any charged number.
  - [ ] Why 2.4× short even uncharged? Candidates: 8TKF is not maximally
    open (a subconductance state?), the continuum fails at 3 Å, or the
    cytosolic exit is not where the profile window ends. Compare RyR1 open
    deposits, whose ~750 pS a structure-based model should also meet.
  - [ ] Ca²⁺ current under physiological ions (Vais 2010: 0.30 pA/mM
    [Ca²⁺]_ER; P_Ca:P_K = 15). The solver takes asymmetric baths already.

## Round 5 — usability

- [ ] Parameter editing in the GUI with a visible "modified" banner (checks
  already refuse to confirm against a modified registry).
- [ ] Session save/restore (structure, style, camera).
- [ ] Variants painted as spheres on all four subunits with class colours;
  ClinVar VUS stratified by conservation layer.
- [ ] AlphaFold models for the unresolved stretches, seams shown.

## Deliberately not doing

- All-atom MD in the interactive loop.
- Editing `ip3r_genes`. A discrepancy found here is reported to the user
  and recorded in `SESSION_LOG.md`; fixing it is that project's work.
