# SESSION_LOG

Newest last. What was done, and why.

## 2026-09-23 — Session 1: the port, the physics, and 26 checks

**Asked:** an IP3R visualisation/simulation app like the PIEZO1 simulator,
able to illustrate, demonstrate and confirm the `ip3r_genes` results, in its
own repository, easy to continue in later sessions.

**Built.**
- Project skeleton, `ip3r_sim` conda env (cloned from `piezo1`), parameter
  registry with provenance gate (33 parameters; 14 references).
- Ported unchanged from PIEZO1: mmCIF reader, `Structure`, camera, impostor
  primitives, scene, splines, mesh builders, shaders, GL widget (HUD removed).
- `scripts/sync_genes.py` imports sequences, domain map, functional sites,
  per-residue conservation, variants and the structure registry from
  ip3r_genes, recording each source's SHA-256; `--check` reports drift.
- Structure measurement, ANM with C4 irreps, DYK/Li–Rinzel gating, closed-cell
  oscillations, stochastic puffs, state-panel comparison.
- 26 checks, GUI with six panels, CLI, 79 tests.

**Measured.**
- 6DQN: our superposition axis agrees with S0's centroid axis to 0.02°; C4
  residual 0.051 Å (S0 0.058); filter 5.08 Å (5.06), gate 2.55 Å (2.52), the
  same lining residues and the same ten IP3 contacts at all four sites.
- ITPR3 state panel: gate 1.95–2.73 Å in every state except activated 8TKF,
  5.85 Å — the gating transition visible at the pore.
- Li–Rinzel oscillation window 0.36–0.63 µM IP3 (sustained only).
- ANM of 6DQN: 2,916 sites, solved in 0.5 s, every mode a clean A/B/E.

**Checker bugs found and fixed (the PIEZO1 lesson, again).**
- Structure validity check gave up after 20,000 lines and rejected 7LHF and
  8TKG, whose headers are longer.
- Numbering check failed 6DQN (94.8 %) because of a 113-residue stubbed
  segment (1434–1546, backbone + CB, unassigned sequence); stubs are now
  excluded and mismatch segments reported.
- Variant AUCs: counted variants, not positions, and included curated UniProt
  rows; S17 scores ClinVar *positions* with an occupancy floor.
- Ligand shells: excluded hydrogens; S22 counts them. Reproduced exactly
  once matched — and the difference became its own check (below).
- P4: counted the RyR sister cell as a paralog.
- Oscillation detector counted damped spirals past the Hopf point (window
  first read 0.35–0.73).
- GL viewport cropped under `QWidget.grab()` at a different pixel ratio.
- My own prose overclaimed that IP3 tunes only the inhibitory flank in the
  DYK model; the test measured 2.4× vs 1.8× and the prose was corrected.

**Genuine discrepancies in ip3r_genes (reported, not edited).**
1. `P6.contacts_heavy_atom`: under S0's heavy-atom definition R503 is 4.78 Å
   (8TKG) and 4.83 Å (8TKH) from IP3; S22's "ten contacts in all six" holds
   only because S22 counts hydrogens.
2. `P5.report_both_metrics`: the S17 report's "on both metrics and in all
   three paralogs" is contradicted by its own table for ITPR1 and ITPR2 on
   the JSD; the paper's Results state it correctly.

**Next:** ROADMAP Round 2 (the resting → activated transition on the
structure, and whether the lowest A mode points towards it).

## 2026-09-23 — Session 2: Round 2, the resting → activated transition

**Sync.** Both repos were up to date and `make sync-check` was clean, so no
verdicts moved (24 confirmed, 2 discrepancies, unchanged).

**Built.**
- `structure/transition.py`: a residue-matched basis for two deposits of one
  paralog. It requires both deposits in human numbering (7LHF is refused), and
  keeps only unstubbed, sequence-matching residues resolved on all eight
  chains. The end is superposed onto the start *as deposited*, so the path
  can be drawn over what is on screen (the PIEZO1 36 Å lesson). Off-basis
  atoms (IP3, lipids) ride the nearest site in space.
- `structure/morph.py`: ported from PIEZO1 for C4, with restrained and
  linear methods. The modal method was not ported; here the ANM *scores*
  the move instead of driving it.
- `physics/transition_modes.py`: overlap, cumulative overlap, exact
  rigid-body removal (makes the result fit-independent; tested), C4 isotypic
  split, and a null for a random direction of the same irrep make-up.
- `ModeSet.collectivity` (Brüschweiler κ). `first(irrep)` now skips local
  artefacts. The Modes tab shows κ.
- Transition tab (preset 8TKG → 8TKF, slider/play, fixed-scale displacement
  colouring, element and overlap plots), plus a `transition` CLI command.
- 6 new parameters and 2 new references (DOIs checked against Crossref).
  Tests: 79 → 95.

**Measured (8TKG → 8TKF, pore fit, stride 3).**
- RMSD 3.14 Å on the pore domain, 16.11 Å overall. The cytosolic RIH_N and
  MIR domains move 22 and 19 Å on average; the channel moves 2.5 Å and the
  filter 0.8 Å.
- The lowest collective A mode (#5) has overlap 0.415, against a null of
  0.021; 20 modes reach 0.638 (null 0.048). The resting network points
  towards activation. From the 8TKF end the lowest A mode reaches only 0.17.
- The displacement is 100.0 % A-symmetric. That is inherited from
  C4-imposed maps, and the docs say so.

**Where I was wrong first.**
- The null was a random direction in all 3N−6 dimensions. Because the move
  is purely A, the right null is a random A direction; it is now
  symmetry-matched.
- **The lowest-A-mode overlap was not robust**: 0.42 at stride 3 but 0.008
  at stride 4. The cause was a weakly attached network fragment (κ 0.01,
  near-zero eigenvalue) taking the "lowest A" slot. κ separates artefacts
  (≤ 0.11) from collective modes (≥ 0.27), so the threshold was set at 0.2.
  With that guard the result holds over strides 1–4 (0.39–0.49). At stride 5
  the network fails outright: it reaches only 0.03, at the null, and the
  report shows that.
- Registered the colour-scale top at 20 Å with a note claiming it was above
  the 95th percentile, but the measured 95th percentile is 24.1 Å. The top
  is now 25 Å and the note carries the measured values.
- Two of my own synthetic tests were wrong (an arbitrary chord threshold;
  expecting a random A field to have no rigid part, when z-translation and
  z-rotation are themselves A). The code was right in both cases.

**Next:** Round 3. Two emergent items were added to ROADMAP: the source of the
stride-3 local modes, and why the gate radius along the morph is not yet
measurable.

## 2026-09-23 — Session 3: Round 3, item 1 — Paper 6's module contrast

**Sync.** Both repos were up to date and `make sync-check` was clean, so no
verdicts moved.

**Built.**
- `core/modules.py` rebuilds the ligand core (the span of the ten contacts)
  and the pore module (PF00520 less the luminal loop, and the whole of
  PF00520 as the sensitivity case) from the imported sites and domain map.
  Each module is validated against what it must and must not contain, as S22
  does, so a wrong span raises instead of being drawn.
- `analysis/module_contrast.py` computes per-tip identities in each module
  from `aln_ITPR*.fasta`. The residue → column map is built by walking the
  reference row, and is refused unless its ungapped sequence is our UniProt
  sequence. S17's `deep_col` is not used, so the column mapping is
  independent too.
- `stats.sign_test` (exact, log-space) and `stats.signed_rank_test` (normal
  approximation, tie-corrected, no continuity correction). Both are
  calibrated: the variance was checked against a permutation null.
- Three checks: `P6.module_map`, `P6.module_contrast` and `P6.loop_reverses`.
  Each tests the *pattern* the prose states (which paralogs are
  significant, and in which direction) as well as the table's numbers.
  `loop_reverses` is calibrated with an *input* plant: every ITPR3 tip is
  given the reference's channel domain.
- Two new parameters: `ligand.module_min_coverage` (S22's 0.5) and
  `check.log_p_tol` (0.02 decades).
- GUI: "Show on structure" now works for residue-keyed checks with no
  structure of their own (drawn on the displayed structure, in human
  numbering only). Modules are drawn as Cα traces with per-atom highlight
  colours. The smoke test asserts that two colours are drawn.
  Tests: 95 → 112.

**Measured.** Every field of the 6 primary and loop-included rows
reproduces: tip counts, drop counts, means, sign counts and ties exactly,
and p to six figures. The published result stands. Paired per orthologue,
the pore leads in ITPR1 (223 vs 32) and ITPR3 (218 vs 42); ITPR2 is level
(p 0.134); and all three reverse when the luminal loop counts as pore. The
answer depends on one boundary, as the paper says.

**Where I was wrong first.**
- My first module colours (blue core, orange pore) disappeared into the
  blue → red conservation ramp in the screenshot. They are now green and
  magenta, outside that ramp.
- The alignments are memoised, so without an explicit cache clear the
  calibration test's "declared sources are complete" proof would have passed
  on the real alignments. The test now clears the cache (`clear_caches`).

**Not done.** The `ibc_literature` sensitivity definition (ITPR1 224–604,
transferred through S17's pairwise alignment) is a literature boundary, not
a rule, and is not rebuilt here. `P6.module_map` skips it, saying so.

**Next:** Round 3, item 2 (ligand shells: conservation against all-atom
distance to IP3).

## 2026-09-23 — Session 4: Round 3, item 2 — ligand shells

**Sync.** Both repos were up to date and `make sync-check` was clean, so no
verdicts moved.

**Built.**
- `structure/shells.py` implements S22's rule on this project's reader: the
  all-atom distance to the IP3 on the residue's own subunit, the best over
  the subunits, and the median over the six depositions. Shell edges are
  registered parameters, and the first edge is the existing contact cutoff.
- `core/pairwise.py` is a Gotoh affine-gap aligner (BLOSUM62, end gaps
  free, vectorised by row, 0.12 s for ITPR3 × ITPR1). It is needed because
  S22 carried the pocket to ITPR1/2 through S17's MAFFT transfer, and
  reusing that transfer would not have been an independent route. It is
  calibrated against a cell-by-cell three-state DP on 24 random cases, and
  the traceback's alignment re-scores to the reported score. A second test:
  ITPR3's ten contacts land on ITPR1's ten.
- `analysis/shell_constraint.py` and `checks_shells.py` add four checks.
  `stats` gains `mann_whitney_greater` and `spearman`, both calibrated
  against scipy in the tests.
- New parameters: `ligand.shell_second_edge`/`shell_third_edge`/
  `shell_radius`, `align.gap_open`/`gap_extend` (EMBOSS needle's; three new
  references), and `check.alpha`. `checks_modules` also used a literal 0.05,
  and now reads the parameter.
- GUI: a "Distance to IP3 (S22 shells)" colouring, with four discrete
  colours sampled from the fixed ramp and grey beyond 15 Å or where a
  subunit has no IP3. The distance is measured on the displayed coordinates,
  so it is valid in any numbering (rat 7LHF has no IP3 and comes out grey).
  "Show on structure" for the shell checks switches to this colouring. The
  exhibits plot JSD against distance (trend) and the per-shell means. The
  smoke test asserts at least 4 shell colours plus grey. Tests: 112 → 137.

**Measured.** Everything S22 publishes about the shells reproduces: the
pocket residue by residue, all 12 per-shell rows (n, means, medians, modal
fraction, whole-protein mean, p), and the three trends (ρ −0.175 / −0.436 /
−0.168). The new result concerns *where* conservation drops: in all three
paralogs the largest fall between adjacent shells is at 11.5 Å
(third→fourth: 0.023 / 0.045 / 0.019). Across 4.5 Å the change is
−0.007 / +0.003 / +0.002. So the "no step at 4.5 Å" claim holds, and more
strongly than a non-significant test of 12 against 14 residues could show
alone. The paper's FEL result puts its drop after the second shell, not the
third. That is a difference between two instruments, recorded in ROADMAP,
not a discrepancy with the paper.

**Where I was wrong first.**
- My first detail text for `P6.no_contact_step` stated "largest at 11.5 Å"
  as fixed prose. A verdict text must not assert a result, so it was
  removed; the found line reports the drop from the data.
- The first shell screenshot showed only one legend row. The panel was
  correct; the grab ran before Qt's deferred layout. The smoke test now
  processes events first.
- A stray `cat >` in a shell command waited on stdin and stalled one run.
  No project effect.

**Next:** Round 3, item 3 (Paper 2 tree viewer).
