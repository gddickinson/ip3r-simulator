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
