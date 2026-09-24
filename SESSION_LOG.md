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

## 2026-09-23 — Session 5: Round 3, item 3 — Paper 2 tree viewer

**Sync.** Both repos were up to date and `make sync-check` was clean, so no
verdicts moved.

**Built.**
- `run_app.command` at the project root, added at the user's request. It
  activates `ip3r_sim` and runs `python -m ip3r`, so it can be
  double-clicked in Finder, and arguments pass through
  (`./run_app.command checks`).
- `analysis/tree.py` asks Paper 2's clade questions of `rooted.nwk` with the
  existing Newick reader. Tip groups come from the census prefix on each
  label, and cyclostomes are recognised by genus, so no ip3r_genes clade
  table enters an answer. A paralog's whole clade is the MRCA of every tip
  whose record names it, and it must contain nothing foreign.
- `analysis/checks_tree.py` adds three rederived checks:
  `P2.paralog_clades`, `P2.cyclostome_lineages` and `P2.support_bar`. Each
  is calibrated by an input plant made in the tree itself (two tips
  swapped, or one support raised), not by editing the answer table.
- `analysis/tree_figure.py` draws the tree as a ladderized phylogram. It is
  shared by the new Tree tab (with a toolbar, tip-label and support
  toggles, a "Vertebrates" zoom and click-to-name) and by the exhibits of
  all four P2 tree checks. "Show" on a tree check opens the Tree tab.
- New parameters: `tree.alrt_min` 80 (Guindon 2010) and `tree.ufboot_min`
  95 (Hoang 2018). `test_tree.py` checks every definition on toy trees with
  known answers. The smoke test opens the tab from a check and asserts
  19/13/19 and cyclostome clades of 4 + 2. Tests: 137 → 148.

**Measured.** Everything the paper says about the tree reproduces from the
Newick alone. The whole clades are 19/13/19 at 100/100 and take in 4/2/1
unnamed shark, chimaera and coelacanth tips, the exact records listed in
`paralog_clades.tsv`. The six vertebrate tips outside all three clades are
all cyclostomes. They form two cyclostome-only clades, each holding hagfish
and lamprey. The 4-tip clade (99.5/100) branches first among the 57
vertebrates, and the 2-tip clade joins ITPR2+ITPR3 at 83.1/77. Of the 131
bipartitions, 107 clear SH-aLRT 80, 97 clear UFBoot 95 and 91 clear both.
A second rule for the whole clade ("expand the named core while the clade
stays pure") gives the same 19/13/19.

**Where I was wrong first.**
- My first support count was 92 of 132 nodes, where the paper has 91 of
  131. The paper was right. Rooting splits one unrooted edge into the
  root's two children, IQ-TREE labels both, and I had counted that one
  bipartition twice. `bipartitions()` now counts it once, and a toy-tree
  test pins that down.
- My first toy Newick in `test_tree.py` had unbalanced parentheses. It is
  now built from named sub-clades.
- In the first vertebrate-zoom screenshot, tip labels spilled above the
  axes and the clade labels were cut off. Labels are now clipped, and the
  zoom sets x as well as y.

**Next:** Round 3, item 4 (Paper 3/4: genome × paralog grid of the character
matrix and the recovery channel).

## 2026-09-23 — Session 6: Round 3, item 4 — Papers 3/4 genome × paralog grid

**Sync.** Both repos were up to date and `make sync-check` was clean, so no
verdicts moved.

**Built.**
- `analysis/genome_grid.py` joins the sweep's three per-cell tables into one
  grid of 309 genomes × ITPR1/2/3 + RyR. The layers are search grade, known
  genes missed, S15a state and protein-record recovery. The recovery channel
  is rebuilt from the table's count columns, not read from its label. The
  contiguity bar is a new registered parameter,
  `genomes.contiguity_bar_bp` = 142,212 (D4's median ITPR span, cited to
  ip3r_genes).
- `analysis/stats.py` gains `fisher_exact` (two-sided, own hypergeometric
  sum), `logistic_fit` (IRLS, Wald p) and `wilson`, and each is calibrated
  on a known answer.
- `analysis/checks_genomes.py` adds three rederived checks:
  - `P3.miss_by_contiguity`: medians, chromosome-level misses, the counts
    above and below the bar, and whether our bar split matches
    `spans_gene`.
  - `P3.contiguity_tests`: all 12 rows of `contiguity_tests.tsv`
    recomputed.
  - `P4.recovery_channels`: per-cell counts against
    `gene_recovery_by_cell.tsv`, and the reason split against Paper 4's
    ledger AR12–AR16.

  Each is calibrated by an input plant: one genome moved just under the
  bar, one RyR cell missed, one reachable gene's resolving count zeroed.
- `analysis/grid_figure.py` draws the grid and the three exhibits.
  `ui/genomes_panel.py` is the new Genomes tab, and "Show" on any P3/P4
  check opens it on the right layer. The smoke test opens it from
  `P3.miss_by_contiguity` and asserts 309 genomes, 182 misses and the bar at
  row 189, then switches to the recovery layer (292 reachable). Tests:
  148 → 164.

**Why.** Paper 3's no-loss result rests on its method control: the search
misses genes only where assemblies are fragmented. Paper 4's headline is
that most genes are unreachable from protein records. Both are claims about
*which* genomes, and until now this project had only checked their totals.
The grid shows the misses sitting below the bar row by row. The checks now
recompute each statistic and do not only count it.

**Measured.** Everything reproduces:
- Misses are 140/923 (ITPR) and 42/309 (RyR), Fisher p 0.578, odds 1.14.
- Odds of finding the gene rise 8.10× (β 2.092, p 1.29e-22) and 19.99×
  (β 2.995) per tenfold N50.
- The Mann–Whitney U is 99,201 / 10,475. Chromosome-level assemblies miss
  3/512 and 0/172.
- Above the bar: 189 genomes and 5/563 misses. Our split and `spans_gene`
  agree in all 1,236 cells.
- Recovery: 257/309, 260/307, 227/307 and RyR 196/309. The reasons are
  289/186/254/15, with 179 reachable. The rebuilt channel matches the label
  in 1,236/1,236 cells.

39 checks: 37 confirmed, 2 discrepancies (unchanged).

**Where I was wrong first.**
- My Mann–Whitney p came out exactly twice the published one (1.05e-52 vs
  5.26e-53). The published test is one-sided (found N50 greater than
  missed), and I had doubled it. It was the checker's error, not the
  paper's.
- The Paper 4 ledger's `expected` column is not all integers, so my first
  version crashed parsing unrelated rows. It now reads only AR12–AR16.
- My first tolerances on the odds ratios were a bare 1 %, which is not a
  registered number. Values are now compared at the precision the table
  prints them.

**Next:** Round 3, item 5 (Paper 1: presence/absence across eukaryotic
clades).

## 2026-09-23 — Session 7: Round 3, item 5 — Paper 1 range across eukaryotes

**What.**
- `analysis/range_table.py` joins the S20 proteome sweep (presence,
  taxonomy, the six assignment tables, relaxed hits) and the S23 genome
  sweep (manifest, control ledger, copies, copy-number ledger).
- `analysis/checks_range.py` holds seven rederived checks, each calibrated
  by an input plant. `P1.absences` moves there from `checks_evolution.py`
  and is upgraded: it used to count the summary table's columns, and now it
  rebuilds every row of `absence_at_genome.tsv` from the per-genome
  ledgers.
- `analysis/range_figure.py` draws the clade bars and five exhibits.
  `ui/range_panel.py` is the new Range tab, and "Show" on any P1 check
  opens it. The smoke test asserts 70 rows, 45 with a call, 662/6,928
  proteomes and 35/35 absences held.
- Six registered parameters (`range.*`), all cited to the ip3r_genes
  thresholds they reproduce. Tests: 164 → 170.

**Why.** Paper 1 is the one the rest of the series takes its family call
from. Until now this project checked a single number of it (the 35
absences), and only by reading a summary. The paper's structure is
"absent here, present beside it in the same kingdom, and the instrument
demonstrably works in both". That can only be seen per clade, and checked
from the rows.

**Measured.** Everything reproduces:
- 662/6,928 proteomes and 45/135 clades. Presence from 2,012 ITPR
  assignment records agrees with the presence table in 6,854/6,854 taxa.
- All 13 named phylum counts match. All 24 cells of the relaxed table
  match, and the one non-MIR exception (PF08454, *Oryza barthii*) clears
  neither PF08709 nor a full-length profile.
- G3, rebuilt from its sentence (≥ 10 proteomes, no call, eukaryotic
  phylum or class), gives the same 35 targets. All 35 genome rows rebuild
  from the ledgers, and neither weakly controlled ciliate genome
  (*Tetrahymena*, *Ichthyophthirius*) sits in an absence target.
- Copies: 117 genomes with none and 43 with one; the top three are
  Macrostomum 18, Stentor 13 and Dysidea 8; Cymbomonas carries 3.
- Chase: 47 real genes and 52 fragments; no contaminant; identity 19.9–39.9
  % (plants) and 20.5–33.5 % (fungi); maximum 45.8 %.

46 checks: 44 confirmed, 2 discrepancies (unchanged). No verdict about the
publication moved.

**Where I was wrong first.** `P1.record_chase` first found 50 real genes
against the published 47. Rule R5 also calls a record a fragment when
UniProt flags it, at any length. Three flagged records are 2,366–2,858 aa,
so my length-only rule passed them. It was the checker's error, and the
fourth time in this project that a discrepancy was the checker's.

**Noted about ip3r_genes (not a discrepancy).** Paper 1 ledger row C20
("Macrostomum lignano carries 18 IP3 receptor genes") verifies
`copy_number.18 == 1`, i.e. that *some* genome carries 18. It does not
verify which genome. `P1.copy_number` checks the species, and it is
Macrostomum.

**Next:** Round 3's scheduled items are done. Round 4 (Mak 1998 gating
model) is next, unless an emergent Round 3 item is preferred.

## 2026-09-24 — Session 8: Round 4, item 1 — Mak et al. 1998 gating model

**What.**
- `physics/gating_mak.py` implements the biphasic Hill model: Eq. 1 has
  one denominator (not a product of two Hill factors), and in Eq. 2 K_inh
  depends on IP3.
- `physics/bell.py` is a model-agnostic ruler: peak, and the Ca²⁺ at half
  that peak on each flank. The DYK `bell_peak` and `hill_fit_left_flank`
  now use it, so both models are measured the same way.
- Nine registered parameters: seven `mak.*` constants from the paper and
  two method parameters for the comparison span (33 nM, 10 µM). There is a
  new reference, `mak1998`.
- The Gating panel has a model selector, and its second plot shows K_inh
  against IP3 with K_act dotted. `python -m ip3r gating --model mak` prints
  the same.
- The smoke test switches the panel to Mak, asserts on it, and writes
  `docs/img/gui_gating_mak.png`. Tests: 170 → 178.

**Why.** Up to now the app only said "experiment finds inhibition alone is
tuned; DYK doesn't". With this round the app shows that experimental
result and measures the difference with an instrument that can fail.

**Verified against the paper (PMC28128 text).** P_max 0.81, K_act
210 ± 20 nM, H_act 1.9 ± 0.3, H_inh 3.9 ± 0.7, K_∞ 52 ± 4 µM, K_IP3
50 ± 4 nM, H_IP3 4 ± 0.5. Uncertainties are kept in the source notes.
Eq. 1's form was confirmed separately: a single denominator 1 + a + b.

**Measured.**

| model | half-activation (33 nM → 10 µM IP3) | half-inhibition |
|---|---|---|
| De Young–Keizer | 2.01× | 2.76× |
| Mak 1998 | 1.016× | 6.22× |

- The plateau at 10 µM IP3 is 0.77–0.81 over 1–20 µM Ca²⁺ (paper: ≈0.8).
- Below K_IP3 the bell collapses: at 10 nM the peak is 0.11 and K_inh is
  0.08 µM, below K_act. The paper reports this for 10–20 nM.
- The Hill curve's K_inh at 33 nM is 8.3 µM; the paper's measured point is
  9.5 µM.
- Over 0.1 → 10 µM IP3 (the old DYK comparison span), Mak barely moves
  (1.06×). Its IP3 sensitivity is about tenfold higher than DYK's, so the
  span was taken from the paper's own statement: at 33 nM, activation was
  not affected.

**Correction to Session 7.** There are 45 checks (43 confirmed,
2 discrepancies), not 46/44. `P1.absences` moved between modules and was
counted as new. The README said 39; it now says 45. No verdict moved.

**Next:** Round 4, item 2 (a puff model with low resting activity).

## 2026-09-24 — Session 9: Round 4, item 2 — park/drive puffs

**What.**
- `physics/park_drive.py` implements the Siekmann et al. 2012 six-state
  IP3R-1 model with the gating variables of Cao et al. 2013.
  `stationary` uses detailed balance on the tree, and the test proves it
  equals the null space of the full generator.
- `physics/puffs_pd.py` is the same mean-field cluster as `puffs`, with the
  receptor swapped. An open receptor also sees its own mouth (+120 µM), as
  in Cao's scheme.
- `physics/puff_compare.py` is one ruler for both clusters:
  - `recruitment` gives the Fano factor, the open fraction, and blip /
    multi-channel / half-cluster events;
  - `coupling_scan` runs both receptors over the same couplings.
  `coupling_effect` moved here from `puffs`, with a `model` argument.
- Both simulators now record, per 1 ms bin (`puff.record_dt`), the snapshot
  and the most open at once. Events are read from the peak. The DYK
  snapshots are unchanged, so its Fano numbers are too.
- `ui/puffs_panel.py` (moved out of `dynamics_panel.py`) has a receptor
  choice, an event-size histogram and "Scan coupling". The CLI gained
  `puffs --model park-drive` and `puffs --scan`. The smoke test drives
  park/drive and writes `docs/img/gui_puffs_pd.png`.
- 48 new parameters (41 `pd.*`, 7 `puff.*`) and three references. The park/
  drive entries live in `scripts/parameter_table_pd.py`, and
  `scripts/param_entry.py` holds the shared constructor. Tests: 178 → 201.

**Why.** The DYK cluster could only say "coupling raises the Fano factor
a little". The roadmap asked for a receptor whose resting activity is low
enough that blips and puffs separate. The comparison is built so that only
the receptor differs: same cluster, coupling scheme, seed and ruler.

**Where the constants came from.** The rates are in Cao 2013's Table S1,
behind PMC's browser challenge, so they were not reachable. Cao et al. 2014
(PLoS Comput Biol, open access) publishes its model code as Text S1, and
every constant was read from that code. Cao 2013's main text was read for
the equations. Its Eq. 10 writes λ_h42 as a Ca²⁺ switch with V = 100 s⁻¹,
while the 2014 code switches on the open state at 20 s⁻¹ with 0.5 s⁻¹
recovery. The code's values are used and the difference is in the source
notes. The main text says "λh24" where Eq. 10 defines λh42, a typo in the
paper. Two sanity checks against the text pass: the drive-mode P_open is
0.701 ("around 70 %"), and q54/q45 is 303 ("∼300 times").

**Measured** (30 s, 0.2 µM IP3, seed 0; large = peak ≥ 10 of 20):

| coupling | DYK Fano | DYK large | PD Fano | PD large |
|---|---|---|---|---|
| 0 | 0.93 | 0 | 0.98 | 0 |
| 0.09 | 1.32 | 1 | 2.79 | 10 |
| 0.17 | 1.21 | 0 | 2.78 | 17 |
| 0.32 | 1.14 | 0 | 2.23 | 16 |
| 1.08 | 1.25 | 1 | 1.60 | 6 |
| 2.00 | 1.19 | 1 | 1.61 | 7 |

- On seeds 1–3 at 0.1 µM coupling, park/drive has 12, 12 and 14 large
  events; DYK has 0–1 at 0.1 and 0 at its own 1 µM.
- Uncoupled open fraction: DYK 5.2 %, park/drive 1.7 %. The stationary
  formulas give 1.6 % against 1.1 %, so the stochastic DYK subunit cluster
  runs well above its own Li–Rinzel steady state.
- The event-size histogram at 0.17 µM has a trough at 6–7 channels and a
  second mode at 8–11. DYK's decays monotonically at every coupling.
- The park/drive coupling (0.1 µM) was chosen by the rule that chose DYK's
  1 µM (maximum Fano). Cao's microdomain gives ≈ 0.11 µM per open channel.

**The instrument, calibrated.** The first "must fail" case I wrote (a 2 ms
step) did not fail. The split step is more accurate than I assumed. Measured
against the clamped stationary P_open at 0.5 µM Ca²⁺, the error is 0.2 % at
0.1 ms, 1.6 % at 0.5 ms, 4.5 % at 2 ms and 10.5 % at 5 ms. The test now
uses 400 channels and a 3 % tolerance, and the must-fail case uses 5 ms.

**Incidental.** `mak1998` was declared twice in `reference_table.py`
(Session 8). Nothing checked reference keys for duplicates; the build now
does.

**Not changed.** No `ip3r_genes` table moved (`make sync-check` clean). The
45 checks are untouched.

**Next:** Round 4, item 3 (unitary current from the pore profile).

## 2026-09-24 — Session 10: Round 4, item 3 — unitary conductance

**What.** Each ITPR3 deposit's pore profile now gives a K+ conductance in
symmetric 140 mM KCl.
- `physics/permeation.py` and `physics/_pnp_kernels.py` are PIEZO1's 1-D
  drift-diffusion solver, ported. The kernels are unchanged; the wetting
  gate and sodium protocol are dropped.
- `physics/pore_charge.py` builds the wall charge from the deposit's own
  side chains. PIEZO1's open structure had no side chains, so it placed
  charges at the Cα. Here stubbed residues are counted as unplaced, never
  guessed.
- `physics/unitary.py` measures one deposit three ways (series, neutral,
  charged) and runs the state panel with a sensitivity sweep.
- New: the CLI `unitary`, a Channel-tab button with a bar plot against the
  measurements (`docs/img/gui_unitary.png`, a smoke-test step), 19
  parameters (`scripts/parameter_table_perm.py`), 3 references, and 7 tests
  (201 → 208 passed, plus the smoke test).

**Why.** The state panel said where the gate opens. This asks whether the
opening is enough, which is a number the single-channel literature has
measured.

**The measurements, read before registering.** Vais et al. 2010
(PMC2995152) gives 545 ± 7 pS for rat ITPR3 in DT40 nuclei (symmetric
140 mM KCl, 0 Mg²⁺, room temperature). Vais cites Mak et al. 2000 as
"370 ± 8 pS", but Mak 2000's own text (PMC2217211) says **358 ± 8 pS**. The
primary value is registered. I also first gave vais2010 a wrong title from
memory; it was corrected from Crossref before the build.

**Measured.**
- Shut (r_free below the K+ radius 1.38 Å): 8TKH 0.25, 8TLA 0.74,
  7T3P 0.83, 6DQN 0.85, 6DQJ 0.99, 8TKG 1.03 Å.
- 8TKF: series 64.2, neutral 65.3, charged 33.0 pS. Over diffusivity
  0.25–1× and ion radius 1–2 Å the neutral range is 25–150 pS and the
  charged 13–75 pS. Measured is 358/545 pS, so the model is 2.4× short at
  best. Our PIEZO1 model was 1.5× high by the same method, so this is not
  a systematic undercount by the method.
- Where the resistance sits: along the ~50 Å pore, with a third at the
  filter slice (r_free 3.08 Å).
- Charge rings (4 copies each): E2398, K2482, D2478, D2518, D2522, K2529.
  The same rings, plus R2524/E2532, line every state.
- Removing one ring at a time gives: E2398 15.5, K2482 63.6, D2478 23.1,
  D2518 29.9, D2522 30.6, K2529 48.9 pS. Acidic only: 174 pS.
- Sensitivity of the charged value: lining margin 0/1.5/3/5/8 Å →
  58/33/33/7/8 pS; smoothing 1.5/3/6 Å → 16/33/151 pS.
- Voltage: 25.6 pS at −20 mV, 31.7 at +20 mV (rectifying; measured I–V is
  linear).
- Grid: step 1/0.5/0.25 Å → neutral 69.5/65.3/65.4, charged 35.1/33.0/31.7.

**Interpretation.** A positive ring next to a negative one is a junction.
Each carrier must cross the zone where it is the excluded co-ion, so small-
signal conductance falls. The Donnan and junction tests show the solver
does this correctly. Whether the protein does is another matter: D2478 is
salt-bridged (2.5 Å) to R2471 of the neighbouring subunit, which sits
outside the lumen and is not counted. So the charged number should not be
quoted; the neutral bound should. The finding pinned by a test: no
setting of the unmeasured constants brings the only open deposit to
either measured value.

**Not changed.** `make sync-check` was clean; no `ip3r_genes` table moved,
and the 45 checks are untouched.

**Next:** Round 5, item 1 (parameter editing), or the salt-bridge / RyR
emergent items under Round 4.



## 2026-09-24 — Round 5.1: parameter editing in the GUI

**What.** Help → Parameters (`Ctrl+Shift+P`) is now an editor, ported from
PIEZO1's. It shows value, default, unit, bounds, kind and source, with the
full reference on the tooltip, and offers a filter, "Show only modified",
edit with clamp reporting, reset selected/all, and import/export. An amber
banner (`ui/params_banner.py`) runs across the window top whenever any value
differs from its default, including values set by `IP3R_PARAMETERS` at
start-up. The registry gained change listeners (`subscribe`), a filter
(`matches`) and `write_overrides`/`read_overrides`; import replaces the
overrides rather than merging, and refuses a malformed file without touching
anything. 7 tests (208 → 215), plus a smoke-test step and two screenshots.

**Why the listeners matter beyond the banner.** Checks refuse to *run* when
the registry is modified, but two memo caches they read
(`checks_structure._summary` and `shell_constraint.measured_shells`/`pocket`)
could be filled by an edited computation and then served after "Reset all".
A check would then confirm on numbers made with the edited value. Those
caches now clear on every change. The test proves the hazard is real: it
fails with the subscription removed. `pairwise._cached` was already keyed on
the gap costs, so it needed nothing.

**Design choices.** Nothing typed in the dialog is persisted: the registry
docstring already said a value typed once must not silently change what a
later run computes. Export is the route to reproduce a set. The first banner
lived above the viewport and wrapped into a tall column in the narrow
centre, so it moved to a full-width toolbar holding one line, with the
detail on the tooltip.

**Not changed.** `make sync-check` was clean; no `ip3r_genes` table moved,
and no check verdict changed.

**Next:** Round 5, item 2 (session save/restore).


## 2026-09-24 — Round 5.2: session save/restore

**What.** File → Save session… (`Ctrl+Shift+S`) and Open session… (`Ctrl+O`)
write and read the view as JSON, and `python -m ip3r --session FILE` starts
on one. The headless half is `io/session.py` (`Session`, `save_session`,
`load_session`, `parameter_differences`), ported from PIEZO1's. The GUI half
is `ui/session_controller.py`. The registry gained `replace` (the whole
override set, one notification); `read_overrides` now uses it.
21 tests (215 → 236), and three smoke-test steps.

**What a session holds, and why.** It holds the view: deposit, style,
colour, layer, subunits, sites, pore, camera (rotation, pivot, distance, pan,
projection), tab, and the transition *spec* (end, fit, method, frame,
paint). It holds no coordinates and no results. As in PIEZO1, a file with
its own copy of the numbers would drift from the code. The field set is
pinned by a test so adding one is a decision. Unlike PIEZO1, it also records
the **parameter overrides**. Every number the view shows depends on them, so
reopening a session under a different set would show different numbers
under the same name. Restoring compares the two sets. If they differ, the
user chooses to apply the session's set or keep their own, and applying
always re-measures the deposit, even the one on screen. Nothing is applied
silently, which is the registry's rule.

**Asynchrony.** Both the load and the morph run on workers, so a restore is
held as pending and finished from `_loaded` and `_transition_built`. A
restore is dropped (and says so) if another deposit arrives first, so style
and camera are never applied to the wrong structure. The atom count is
stored, and a changed file is reported because the camera addresses
coordinates.

**Calibration.** The smoke test saves 8TKG→8TKF at frame 5 (backbone, gate
lining, chain D hidden, orbited camera, Channel tab, one parameter edited).
It then resets parameters, loads 6DQN, restores, and compares every field.
It passed first time, so a planted fault (camera restore removed) was run
and caught (rotation and distance reported). Looking at the screenshot
found a real bug the field comparison could not see: the deposition list
still highlighted 6DQN. The restore now loads through the list, and the
smoke test checks the list selection.

**Not changed.** `make sync-check` was clean; no `ip3r_genes` table moved,
and no check verdict changed.

**Next:** Round 5, item 3 (variants as spheres; VUS by conservation layer).


## 2026-09-24 — Round 5.3: variants as spheres; VUS by conservation layer

**What.** The Variants tab can draw the S17 harvest on the structure. "Draw
on structure" puts one sphere per variant residue on the Cα of every visible
subunit, in class colours. Where a residue carries alleles of several
classes, the most decisive class wins. With a layer chosen under "VUS by
layer", each VUS is placed against its own gene's labelled medians (Paper 5
§8). The table gains a stratum column, a strip plot shows the three classes
and both medians on a fixed 0–1 axis, and the VUS spheres take the stratum
colour (unscored = grey). The tab now follows the deposit's paralog on load.
New modules: `analysis/vus_strata.py` (the rule, headless),
`analysis/vus_figure.py`, `analysis/checks_variants.py` and
`render/variant_spheres.py`. There are 6 new tests (236 → 242) and one smoke
step.

**The check.** `P5.vus_stratification` (rederived) rebuilds all 12 rows of
`vus_stratification.tsv` from `variants.tsv` and the per-residue tables,
using a rule written from S17's description: positions not alleles, the
occupancy floor, and ties counting. Every field agrees. The planted relabel
of one ITPR3 P/LP variant in `variants.tsv` flips it, which proves the rows
are recomputed and not just read.

**Two things checked because they could have been checker bugs, and were
not.** (1) Our `constraint.json` masks the deep layer by `deep_reliable`,
but S17's stratification masks by occupancy ≥ 0.5. The two agree at every
residue of all three paralogs, so the viewer (which uses the resource) and
the check (which uses the tables) give the same numbers, and a test pins
this. (2) The stratification takes every source, so the 13 curated UniProt
P/LP records are in its P/LP median, while the AUC test is ClinVar only. Run
ClinVar-only, the table is identical, so the inconsistency moves no number.
This is recorded in `docs/SCIENCE.md` and is not a discrepancy.

**The spheres follow the view, not the deposit.** The first version placed
them at deposited coordinates, which is wrong mid-morph. The smoke test runs
on the restored 8TKG→8TKF session at frame 5 with chain D hidden. It
requires the spheres on exactly the visible chains (930 spheres, 3 subunits)
and displaced from the deposit, and it requires a refusal when the tab is
switched to ITPR1.

**Record corrected.** HEAD had 45 registered checks, not the 46 the Round 3
(Paper 1) entry stated. The count is now 46: 44 confirmed, 2 discrepancies,
both unchanged.

**Not changed.** `make sync-check` was clean; no `ip3r_genes` table moved,
and no check verdict changed.

**Next:** Round 5, item 4 (AlphaFold models for the unresolved stretches).


## 2026-09-24 — Round 5.4: AlphaFold models for the unresolved stretches, seams shown

**What.** Representation → Completeness draws an AlphaFold model's residues
where a deposit has none. New modules: `io/predictions.py` (AlphaFold DB
download, with entry and version discovered from the API),
`structure/graft.py` (the fill), `structure/graft_calibration.py` (how good
a fill is), `ui/fill_overlay.py` and `ui/fill_controller.py`. There is a
`graft` CLI command, `fetch` also downloads the models, and there are 7 new
registered parameters and 2 references (Jumper 2021, Varadi 2022). A session
gains a `completeness` field (the choice, not the fill). There are 12 new
tests (242 → 254) and two smoke steps.

**What AlphaFold DB actually holds** (the first thing measured). Canonical
ITPR3; ITPR1 isoform 4 only; rat ITPR1 isoform 8 only; for ITPR2 a
181-residue isoform. So the model is never matched by name. It is chosen by
identity by number against the variant-painting bar. All seven ITPR3
deposits match at 99.1–100 %. 7LHF (22.1 %) and 9YKK (17.7 %) are refused,
with every model's identity in the message.

**Method: port of PIEZO1's local-anchor graft, per stretch.** Anchors must be
resolved, unstubbed and the same amino acid in both models. That makes 6DQN
skip 8 stretches per subunit around its unregistered 1434–1546 segment
instead of fitting them to it. A terminus is extrapolated and filled only in
"+ gaps and ends". The fill is a *separate* structure. Nothing measures on
it, and `place()` re-fits it on its anchors in any coordinates, so it follows
morph and mode frames.

**Calibration, and the parameter it moved.** Stretches that another ITPR3
deposit leaves unresolved, but 8TKG or 8TKF resolves, were hidden, filled and
scored. Over 31, the median was 1.46 Å against 5.59 Å for a straight line
and 10.5 Å for a whole-chain fit. The fill beats the line 28/31, and pLDDT
tracks error (ρ −0.79). The true seams reached 5.42 Å, so the first join
tolerance (4.5 Å) would have called 9 of 62 correct seams broken. It is now
5.5 Å, cited as measured here. The honest limit: those stretches sit at
pLDDT 53–77, while 8TKG's real gaps average pLDDT 38. The calibration does
not reach them, and the docs say so.

**Found by the smoke test.** On its first run the fill for a deposit showing
morph frame 5 was drawn at the deposited coordinates, because `show()` did
not place it. Fixed; the step now requires the fill to move with the frame
and to return to its built position at frame 0. A hidden number was also
removed before it shipped: a minimum chain length (`anchor_window × 10`) was
replaced by the rule already in use, that a chain must itself be in the
prediction's numbering.

**Not changed.** `make sync-check` was clean; no `ip3r_genes` table moved,
and no check verdict changed (46 checks: 44 confirmed, 2 discrepancies).

**Next:** Round 5 is complete. The emergent items are listed in `ROADMAP.md`.



## 2026-09-24 — Round 4 emergent: salt bridges cancelled before the wall charge counts

**What.** New `physics/salt_bridges.py`. `pore_charge(..., pair_bridges=True)`
drops every lining group that is half of an ion pair. `Unitary` gains a
fourth reading, *paired*, which appears in `python -m ip3r unitary`, in the
Channel panel (a third bar) and in the sensitivity sweep. One new registered
parameter (`pore_charge.salt_bridge_cutoff`, 4 Å) and one reference (Barlow &
Thornton 1983, whose abstract states the ≤ 4 Å criterion; read on PubMed
before registering). There are 4 new tests (254 → 258), and the smoke test
now requires the paired reading to be drawn.

**Why this rule.** A carboxylate held by a guanidinium is not a free charge
acting on the lumen, so the round asked for ion pairs to be neutralised
before any charged number is trusted. Pairs are matched one-to-one, closest
first, so an Arg reaching two Asps cancels one: every bridge removes exactly
+1 and −1. Partners are searched in the whole deposit, because D2478's
partner (R2471 of the next subunit) does not itself line the pore. The tests
were checked by planting many-to-one pairing, which the first test catches.

**Measured on 8TKF.** At 4 Å only the four D2478–R2471′ bridges qualify
(2.45–2.58 Å). The wall goes from −8 to −4 e and the conductance from 33 to
23 pS (8–54 pS over diffusivity × ion radius). The cutoff matters only near
4.25–4.5 Å: D2518–R2524′ sits at 4.23–4.34 Å and K2482–D2400 at 4.31–4.40 Å.
Over cutoffs of 3–6 Å the reading is 23 / 23 / 23 / 23 / 38 / 38 / 38 pS
(3, 3.5, 4, 4.25, 4.5, 5, 6 Å), always below the neutral 65 pS. The 18.6 M
peak partition density does not move, because it sits on the D2518/D2522
rings at z ≈ −63 Å (lumen 4.4 Å), not at the filter.

**What it means.** Cancelling ion pairs does not rescue the charged model,
and the gap to 358–545 pS is not a charge-counting artefact. The least
plausible assumption left is full ionisation of eight aspartates in a 4.4 Å
lumen, a pKa question, now an emergent item.

**Not changed.** `make sync-check` was clean; no `ip3r_genes` table moved,
and no check verdict changed (46 checks: 44 confirmed, 2 discrepancies).

**Next:** the gate radius along the morph (Round 2), or the uncharged 2.4×
shortfall against RyR1 open deposits (Round 4).


## 2026-09-24 — Round 6.1: RyR1 structures, and the test that refuted this morning's pairing

**What was asked.** The user asked for the app to model RyRs. Asked about
scope, they chose structures plus gating kinetics, RyR1 only. This is the
structural half. Kinetics is Round 6.2.

**What.** RyR1 is a *numbering* (`config.RYR_ACC`, `NUMBERINGS`), not a
publication paralog, so every `ip3r_genes` check and import is untouched.
`scripts/curate_ryr.py` (`make ryr`) writes `resources/ryr1.json`: the P11716
sequence, InterPro Pfam domains, and a six-deposit panel. It records the
SHA-256 of every response and, for each of the 140 rejected entries, the
first rule it fails. The registry, annotations (grey where the publication
has no data), numbering, state panel, unitary (bath by family), transition
preset and Channel panel all follow the loaded deposit's family. New pieces:
`physics/ryr_mutants.py`, the `mutants` CLI, `parameter_table_ryr.py` (the
bath and six measured conductances, read from Xu et al. 2006 Table 2 and
Methods before registering) and `tests/test_ryr.py` (7 tests; 265 total).
The smoke test gains two RyR steps.

**Why these deposits.** The PDB holds 158 RyR1 entries, so the panel had to
come from rules, not preference. The rules: full-length EM tetramers, wild
type, activators only, ≤ 4 Å (side chains carry the wall charge), and the
best resolution per title-stated state. They picked the best primed deposit
(8RRX, nanodisc) and the best open one (9HEO, micelles) from different
preparations, so a sixth rule adds the open deposit's primed partner from
the same paper (9R8O) for the morph.

**Measured.** All six deposits are 100 % in P11716 numbering. The shut states
gate at I4937, and only 9HEO opens (5.05 Å). The morph 9R8O → 9HEO has
RMSD 2.72 Å, and the ANM overlap is 0.174 over 20 modes (null 0.031),
against ITPR3's 0.64. The 9HEO conductance is 136 pS neutral and 180 pS
charged against 801 pS measured: 5.9× short. So the 2.4× ITPR3 shortfall
is the method's, not ITPR3's.

**The finding.** On 9HEO, Xu et al.'s five charge-neutralising mutants give
the right direction for all four lining residues and the right null for
E4955Q. But D4899Q (measured 0.20×) is modelled at 0.90×, and at 1.00× with
salt bridges cancelled. D4899 is the homologue of ITPR3's D2478 (own
alignment; R4892 ≡ R2471 too) and is bridged the same way, so the "paired"
reading added this morning is refuted by measurement. Its ROADMAP and
SCIENCE entries now say so, and *charged* is the reading to trust. Found on
the way: the charged solve for primed 9R8O did not converge but printed
0.0 pS. It now prints `n.c.`, and the sweep skips unconverged corners.

**Not changed.** `make sync-check` was clean. The 46 checks give the same
verdicts (44 confirmed, 2 discrepancies).

**Next:** Round 6.2, RyR1 gating kinetics and sparks.


## 2026-09-24 — Round 6.2: RyR1 gating and sparks

**What.** `physics/ryr_gating.py` (Stern 1997 four-state scheme and
Murayama 2015 bell), `physics/sparks.py` (RyR1 cluster), a third model in
Dynamics → Gating and a third receptor in Puffs, `cli_ryr.py` (`mutants` moved
there to keep `cli.py` well under budget, plus `ryr-gating` and `sparks`),
20 registered constants and 3 references, `tests/test_ryr_gating.py`
(13 tests; 278 total), and two smoke steps.

**How the sources were chosen.** A research agent searched for RyR1 models
whose constants can be read from the source's own text. JBC, Biophys J
PDFs and several classics (Meissner 1997, Keizer–Levine, Copello 1997) are
behind a script challenge. I re-checked every constant used against the
downloaded text myself. Two corrections to the agent's report came out of
that. "180 nM" in Stern's Table I is fura-2's Kd, not resting Ca²⁺. And
"2N = 60" counts the voltage-coupled channels too, so the Ca²⁺-gated cluster
is 30.

**Decisions and why.**
- Stern's k_i is printed as 2 × 10⁻⁶ M⁻¹ s⁻¹. The text's 10 µM Kd makes it
  2 × 10⁶, so it is registered that way, with both the printed and the
  corrected value in the source note. A test shows the printed value would
  remove inhibition.
- The coupling is derived from Stern's current, diffusivity and spacing,
  not chosen, and a scan shows the answer holds down to about 0.2× of it.
- `spark.dt` was first 1e-4 s. Measuring convergence showed 15–20 % longer
  sparks than at 1e-5 or 2.5e-5 s, so it is 2.5e-5 s.

**Measured.** The bells agree on half-activation (3.9 vs 4.4 µM). On
half-inhibition the scheme is 6.7× too sensitive (48 vs 320 µM). Sparks:
uncoupled, only blips; coupled, 1.5 per second reaching 25–30 of 30, Fano 4,
switched on between 0.44 and 0.80 µM per open channel. They last ~120 ms
against a measured release of 6.3 ms (frog). The trace and a direct
calculation show why: the mean-field cluster has a self-sustaining point
between 5 and 6 open. So the next physics step is a spatial Ca²⁺ field, not
a new constant.

**Found on the way: the GUI smoke test outgrew its budget.** It timed out
at 540 s. Per-step timestamps (now printed) showed the RyR1 load as the
largest step, and a profile put 14 of 15 s of a cartoon build in
`parallel_transport_frames`, where `np.cross` ran on 436,000 single
3-vectors. The loop is now plain-float arithmetic: identical to 4e-15 on
random paths, 29× faster, and it speeds every structure. The smoke test
dropped from 771 to 344 s, and its budget is 600 s. There was no test of
the frames, so `tests/test_spline.py` adds three (orthonormal, a plane curve
keeps its out-of-plane normal, a helix normal never turns faster than the
tangent). Profiling also exposed a latent crash: ticking the Paper 6 modules
on a RyR deposit raised `ModuleRefusal` inside a Qt slot, and it now draws
nothing. 281 tests.

**Not changed.** No `ip3r_genes` table moved, and the 46 checks are
unchanged (44 confirmed, 2 discrepancies).

**Next:** a spatial Ca²⁺ field for the clusters (it would test the IP3R
puffs too), or the gate radius along the morph.


## 2026-09-24 — Round 6.3: RyR1 sparks in the junctional cleft

**What.** `physics/cleft.py` (Stern 1997's cleft, a steady finite-volume
solve giving each C channel's Ca²⁺ per open C channel), and
`physics/sparks_cleft.py` (the array simulated exactly by Gillespie's
method). The cleft is a fourth receptor in Puffs and `sparks --cleft` in
the CLI. `spark_ends` measures each spark's duration and how many channels
are inactivated at its start and end; both spark simulators now record
`n_inactivated`. There are five registered geometry parameters,
`tests/test_cleft.py` (9 tests; 290 total) and a smoke step.

**Why this geometry.** The roadmap asked for a spatial field. Rather than
invent a cluster layout, I read Stern et al.'s Appendix and ported their
method: Table I dimensions, the Fig. 7 B chessboard, the 30 nm diffuse
source and their Eq. 13 edge coefficient. Eq. 13 is an image on PMC and was
read from the GIF. Fig. 9 was digitised (the scale bar is 92 px = 500 µM)
to give the solve something to be held to.

**Decisions and why.**
- Gillespie rather than a fixed step, because the steady field makes the
  array a true Markov process (as Stern argued), so there is no step to
  converge. It is also fast: 0.1 s for 20 s simulated.
- The GUI's coupling knob scales only the off-diagonal. A channel's own
  release belongs to the channel, so "uncoupled" still means independent
  channels that see themselves.
- No buffer. A steady field is unchanged by fixed fast buffers. A mobile
  buffer needs time-dependent diffusion, which Stern also treated
  separately, so it is recorded as emergent rather than faked with the
  excess-buffer approximation.
- `docs/SCIENCE.md` had reached 605 lines, so the RyR1 sections moved to
  `docs/SCIENCE_RYR.md` to meet the 500-line rule.

**Measured.** Per pA the solve gives 171/65/56/15 µM (at the source, across
the row, 30 and 60 nm along it) against ~185/75/73/25 from Fig. 9. The
shapes agree; this solve falls off faster. Sparks last a median 21 ms
against ~130 ms mean-field, and they end with 16 of 30 channels
inactivated (0 at the start). The mean-field sparks end with 27 of 30
inactivated. The duration holds at 16–22 ms over coupling 0.75–2× and
source diameter 15–45 nm. At 0.5× coupling there are no sparks. The result
is still 3× the 6.3 ms frog release. A planted wrong transition table
fails two of the new tests (checked, then restored).

**Not changed.** `make sync-check` was clean, no checks were touched, and
the verdicts are unchanged (44 confirmed, 2 discrepancies).

**Next:** refit Stern's inactivation to Murayama's bell and rerun the
cleft. The geometry cannot close the remaining 3×.


## 2026-09-24 — Round 6.4: what ends a cleft spark

**What.** `ryr_gating.fit_to_bell` fits Stern's Ka and Ki to both half-peak
flanks of Murayama's measured bell, moving the off rates and keeping the on
rates. `ryr_gating.with_constants` builds a scheme from Ka, Ki and an
inactivation rate scale. `physics/spark_termination.py` runs the cleft
array under Stern, the 25 °C fit and the 37 °C fit, and scans Ki and the
inactivation rate. `spark_ends` now flags an `unterminated` spark (still
running at the trace's end). The CLI gains `spark-termination` and
`sparks --cleft --fit`, and the GUI gains a fifth Puffs receptor and the
fitted bell (dashed) in Gating. 8 new parameters: Murayama's 37 °C WT row,
read from the S1 Table .doc via the PLoS supplement URL, and the scan
grids. `tests/test_spark_termination.py` has 9 tests (299 total). A planted
fit that returns Stern unchanged fails 6 of them (checked, then restored).

**Why a scan and not just a refit.** A steady-state bell fixes only Ka and
Ki. How fast inactivation is stays free, and a refit at one arbitrary rate
could have hidden a rate that works. So Ki and the rate were both scanned.

**Measured** (4 seeds × 10 s). Stern: median 19 ms, all ended. Fitted to
25 °C (Ka 4.9, Ki 249 µM): the array never shuts once it fires (4/4
unended, 70 % of channel-time open). Fitted to 37 °C (Ka 16.9, Ki 358 µM):
no sparks at all. Ki scan: 19/18/17 ms at 3/5.5/10 µM, 25/33/89 ms at
19/35/64 µM, 1.4 s at 118 µM, unended from 217 µM. Rate scan at Ki 10:
162 ms at 1/30×, down to 13–15 ms at 3–10×, where sparks nearly stop
starting (none at 30×). At the fitted Ka and Ki, no rate ends a spark.

**What it means.** The roadmap's suspect was wrong in direction. Stern's
over-sensitive inactivation did not lengthen sparks. It is the only
thing in the scheme that ends them. With the measured sensitivity, the
tens of µM a channel sees in the cleft cannot inactivate enough of the
array. The measured bell has no Mg²⁺, which is the obvious missing
terminator (with luminal depletion and the V channels). All three are
recorded as emergent rather than assumed.

**Not changed.** `make sync-check` was clean, no checks were touched, and
the verdicts are unchanged (44 confirmed, 2 discrepancies).

**Next:** source a Mg²⁺-inclusive RyR1 bell or scheme with readable
constants before modelling Mg²⁺.


## 2026-09-24 — Round 6.5: Mg²⁺ ends a triggered cleft spark

**What.** `SternParams` gains free Mg²⁺ (`mg`, default 0) with two
sourced actions. At the activation site it competes in rapid
equilibrium: the on rate is divided by (1 + Mg/K_Mg,A)², so Ka becomes
Ka (1 + Mg/K_Mg,A). At the inactivation site it binds like Ca²⁺
(`mg_i` = 1). `with_mg` builds a scheme under Mg²⁺. `fit_to_bell` now
fits at Mg²⁺ 0, because the bell was measured without it.
`simulate_sparks_cleft(trigger=True)` opens every channel that is not
inactivated at t = 0. `physics/spark_mg.py` times triggered sparks,
dissects the two sites and scans Mg²⁺, and the CLI gains `spark-mg`.
There are 7 new parameters and 3 new references.
`tests/test_spark_mg.py` has 10 tests (309 total). A planted fault that
disables the activation-site competition fails 3 of them (checked, then
restored).

**Sources, and why these.** The roadmap asked for a Mg²⁺-inclusive RyR1
bell with readable constants first. The best match is Meissner et al.
1997 (JBC), which used the same [³H]ryanodine assay as Murayama. It sits
behind a Cloudflare challenge, and Chrome was not connected. The open
sources give the two Mg²⁺ actions separately, which a two-gate scheme
needs anyway:
- Laver 2004 Table I: K_Mg,A 54 µM, beside a Ca²⁺ affinity of 0.51 µM,
  measured with ATP.
- Laver 1997: Ca²⁺ and Mg²⁺ inhibit identically at the I1 site.
- Laver 2018: 1 mM free Mg²⁺ in the fibre.

Nothing was fitted. Because K_Mg,A's ATP condition does not match
Murayama's, there is a second reading: Laver's selectivity carried onto
the fitted Ka, 521 µM. Both are reported.

**Why a trigger.** Under fibre Mg²⁺ no spark starts by itself, so the
spontaneous ruler can only say "none". In the fibre, the V channels start
sparks, and they are not simulated. Opening the available channels at
t = 0 is the simplest stand-in that asks the question at hand: once lit,
does the array shut?

**Measured** (fitted to 25 °C, 1 mM Mg²⁺, 20 trials):
- No Mg²⁺: never shuts (0/20 in 2 s).
- Activation site only: 20/20 shut, 17.5 ms (54 µM) or 49.5 ms (521 µM),
  with 0 channels inactivated.
- Inactivation site only: 6 open, 25 ms.
- Both sites: 6 open, 4–6 ms.
- Scan (both sites, 54 µM): the array never shuts up to 25 µM, then
  131 / 15 / 10.5 / 4 ms at 63 / 158 / 398 / 1,000 µM.

**What it means.** Mg²⁺ is the missing terminator, but it does not work
through inactivation. Competition at the activation site drops the
array's positive feedback below one, so a spark decays stochastically
(induction decay). The I1-site action mostly removes channels before
the spark. Laver 2018 proposes that the voltage sensor lifts that block
during E–C coupling. The activation-site row (17.5–49.5 ms) would then
be the physiological one, still 3–8× the 6.3 ms frog release. The
both-sites row matches the duration only by leaving 24 of 30 channels
shut.

**Process note.** The Unpaywall lookup for Meissner 1997 was made with
the user's email address as its required `email` parameter. That should
not have been sent. Future lookups use Crossref / Europe PMC / Semantic
Scholar, which need no address.

**Not changed.** `make sync-check` was clean, no checks were touched, and
the verdicts are unchanged (44 confirmed, 2 discrepancies). No UI files
changed, so no screenshots were needed.

**Next:** pin K_Mg,A in Murayama's condition (Meissner 1997; the user
may need to supply the PDF), then the V-channel trigger.


## 2026-09-24 — Round 6.5 closed; Mg²⁺ in the GUI

**Closed 6.5.** The Round 6.5 work had been written up but not committed.
Tests (309), lint and sizes passed as left, so it was committed and pushed
as it stood.

**Why not the logged next step.** Pinning K_Mg,A needs Meissner et al.
1997 (JBC 272:1628). Its abstract confirms Mg²⁺ competes at the Ca²⁺
activation site in the [³H]ryanodine assay, and the article is CC-BY. But
the Elsevier API returns metadata only, and jbc.org, ScienceDirect and the
Semantic Scholar PDF link all sit behind Cloudflare. Chrome was not
connected. The V-channel trigger needs the 10-state allosteric model of
Ríos, Karhanek, Ma & González 1993 (J Gen Physiol 102:449), whose
constants Stern 1997 cites but does not print (Table I is C channels
only). That PMC record is a page scan whose PDF sits behind the
download challenge, and so are its 1993–94 companions. Both items
therefore need the user to supply PDFs, so this session took the next
self-contained item: Mg²⁺ in the GUI.

**What.**
- `puff_compare.params_for(mg=, k_mg_a=)` puts any spark receptor under
  Mg²⁺ and refuses an IP3R receptor.
- `spark_mg.READINGS` / `k_mg_a_reading` resolve the two K_Mg,A readings
  in one place, now used by both the CLI and the GUI.
- Puffs: a free-Mg²⁺ spin box, the reading, and "Triggered sparks vs Mg²⁺"
  (cleft receptors; worker; median duration, opened and inactivated
  against Mg²⁺). Values are read on the main thread and passed to the
  worker.
- Gating: the fitted bell under 1 mM Mg²⁺ is drawn relative to the
  Mg²⁺-free peak, and the text gives both readings. The peak falls to 17 %
  with half-activation at 77 µM (54 µM reading), or 21 % and 13 µM (521 µM).
- 2 new tests (311 total).

**Found by the smoke test.** A cluster with no events (the fitted cleft
under 1 mM Mg²⁺: nothing starts by itself) crashed the Puffs panel. The
event-size axis was set to log with no positive data. The panel now says
"no events". That was a real bug for any silent cluster, not only under
Mg²⁺.

**Structure.** `scripts/screenshot_app.py` was at 485 lines, so the spark
steps moved to `scripts/screenshot_sparks.py` (`spark_step`). That file
adds two steps: the fitted cleft under Mg²⁺ must be silent, and the
triggered scan must never shut at Mg²⁺ 0 and always shut at 1 mM.
New screenshot: `gui_sparks_mg.png`. The overall hang timeout rose from
600 to 900 s. Unloaded, the full run takes ~130 s, but this afternoon the
GUI process got ~35 % of a core (macOS throttling a background window;
`caffeinate` did not help) and reached the new scan at ~560 s. Leftover
instances from killed runs had also been competing: check `pgrep -f
screenshot_app` before rerunning. The passing run took 620 s.

**Not changed.** `make sync-check` was clean, no checks were touched, and
the verdicts are unchanged.

**Next:** the user supplies Meissner 1997 (K_Mg,A without ATP) and Ríos
et al. 1993 (the V-channel model). Without them, the next open item is
Round 2's gate radius along the morph.


## 2026-09-24 — Round 2: the gate along the morph

**Why this item.** Meissner 1997 and Ríos et al. 1993 have still not been
supplied, so this session took the logged fallback: Round 2's unmeasured
gate radius along the morph.

**The problem.** `displaced_coords` carries each side chain rigidly with its
Cα. So the last frame of 8TKG → 8TKF was 8TKF's backbone lined by 8TKG's
rotamers, and any pore profile of a morph frame was wrong at exactly the
gate.

**What.**
- `structure/morph_pore.py`:
  - `atom_path` matches every heavy atom of the basis residues by chain,
    residue and name, and carries the end through the transition's own
    superposition (now kept as `meta["end_transform"]`). Each atom's offset
    from its Cα is interpolated linearly while the Cα follows the morph, so
    both ends are the deposits atom for atom.
  - `gate_path` re-finds the axis and the pore-domain span on every frame
    and applies the deposit rule. It measures both ways: interpolated and
    rigid.
  - `offset_error` is the chord a swinging side chain takes. It is up to
    3.9 Å over all atoms, but ≤ 0.23 Å on the gate's lining.
- `pore.py`: `radial_profile` and `constriction_indices` were factored out
  of `pore_profile` and `find_constrictions` (S0's windows are now named
  constants), so the morph uses the same rule rather than a copy.
- Viewer: `coords_at` now draws the interpolated atoms. The Transition tab
  plots gate, rigid gate and filter along the path beside the overlap plot
  (`PlotCanvas.span_row`), and the frame label shows the gate.
- CLI: `transition --gate`.
- Tests: `test_morph_pore` (7). The smoke test also asserts the gate
  opens from < 3 to > 5.5 Å.

**Measured.**

| | start gate | end gate | rigid at end | half-way |
|---|---|---|---|---|
| ITPR3 8TKG → 8TKF | 2.73 Å (F2513) | 5.85 Å (N2510) | 5.16 Å | 0.42 |
| RyR1 9R8O → 9HEO | 3.32 Å (I4937) | 5.05 Å (Q4933) | 5.24 Å | 0.44 |

- Endpoint frames equal the deposits' own measurements (protein only).
- There is no overshoot.
- In ITPR3 the constriction moves 7 Å down, from I2517 to N2510, at
  t ≈ 0.7.
- The shortcut gets RyR1's end gate wrong in the other direction (too
  wide), so it has no consistent bias.
- The half-way point is a property of a linear path, not of gating order.
  `docs/SCIENCE.md` says so.

**Not changed.** Sync was clean, and no check or verdict moved.

**Next:** still the PDFs (Meissner 1997, Ríos 1993). Without them, the next
item is Round 2's stride-3 local modes.

## 2026-09-24 (2) — Round 2: where the local modes come from

**Why this item.** The PDFs (Meissner et al. 1997, JBC 272:1628,
doi:10.1074/jbc.272.3.1628; Ríos et al. 1993) have still not been supplied.
The user asked which ones are wanted and was told. The logged fallback was
Round 2's stride-3 local modes.

**What.**
- `physics/network_checks.py`:
  - `local_modes` puts every low-κ mode on the sequence: the residue
    carrying most of it (summed over subunits), and the unresolved stretch
    between that residue's sampled neighbours. `describe_local` groups the
    modes by residue.
  - `cutoff_scan` re-solves the network over a registered grid
    (`anm.scan_cutoff_low/high/step`: 12–21 Å in 1.5 Å steps).
  - `rmsip` and `stride_agreement` score a strided network against
    stride 1.
- `transition_overlap(cutoff=)`. The report now names where the local
  modes sit and gives the collective A modes together.
- CLI: `transition --cutoff-scan --stride-check`; `modes` prints κ and
  where the local modes sit. Modes tab: the status line and a selected
  local mode name the residue.
- The `anm.stride` source note claimed the low modes were "insensitive" to
  the stride. That was measured false and has been rewritten with the RMSIP
  values.
- Tests: `test_network_checks` (8). A planted flap after a gap is located,
  and the same slab without it has no local modes. The real-data tests pin
  8TKG's five and the cutoff behaviour.

**Measured.**
- 8TKG stride 3: #11–15 (B E E A B) are all residue 86, after the
  unresolved 77–85. Residue 86 has 5 springs, against a median of 17. At
  stride 4 the same flap takes modes 2–5. Strides 1 and 2 have no local
  modes, and neither does 6DQN.
- RyR1 9R8O has its own local modes: residue 896 at strides 1 and 3,
  residue 1988 at stride 2.
- Sequence-neighbour springs at 1, 10 and 100 × γ change nothing. The flap
  swings as a body, and a chain spring resists only stretching. So bridging
  it is answered: no.
- A constant-coordination cutoff (15 Å × stride^⅓) removes the flap. But it
  changes the model: the headline number moved.
- RMSIP of 20 modes against stride 1: 0.970 at stride 2, 0.892 at stride 3,
  0.769 at stride 4.

**The finding.** At stride 1, where there is no sampling artefact, the
lowest-A-mode overlap is 0.48 / 0.45 / 0.43 / 0.38 / 0.33 / 0.24 at cutoffs
of 12 / 13.5 / 15 / 16.5 / 18 / 21.6 Å. Over the same range the A modes
among the first 10 together give 0.51 / 0.66 / 0.67 / 0.67 / 0.67 / 0.68.
- The cutoff only redistributes the move among near-degenerate A modes.
- The logged "0.415, one A mode carries two-thirds of 20 modes" held only at
  15 Å. The A subspace is the result, and README and SCIENCE.md now say so.
- Round 2's robustness test had varied the stride but never the cutoff.
- For RyR1 the lowest A mode (0.12) is flat in the cutoff.

**Not changed.** Defaults (cutoff 15 Å, stride 3, κ 0.2) are unchanged.
Sync was clean. No check or verdict moved. `make test` gave 326 passed
before the new tests (334 with them), lint and sizes were clean, and
`make screenshots` passed.

**Next:** the PDFs have now been supplied (`pdfs/meissner_1997.pdf`,
`pdfs/rios_1993.pdf`, git-ignored; both have text), so do Round 6.5's
emergent items: pin K_Mg,A from Meissner 1997, then the V-channel
trigger from Ríos 1993. The user's Mendeley library
(`~/Documents/Mendeley Desktop`) has about 2,100 more PDFs.
