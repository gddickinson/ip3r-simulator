# ROADMAP — what is not done

Work proceeds in rounds of one focused task each. Every round: implement,
test (`make test`, `make screenshots` if the UI changed), update the docs,
commit, push. `[ ]` planned, `[x]` done (with what it measured). The
completed Round 1 is recorded in `SESSION_LOG.md`.

**Next: Round 3.**

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

- [ ] Paper 6 module contrast (`module_contrast.tsv`): re-derive the paired
  per-orthologue core-vs-pore identities from the deep alignments
  (`aln_ITPR*.fasta`) and the module map; highlight both modules.
- [ ] Ligand shells: colour residues by all-atom distance to IP3 (S22's
  shells) and plot conservation against distance (the "no step at 4.5 Å"
  result).
- [ ] Paper 2: a tree viewer for `rooted.nwk` with the paralog clades boxed
  and the cyclostome tips marked.
- [ ] Paper 3/4: a genome × paralog grid of the character matrix and the
  recovery channel, sortable by contiguity.
- [ ] Paper 1: presence/absence across eukaryotic clades from the S20/S23
  tables.

## Round 4 — better physics

- [ ] Mak et al. 1998 Hill-type gating model, in which IP3 tunes Ca²⁺
  inhibition alone (the DYK model shifts both flanks; see SCIENCE.md) —
  every constant verified against the paper before it is registered.
- [ ] A puff model with low resting activity (e.g. a Siekmann-type park/drive
  scheme) so blips and puffs separate cleanly; the DYK cluster gives only
  modest clustering (Fano ~1.4).
- [ ] Unitary current from the pore profile (1-D drift–diffusion, as PIEZO1's
  `physics/permeation.py`), compared across the state panel.

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
