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
