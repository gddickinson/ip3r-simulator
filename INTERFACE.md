# INTERFACE.md — navigation map

Read this before opening any source file. Every `.py` is under 500 lines.

## Dependency direction

```
io ──▶ core ──▶ structure ──▶ physics ──▶ analysis
                                  │
                   render ◀───────┴───────▶ ui
```

`physics` and below never import `render` or `ui`, so the science runs
headless (CLI, tests, notebooks).

## Repository layout

| Path | Contents |
|---|---|
| `ip3r/` | the application package |
| `scripts/` | maintenance: parameter build, ip3r_genes import, GUI smoke test |
| `tests/` | pytest suite; real-data tests skip when data is absent |
| `docs/SCIENCE.md` | the models, their equations and sources; what each check establishes |
| `docs/img/` | screenshots written by `scripts/screenshot_app.py` |
| `ROADMAP.md` | what is not done, in rounds |
| `SESSION_LOG.md` | what was done and why, per session |
| `Makefile` | every task (`make help`) |
| `run_app.command` | launcher: activates `ip3r_sim` and runs `python -m ip3r` (double-clickable on macOS; arguments pass through) |
| `ref/`, `data/` | downloads and derived output — git-ignored, regenerable |

## `ip3r/` — top level

| File | Purpose |
|---|---|
| `config.py` | paths (`RESOURCE_DIR`, `REF_DIR`, `GENES_DIR` = `../ip3r_genes` or `$IP3R_GENES_DIR`, `genes_results()` resolved at call time), `PARALOG_ACC`, `DEFAULT_STRUCTURE` (6DQN), `RenderSettings` |
| `parameters.py` | the parameter registry (`PARAMETERS.value(key)`, overrides tracked, `IP3R_PARAMETERS` override file). Ported from PIEZO1. |
| `cli.py` | `python -m ip3r <fetch|info|checks|states|modes|transition|gating|oscillate|puffs|params>`; no argument launches the GUI |
| `__main__.py` | entry point |

## `ip3r/io/`

| File | Key names |
|---|---|
| `cif_reader.py` | `read_structure_file()` — fast mmCIF/PDB → numpy arrays (ported) |
| `registry.py` | `StructureEntry`, `load_registry()`, `get_entry()`, `local_path()` — the 9 curated depositions from `resources/structures.json` |
| `fetch.py` | `fetch_structure()`, `fetch_all()`, `is_valid()` — RCSB `.cif.gz` into `ref/structures`; `IP3R_STRUCTURE_MIRROR` copies from a local mirror (e.g. the ip3r_genes data root) |
| `loader.py` | `load(pdb_id)` memoised; `ALLOW_FETCH` switch (off by default; the GUI and `--fetch` turn it on); `StructureUnavailable` |

## `ip3r/core/`

| File | Key names |
|---|---|
| `structure.py` | `Structure` — structure-of-arrays container, masks, residue index (ported) |
| `annotations.py` | per paralog: `reference_sequence`, `elements` (Pfam + 6DQN structural elements), `residue_elements` (one element per residue, specific wins), `element_of/element_array`, `functional_sites` (10 IP3 contacts, filter/gate lining), `residue_constraint/constraint_at` (S17 JSD, 4 layers), `variants`; colour/label tables |
| `modules.py` | Paper 6's modules rebuilt from sites + domains: `module(paralog, definition)` (`contact_span`, `channel_minus_luminal`, `channel_all`) → `Module`, `modules()` (the primary pair, disjoint), `ModuleRefusal`, `MODULE_COLORS`, `MODULES_KEY` (the GUI's site toggle) |
| `pairwise.py` | `align` (Gotoh affine-gap global, BLOSUM62, end gaps free; gap costs `align.*`), `transfer_map`, `paralog_transfer(src, dst)` (memoised) — carries residue numbers between paralogs independently of S17's MAFFT |
| `genes_data.py` | live read-only access to `ip3r_genes/results`: `read_tsv`, `read_json`, `read_text`, `available`; raises `GenesDataMissing` (→ check `not_run`) |

## `ip3r/structure/` — measurement

| File | Key names |
|---|---|
| `symmetry.py` | `kabsch`, `rotation_matrix`, `rotation_axis_angle`, `axis_by_superposition` (our method), `axis_by_centroids` (S0's), `tetramer_frame` → `Frame` (z on the axis, cytosol +z, chains right-handed), `c4_residual` |
| `pore.py` | `pore_profile` (`r_min` = S0's quantity; `r_free` = less vdW), `tm_span`, `find_constrictions` (filter = luminal min, gate = cytosolic min), `lining_residues`; `include_hetero` reproduces S0 |
| `ligand.py` | `ligand_sites` (IP3 copies, subunit by proximity), `contacts` (≤ cutoff, own vs other subunit), `residue_distances` (S22 shells; `heavy_only=False` = S22's all-atom rule) |
| `shells.py` | S22's ligand shells: `shell_edges`/`shell_of` (registered edges, half-open), `chain_distances` (own-subunit IP3, all atoms), `deposit_distances`, `consensus_shells` → `ShellResidue` (median over deposits), `atom_ligand_distance` (per atom, for painting) |
| `numbering.py` | `check_numbering`, `best_numbering` — S24's rule, stubbed (backbone+CB) residues excluded, mismatch segments reported |
| `channel.py` | `measure_channel(st)` → `ChannelSummary` (axis two ways, residual, numbering, span, profile, constrictions, IP3 contacts) — shared by GUI, CLI and checks |
| `transition.py` | `prepare_transition(start, end, fit)` → `Transition` (residue-matched basis: unstubbed, sequence-matching, all 8 chains; cyclic subunit correspondence; end superposed onto the start *as deposited*; `residue_distance`, `element_means`), `atom_site_index` (own residue, else nearest site in space), `displaced_coords`, `atom_displacement` (NaN off basis), `TransitionUnavailable` |
| `morph.py` | `morph(start, end, method)` → `MorphTrajectory` (`restrained` / `linear`, `bond_error`, `nearest`), `peptide_pairs`, `NOTE` (the "interpolation, not trajectory" sentence) |
| `states.py` | `state_panel(paralog)` → `StateRow`s: every human deposit measured the same way (the gating transition at the pore) |

## `ip3r/physics/` — simulation

| File | Key names |
|---|---|
| `anm.py` | `build_hessian` (inverse-square springs), `ANM.calc_modes` (drops 6 × components), `ANM.label_symmetry` (C4 irreps A/B/E), `apply_generator`, `ModeSet` (`collectivity` κ, `is_collective`, `first(irrep)` skips local artefacts), `tetramer_sites`, `atom_displacements` |
| `transition_modes.py` | `transition_overlap(tr, reference)` → `TransitionOverlap` (overlap, cumulative, symmetry-matched null, irrep fractions, `report()`), `remove_rigid_body`, `irrep_fractions`, `null_cumulative` |
| `gating.py` | De Young–Keizer / Li–Rinzel: `m_inf`, `n_inf`, `q2`, `h_inf`, `tau_h`, `open_probability`, `bell_peak`, `hill_fit_left_flank`, `GatingParams` |
| `calcium.py` | closed-cell Li–Rinzel: `simulate` → `Trace`, `fluxes`, `oscillation_metrics` (sustained only), `oscillation_window` (0.36–0.63 µM measured), `steady_state`, `CellParams` |
| `puffs.py` | stochastic DYK cluster: `simulate_cluster` → `PuffTrace`, `detect_events`, `fano`, `coupling_effect`, `PuffParams` |

## `ip3r/analysis/` — the findings checks

| File | Key names |
|---|---|
| `checks.py` | framework: `register` decorator, `Check`, `Outcome` (`confirmed/discrepancy/not_run/error`), `run_check(s)`, `all_checks`, `PAPERS`; refuses to confirm against a modified registry |
| `checks_structure.py` | `S0.*` (C4, filter, gate, pore profile, IP3 contacts — recomputed from coordinates) and `P6.shell_agreement`, `P6.contacts_heavy_atom` |
| `checks_constraint.py` | `P5.*` (element means, rankings, gate identity, variant AUCs under S17's position rules, deep-ranks-third, VUS count, ω) and `P6.contacts_vs_core` |
| `checks_evolution.py` | `P1.absences`, `P2.sister_pair` (own Newick reader), `P2.au_test`, `P2.teleost_itpr1`, `P3.*`, `P4.unreachable`, `LEDGER.claims` |
| `checks_modules.py` | `P6.module_map`, `P6.module_contrast`, `P6.loop_reverses` |
| `checks_shells.py` | `P6.shell_distances` (recomputed), `P6.shell_constraint`, `P6.shell_trend`, `P6.no_contact_step` |
| `shell_constraint.py` | `measured_shells` (the 6-deposit pocket), `pocket(gene)` → `Pocket` (carried by own alignment, joined to deep JSD), `shell_rows`, `trend`, `contact_step`, `clear_caches` |
| `module_contrast.py` | `read_alignment`, `reference_columns` (own residue → column map, refused on a sequence mismatch), `tip_identities`, `paired_contrast` → `PairedContrast`, `clear_caches` |
| `stats.py` | `auc` (Mann-Whitney, ties ½), `rank_average`, `mean_by_group`, `sign_test` (exact), `signed_rank_test` (normal approx., tie-corrected), `mann_whitney_greater` (one-sided, tie + continuity corrected), `spearman` (t-distribution p), `fisher_exact` (two-sided), `logistic_fit` (IRLS, Wald p), `wilson` |
| `newick.py` | `parse`, `leaves`, `mrca`, `smallest_clade_containing` |
| `tree.py` | Paper 2's clade questions on `rooted.nwk`: `load_tree`, `group_of` (census prefix), `is_cyclostome`/`genus_of`, `support`, `is_supported` (registered bars), `bipartitions` (root's twin edge counted once), `paralog_clades` → `Clade` (MRCA of named tips; named/unnamed/foreign), `outgroup_clade`, `cyclostome_clades` → `CyclostomeClade` |
| `tree_figure.py` | `layout` (ladderized phylogram), `draw_tree(ax, root, labels, supports)` — clades boxed, cyclostomes marked, supported nodes dotted; `CLADE_COLORS` |
| `checks_tree.py` | `P2.paralog_clades`, `P2.cyclostome_lineages`, `P2.support_bar` (rederived from the tree; the Newick travels in the outcome for the exhibit) |
| `genome_grid.py` | Papers 3/4 per cell: `load_grid(layers)` → `GenomeGrid` (309 genomes × `CELLS` ITPR1–3 + RYR; layers `search`/`miss`/`state`/`recovery`; `counts`, `row`, `subset`), `Genome` (N50, level, source, `above_bar`), `recovery_channel(row)` (rebuilt from counts), `above_bar` (registered `genomes.contiguity_bar_bp`), `order`/`ORDERS` |
| `grid_figure.py` | `draw_grid(ax, grid, layer)` (fixed category colours, N50 strip on a fixed log scale, bar line, class blocks), `LAYER_STYLE`, `LAYER_TITLES`; exhibits `draw_misses`, `draw_logistic`, `draw_recovery` |
| `checks_genomes.py` | `P3.miss_by_contiguity`, `P3.contiguity_tests` (every row of contiguity_tests.tsv), `P4.recovery_channels` (rederived) |
| `exhibits.py` | `draw(ax, check_id, outcome)` — figures from a check's own numbers |

## `ip3r/render/` (moderngl, OpenGL 4.1)

`camera.py`, `primitives.py`, `scene.py`, `spline.py`, `geometry_builders.py`,
`shaders/` — ported unchanged from PIEZO1 (impostor spheres/cylinders,
cartoon sweeps, trackball camera). `colormaps.py` — chain, element, fixed
conservation ramp (0.50–0.95 JSD; grey = not scored), fixed displacement ramp
(0–25 Å), `SHELL_COLORS`/`shell_colors` (S22's four shells, grey beyond). `representations.py` —
`MolecularView` (styles × `ColorBy`, highlight (uniform or per-atom `highlight_rgb`), chain filter, `update_coords`
for animation; `ColorBy.DISPLACEMENT` from a built transition; `ColorBy.LIGAND_SHELL` from `structure.shells`).

## `ip3r/ui/` (PyQt6)

| File | Purpose |
|---|---|
| `app.py` | `main()` — surface format, theme, window, initial load |
| `main_window.py` | layout and wiring; menus; loads on workers; `CHECK_SITES` maps a check to what "Show on structure" highlights, `CHECK_COLOURS` to a colouring (the shell checks), `CHECK_TREE` to the Tree tab, `CHECK_GENOMES` to the Genomes tab and its layer |
| `scene_controller.py` | what the viewport draws: `MolecularView`, pore spheres, site/variant highlights, side/top views, mode animation |
| `gl_widget.py` | `ViewportWidget` (ported; viewport sized from the bound FBO every frame) |
| `structure_panel.py` | deposition list, style, colour, layer, subunits, measured sites, legend |
| `channel_panel.py` | `ChannelSummary` text, pore profile vs S0's, ITPR3 state comparison |
| `modes_panel.py` | ANM table with irreps and κ, animation controls |
| `transition_panel.py` | Transition tab: end state, fit, method, 8TKG→8TKF preset, frame slider/play, displacement colouring, element and overlap plots |
| `transition_controller.py` | `build_transition` (worker), `TransitionController` (install, `coords_at`, `show_frame`, `play`, `reset` — path built from the displayed structure) |
| `dynamics_panel.py` | Gating (bell), Oscillations (+ window scan), Puffs (coupled vs uncoupled) |
| `findings_panel.py` | checks by paper, run on a worker, claim/method/verdict, exhibit, show on structure (`showable`: residue-keyed checks drawn on the displayed structure) |
| `tree_panel.py` | Tree tab: `draw_tree` on a toolbar canvas, tip labels / support toggles, "Vertebrates" zoom, click names a tip; loaded on a worker when first shown |
| `genomes_panel.py` | Genomes tab: `draw_grid` with layer / sort / class / above-bar controls, click names a genome; loaded on a worker when first shown; `show_layer` (from a check's "Show") |
| `variants_panel.py` | S17 variants per paralog/class; highlight only in matching numbering |
| `params_dialog.py`, `plot_canvas.py`, `workers.py`, `theme.py` | helpers |

## `ip3r/resources/` (committed)

`parameters.json`, `references.json` (built by `scripts/build_parameters.py`);
`sequences.json`, `domains.json`, `sites.json`, `constraint.json`,
`variants.json`, `structures.json` (built by `scripts/sync_genes.py` from
ip3r_genes, each with the source paths, SHA-256 and ip3r_genes commit).

## `scripts/`

| File | Purpose |
|---|---|
| `parameter_table.py`, `reference_table.py`, `build_parameters.py` | the registry and its provenance gate |
| `sync_genes.py` | import resources from ip3r_genes; `--check` reports drift |
| `screenshot_app.py` | scripted GUI smoke test + README screenshots |
| `create_env.sh` | the `ip3r_sim` conda env |

## `tests/`

Transition (`test_transition` synthetic calibrations; `test_transition_real`
— the drawn end must be 8TKF as a shape, with a case that must fail), physics (`test_anm`, `test_gating`, `test_calcium`, `test_puffs`), geometry
(`test_symmetry`, `test_pore`, `test_structures_real`), statistics
(`test_stats` — Fisher vs scipy and the tea-tasting value, logistic slope = log OR for a binary predictor, Wilson vs Newcombe; `test_newick`), grid (`test_genome_grid` — the channel rule, the bar, ordering, the real grid's counts), alignment (`test_pairwise` — score equals a cell-by-cell reference DP), shells (`test_shells`), tree (`test_tree` — toy trees with known clades, root twin edge, a misplaced tip is foreign), modules (`test_modules` — spans hold their sites, column map lands on the residue, a tampered reference is refused), provenance (`test_parameters`,
`test_resources`), rules (`test_sizes`), CLI, and
**`test_checks_calibration.py`** — every check flipped by a planted input
(`PLANTS`), sources proven complete, `not_run` without data, refusal under
modified parameters.
