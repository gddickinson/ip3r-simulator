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
| `docs/SCIENCE_RYR.md` | the same for the ryanodine receptor: structures, conductance, mutants, gating, sparks (mean-field and in the cleft) |
| `docs/SCIENCE_EC.md` | the couplon under voltage clamp: V channels, Stern's calibration, SR depletion |
| `docs/img/` | screenshots written by `scripts/screenshot_app.py` |
| `ROADMAP.md` | what is not done, in rounds (Round 6, the ryanodine receptors, in `ROADMAP_RYR.md`) |
| `SESSION_LOG.md` | what was done and why, per session |
| `Makefile` | every task (`make help`) |
| `run_app.command` | launcher: activates `ip3r_sim` and runs `python -m ip3r` (double-clickable on macOS; arguments pass through) |
| `ref/`, `data/` | downloads and derived output — git-ignored, regenerable |

## `ip3r/` — top level

| File | Purpose |
|---|---|
| `config.py` | paths (`RESOURCE_DIR`, `REF_DIR`, `GENES_DIR` = `../ip3r_genes` or `$IP3R_GENES_DIR`, `genes_results()` resolved at call time), `PARALOG_ACC` / `PARALOGS` (the publication's three), `RYR_ACC` (rabbit RyR1 P11716), `NUMBERINGS` (both: what a deposit can be numbered in), `DEFAULT_STRUCTURE` (6DQN), `RenderSettings` |
| `parameters.py` | the parameter registry (`PARAMETERS.value(key)`, overrides tracked, `IP3R_PARAMETERS` override file; `subscribe` for change listeners — caches of parameter-dependent results and the GUI banner; `matches` (editor filter), `write_overrides`/`read_overrides`, `replace` (whole set, one notification); `references()` for tooltips). Ported from PIEZO1. |
| `cli.py` | `python -m ip3r <fetch|info|graft|checks|states|unitary|modes|transition [--cutoff-scan --stride-check]|gating|oscillate|puffs|params>`; no argument launches the GUI |
| `cli_ryr.py` | the RyR1 commands, `register(sub)`: `mutants`, `ryr-gating`, `sparks [--scan] [--cleft [--fit]]`, `spark-termination [--scan fit\|ki\|rate] [--fitted]`, `spark-mg [--scan] [--reading R] [--spontaneous] [--two-site]`, `ec [--scan] [--reading R] [--only S] [--trials N] [--two-site] [--depletion [--pool-scan] [--large]]` |
| `__main__.py` | entry point |

## `ip3r/io/`

| File | Key names |
|---|---|
| `cif_reader.py` | `read_structure_file()` — fast mmCIF/PDB → numpy arrays (ported) |
| `registry.py` | `StructureEntry` (`family`: IP3R / RyR), `load_registry()`, `get_entry()`, `local_path()` — the 9 IP3R depositions from `resources/structures.json` plus the 6 RyR1 ones from `resources/ryr1.json` |
| `fetch.py` | `fetch_structure()`, `fetch_all()`, `is_valid()` — RCSB `.cif.gz` into `ref/structures`; `IP3R_STRUCTURE_MIRROR` copies from a local mirror (e.g. the ip3r_genes data root) |
| `loader.py` | `load(pdb_id)` memoised; `ALLOW_FETCH` switch (off by default; the GUI and `--fetch` turn it on); `StructureUnavailable` |
| `predictions.py` | AlphaFold DB models into `ref/alphafold` (`fetch_prediction` discovers entry + version from the API; `ACCESSIONS`, `UNAVAILABLE` with why), `local_predictions`, `load_prediction` (memoised; pLDDT in `b_factor`) |
| `session.py` | `Session` (the view: deposit, style, colour, layer, subunits, sites, pore, completeness, camera, tab, transition spec, and the parameter overrides it was saved under — never coordinates or results), `from_dict` refuses a wrong type or newer format by name, `save_session`/`load_session`, `parameter_differences` |

## `ip3r/core/`

| File | Key names |
|---|---|
| `structure.py` | `Structure` — structure-of-arrays container, masks, residue index (ported) |
| `annotations.py` | per paralog (or `RYR1`, `is_ryr`: sequence and domains from `ryr1.json`; no constraint/sites/variants, so grey): `reference_sequence`, `elements` (Pfam + 6DQN structural elements), `residue_elements` (one element per residue, specific wins), `element_of/element_array`, `functional_sites` (10 IP3 contacts, filter/gate lining), `residue_constraint/constraint_at` (S17 JSD, 4 layers), `variants`; colour/label tables |
| `modules.py` | Paper 6's modules rebuilt from sites + domains: `module(paralog, definition)` (`contact_span`, `channel_minus_luminal`, `channel_all`) → `Module`, `modules()` (the primary pair, disjoint), `ModuleRefusal`, `MODULE_COLORS`, `MODULES_KEY` (the GUI's site toggle) |
| `pairwise.py` | `align` (Gotoh affine-gap global, BLOSUM62, end gaps free; gap costs `align.*`), `transfer_map`, `paralog_transfer(src, dst)` (memoised) — carries residue numbers between paralogs independently of S17's MAFFT |
| `genes_data.py` | live read-only access to `ip3r_genes/results`: `read_tsv`, `read_json`, `read_text`, `available`; raises `GenesDataMissing` (→ check `not_run`) |

## `ip3r/structure/` — measurement

| File | Key names |
|---|---|
| `symmetry.py` | `kabsch`, `rotation_matrix`, `rotation_axis_angle`, `axis_by_superposition` (our method), `axis_by_centroids` (S0's), `tetramer_frame` → `Frame` (z on the axis, cytosol +z, chains right-handed), `c4_residual` |
| `pore.py` | `pore_profile` (`r_min` = S0's quantity; `r_free` = less vdW; `radial_profile` on atoms already in frame), `tm_span`, `find_constrictions` (filter = luminal min, gate = cytosolic min; the rule itself is `constriction_indices`, S0's windows `PROFILE_MARGIN`/`GATE_MARGIN`), `lining_residues`; `include_hetero` reproduces S0 |
| `ligand.py` | `ligand_sites` (IP3 copies, subunit by proximity), `contacts` (≤ cutoff, own vs other subunit), `residue_distances` (S22 shells; `heavy_only=False` = S22's all-atom rule) |
| `shells.py` | S22's ligand shells: `shell_edges`/`shell_of` (registered edges, half-open), `chain_distances` (own-subunit IP3, all atoms), `deposit_distances`, `consensus_shells` → `ShellResidue` (median over deposits), `atom_ligand_distance` (per atom, for painting) |
| `numbering.py` | `check_numbering`, `best_numbering` — S24's rule, stubbed (backbone+CB) residues excluded, mismatch segments reported |
| `channel.py` | `measure_channel(st)` → `ChannelSummary` (axis two ways, residual, numbering, span, profile, constrictions, IP3 contacts) — shared by GUI, CLI and checks |
| `transition.py` | `prepare_transition(start, end, fit)` → `Transition` (`meta["end_transform"]` = the end's superposition;residue-matched basis: unstubbed, sequence-matching, all 8 chains; cyclic subunit correspondence; end superposed onto the start *as deposited*; `residue_distance`, `element_means`), `atom_site_index` (own residue, else nearest site in space), `displaced_coords`, `atom_displacement` (NaN off basis), `TransitionUnavailable` |
| `morph_pore.py` | the gate along a morph: `atom_path(start, end, tr)` → `AtomPath` (heavy atoms matched by name; Cα offset interpolated, `coords`, `rigid` = the old shortcut, `offset_error` chord, `index` into the start deposit), `gate_path(path, tr, mt)` → `GatePath` (gate/filter per frame, `rigid_gate`, lining, `half_open`, `overshoot`); axis re-found per frame |
| `morph.py` | `morph(start, end, method)` → `MorphTrajectory` (`restrained` / `linear`, `bond_error`, `nearest`), `peptide_pairs`, `NOTE` (the "interpolation, not trajectory" sentence) |
| `graft.py` | AlphaFold fills: `prediction_for(st)` (the model in the deposit's numbering, else `GraftRefusal`), `unresolved` → `Stretch` (gap / n_term / c_term), `fill_stretches`, `fill_structure(st, mode)` → `FilledModel` (predicted atoms only; `Fill` per stretch: anchor RMSD, seam distances, pLDDT, clashes; `Skip` with reason; `place(xyz)` re-fits to new coordinates; `seam_segments`), `FILL_MODES`, `placed_ca` |
| `graft_calibration.py` | `calibrate(host, others, pred)` → `Trial`s: stretches others miss, hidden in the host, filled, scored vs a straight line and a global fit; `candidate_stretches`, `hide`, `trial` |
| `states.py` | `state_panel(paralog)` → `StateRow`s: every human deposit measured the same way (the gating transition at the pore) |

## `ip3r/physics/` — simulation

| File | Key names |
|---|---|
| `anm.py` | `build_hessian` (inverse-square springs), `ANM.calc_modes` (drops 6 × components), `ANM.label_symmetry` (C4 irreps A/B/E), `apply_generator`, `ModeSet` (`collectivity` κ, `is_collective`, `first(irrep)` skips local artefacts), `tetramer_sites`, `atom_displacements` |
| `network_checks.py` | the network against its own choices: `local_modes` → `LocalMode` (residue a low-κ mode sits on, the unresolved stretch between its sampled neighbours), `describe_local`, `cutoff_scan` → `CutoffRow` (lowest collective A, collective A together, cumulative, per cutoff of the registered grid `cutoff_grid`), `rmsip`, `stride_agreement` (RMSIP vs stride 1) |
| `transition_modes.py` | `transition_overlap(tr, reference, cutoff=)` → `TransitionOverlap` (overlap, cumulative, symmetry-matched null, irrep fractions, `report()`), `remove_rigid_body`, `irrep_fractions`, `null_cumulative` |
| `gating.py` | De Young–Keizer / Li–Rinzel: `m_inf`, `n_inf`, `q2`, `h_inf`, `tau_h`, `open_probability`, `bell_peak`, `hill_fit_left_flank`, `bell_at` (via `bell`), `GatingParams` |
| `gating_mak.py` | Mak et al. 1998 Hill-type steady state: `MakParams`, `k_inh` (Eq. 2), `open_probability` (Eq. 1), `bell_at`, `compare_flanks` (both models, one ruler) |
| `bell.py` | model-agnostic bell measurement: `measure_bell(f)` → `Bell` (peak, half-activation, half-inhibition, `width_decades`), `flank_shifts` |
| `calcium.py` | closed-cell Li–Rinzel: `simulate` → `Trace`, `fluxes`, `oscillation_metrics` (sustained only), `oscillation_window` (0.36–0.63 µM measured), `steady_state`, `CellParams` |
| `puffs.py` | stochastic DYK cluster: `simulate_cluster` → `PuffTrace` (`n_open` snapshot per `puff.record_dt` bin, `n_peak` most open within it, `peaks`), `detect_events` (on peaks), `fano`, `PuffParams` |
| `park_drive.py` | Siekmann/Cao park/drive receptor: `ParkDriveParams` (the 41 `pd.*` constants), `ip3_functions`, `gate_inf`, `mode_rates`, `constant_generator`, `stationary` (detailed balance), `open_probability`, `park_fraction`, `drive_open_probability`, `bell_at`; `STATES`/`OPEN`/`PARK` |
| `puffs_pd.py` | park/drive cluster: `simulate_cluster_pd` (split step: exact `expm` within modes, then mode switch; gate equilibria tabulated per number open; `clamp_ca` for the single-channel condition), `ParkDrivePuffParams` |
| `permeation.py` | 1-D drift-diffusion (ported from PIEZO1): `IonSpecies`, `potassium_species` (symmetric KCl; sweep overrides), `solve_pnp(z, r_free, fixed_charge=)` → `PermeationResult` (ohmic closure uncharged, local electroneutrality charged; Hall access), `series_conductance` (closed-form check), `debye_length`, `blocking_mechanisms` (steric only) |
| `_pnp_kernels.py` | the discretisation, ported unchanged: Scharfetter–Gummel `_nernst_planck`, `_ohmic_potential`, `_donnan_potential`, `_neutrality_step`, row-scaled Dirichlet solve |
| `pore_charge.py` | wall charge from the deposit's side-chain atoms: `charged_groups` (charge centre within `pore_charge.lining_margin` of the lumen; stubbed residues → `unplaced`), `map_charge` (Gaussian, charge-conserving), `pore_charge(pair_bridges=)` → `PoreCharge` (`bridged`: the pairs dropped) |
| `ryr_gating.py` | RyR1 by cytosolic Ca²⁺: `SternParams` (Stern 1997 two gates in series, `k_a`, `k_i`), `generator`, `stationary` (null space), `open_probability` (closed form), `bell_at`, `with_constants` (Ka, Ki, rate scale), `fit_to_bell` (off rates moved to both measured flanks); Mg²⁺ (`mg`, `k_mg_a` competitive at activation, `mg_i` equal to Ca²⁺ at inactivation; `k_a_eff`, `with_mg`); a scheme carries its own tables (`dest`, `open_mask`, `inact_mask`, `trigger_map`) and `exit_rates`/`generator`/`open_probability`, so every simulator runs any scheme; `MurayamaParams` (`at_37`), `murayama_activity` (Murayama 2015 Eq. 1), `murayama_bell`, `compare_bells`; `STATES`, `OPEN` |
| `ryr_two_site.py` | Stern's C channel with a two-site inactivation gate (inactivated only with both Ca²⁺ bound; six states, `s = a + 2n`): `TwoSiteParams` (a `SternParams`: `k_inact2_off`, `k_i2`, `uninactivated`), `two_site(k_a, k1, k2)`, `fit_two_site` (both half-peak points and the bell's log-log slope at half inhibition), `bell_slope`, `gate_hill_slope` |
| `sparks.py` | RyR1 cluster: `SparkParams`, `diffusion_coupling` (point source at the channel spacing), `simulate_sparks` (exact step tabulated per number open; returns a `PuffTrace` with `n_inactivated`), `spark_couplings(base=)` |
| `cleft.py` | Stern 1997's junctional cleft: `CleftGeometry` (`from_parameters`), `couplon` (two rows, V/C chessboard; `c_sites`, `v_sites`), `v_coupling_matrix` (µM at each C channel per open V channel), `edge_transfer` (Eq. 13), `coupling_matrix` (finite-volume steady solve, µM at each C channel per open C channel; memoised per geometry), `field` (the cleft map), `nearest_coupling` |
| `sparks_cleft.py` | RyR1 array in the cleft: `CleftSparkParams` (`ca_per_open` = nearest-neighbour coupling), `couplings_for` (off-diagonal scaled, own release kept), `c_rates` and `initial_states` (shared with the couplon), `simulate_sparks_cleft` (Gillespie, exact; a `PuffTrace`; `trigger` opens the available channels at t = 0), `native_coupling`, `DEST` |
| `spark_termination.py` | what ends a cleft spark: `Termination` (sparks, median duration, `unterminated`, inactivated at end, open fraction), `measure(sp)`, `refit` (Stern vs fitted to Murayama 25/37 °C), `ki_scan`, `rate_scan(base)` (Ki fixed, rates scaled) |
| `spark_mg.py` | Mg²⁺ and the cleft spark: `triggered(sp)` → `Triggered` (available channels opened at t = 0, timed to all shut; opened, ended, inactivated before/after), `dissect` (no Mg²⁺ / activation site / inactivation site / both at `ryr.mg_free`), `mg_scan`, `spontaneous`, `k_mg_a_by_ratio` (Laver 2004's selectivity on the fitted Ka), `READINGS` + `k_mg_a_reading` (the three K_Mg,A readings: measured, selectivity, meissner; shared by CLI and GUI) |
| `allosteric_v.py` | the V channel, Rios 1993's allosteric model: `RiosParams` (fiber 827; `rate_scale` = Stern's ×2), `generator(V)` (C0–C4, O0–O4), `exits` (the table the couplon draws from), `stationary`, `open_probability` (Eq. 6), `charge` (Eq. 7), `step_response` (master equation, Stern's Fig. 11 check) |
| `couplon.py` | the whole couplon: `CouplonParams` (C array + V model + V current), `step_protocol`, `v_trajectory` (one V channel through a voltage protocol, exact), `simulate_couplon` (V channels first, then the C array by Gillespie with V openings as scheduled events; `scale=(times, fractions)` makes every unitary current a piecewise-constant fraction, the SR content) → `CouplonTrace` (`flux()`, `scale`) |
| `lumen.py` | SR depletion (Stern's Fig. 20 pool): `LumenParams`, `mm_per_pa_ms` (couplon density → mM per pA·ms), `content_path` (release out, first-order uptake), `depleted_ensemble` → `Depleted` (path consistent with the ensemble's own release by damped fixed-point iteration; `released`, `peak_to_plateau`, `corrected` = Schneider's correction, `corrected_peak_to_plateau`, `full` = same couplons never depleting), `fig20_couplon` (28 channels), `depletion_panel`, `pool_scan` (pool × factor at Rios's +20 mV); `po_peak_to_plateau` |
| `ec_release.py` | release under voltage clamp: `ensemble(v, scale=)` → `Ensemble`, `summarise` → `Summary` (peak, plateau, after repolarisation, flux peak/plateau), `configurations(reading, two_site=)` (Stern / fitted / fitted + Mg²⁺ at the activation site / both sites; `two_site` uses `fit_two_site`), `with_gating`, `events`, `event_stats` → `EventStats` (single C events at `ec.spark_voltage`), `voltages` |
| `mg_competition.py` | Meissner 1997's Eq. 4 with two competitors at the activation site: `hill_ka` (one competitor), `ka_shift(mg, na=)` (Mg²⁺'s shift of half-activation in Na⁺; default Murayama's 0.17 M), `equivalent_k_mg_a(mg=)` (the K_Mg,A giving that shift in this scheme's form; 769 µM at 1 mM) |
| `ryr_mutants.py` | RyR1 charge mutants vs Xu 2006: `mutants()` (from the registered `permeation.published_ryr1_*`), `mutant_panel(st)` → wild type + `MutantRow`s (measured vs modelled ratio, charged and paired; lining, bridged), `open_deposit`, `open_mutant_panel` |
| `salt_bridges.py` | `salt_bridges(st, cutoff)` → `Bridge`s: Barlow & Thornton ion pairs (charged N–O ≤ `pore_charge.salt_bridge_cutoff`), matched one-to-one closest first over the whole deposit |
| `unitary.py` | `unitary(st, neutralise=)` → `Unitary` (series / neutral / charged / paired = salt bridges cancelled; `n.c.` where unconverged; bath by family, `bath_for`), `published(paralog)`, `unitary_panel(paralog, sweep=)` (the S11 state panel), `sensitivity` (diffusivity × ion-radius corners), `published` (Mak 2000, Vais 2010) |
| `puff_compare.py` | one ruler for all clusters: `MODELS` (the IP3R pair), `ALL_MODELS` (+ `SPARKS`: `SPARK` = `ryr1`, `SPARK_CLEFT` = `ryr1-cleft`, `SPARK_FIT` = `ryr1-cleft-fit`, `fitted_cleft_params`), `spark_scan(model=)`, `spark_ends` (duration, inactivated at start/end, `unterminated`), `MODEL_LABELS`, `params_for` (`mg=`, `k_mg_a=` for a spark receptor; IP3R refuses), `simulate`, `event_sizes`, `recruitment` (Fano, open fraction, blips / multi / large ≥ half the cluster, rate), `coupling_effect(model=)`, `scan_couplings`, `coupling_scan` |

## `ip3r/analysis/` — the findings checks

| File | Key names |
|---|---|
| `checks.py` | framework: `register` decorator, `Check`, `Outcome` (`confirmed/discrepancy/not_run/error`), `run_check(s)`, `all_checks`, `PAPERS`; refuses to confirm against a modified registry |
| `checks_structure.py` | `S0.*` (C4, filter, gate, pore profile, IP3 contacts — recomputed from coordinates) and `P6.shell_agreement`, `P6.contacts_heavy_atom` |
| `checks_constraint.py` | `P5.*` (element means, rankings, gate identity, variant AUCs under S17's position rules, deep-ranks-third, VUS count, ω) and `P6.contacts_vs_core` |
| `checks_evolution.py` | `P2.sister_pair` (own Newick reader), `P2.au_test`, `P2.teleost_itpr1`, `P3.*`, `P4.unreachable`, `LEDGER.claims` |
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
| `range_table.py` | Paper 1's S20/S23 rows: `load_proteomes` → `Proteome` (clade from taxonomy.tsv), `clade_rows` → `CladeRow`, `count_at(rank, name)`, `taxon_calls` (from the 6 assignment tables), `absence_targets` (S23 rule G3, `range.absence_min_proteomes`), `genome_absences` → `GenomeAbsence` (rebuilt from manifest + control ledger + copies + copy ledger; `CONTROLLED`), `copy_numbers`, `substantial_table` (relaxed hits at the registered bar; `LINEAGES`), `load_range` → `Range` (`absences_in(clade)`) |
| `range_figure.py` | `draw_range(ax, clades, absences)` — clade bars on a fixed 0–1 scale, supergroup colours (`GROUP_COLORS`), prokaryotes collapsed, genome absences crossed; exhibits `draw_presence`, `draw_relaxed`, `draw_absences`, `draw_copies`, `draw_chase` |
| `checks_range.py` | `P1.presence_range`, `P1.kingdom_absences`, `P1.relaxed_controls`, `P1.absence_targets`, `P1.absences`, `P1.copy_number`, `P1.record_chase` (all rederived) |
| `vus_strata.py` | Paper 5 §8: `stratify(rows, score, gene, layer, sources=)` → `Stratification` (per-position scores per class, both medians, `strata`, `row()` = S17's fields), `stratum_of`, `STRATA`, `CLASS_COLORS`/`STRATUM_COLORS` |
| `vus_figure.py` | `draw_fractions` (the check's exhibit: VUS shares per gene × layer, fixed 0–1), `draw_strip` (one gene, one layer: classes and medians) |
| `checks_variants.py` | `P5.vus_stratification` (every row of vus_stratification.tsv rebuilt; also ClinVar-only), `table_score` |
| `exhibits.py` | `draw(ax, check_id, outcome)` — figures from a check's own numbers |

## `ip3r/render/` (moderngl, OpenGL 4.1)

`camera.py`, `primitives.py`, `scene.py`, `spline.py`, `geometry_builders.py`,
`shaders/` — ported unchanged from PIEZO1 (impostor spheres/cylinders,
cartoon sweeps, trackball camera). `colormaps.py` — chain, element, fixed
conservation ramp (0.50–0.95 JSD; grey = not scored), fixed displacement ramp
(0–25 Å), `SHELL_COLORS`/`shell_colors` (S22's four shells, grey beyond), `PLDDT_COLORS`/`plddt_colors` (AlphaFold's fixed bands), `SEAM_COLORS`. `representations.py` —
`variant_spheres.py` — `variant_spheres(st, gene, classes, layer, chain_mask)` (Cα of each variant residue on every visible subunit; most decisive class wins; VUS by stratum from the resources), `resource_stratification`, `variant_classes`. `MolecularView` (styles × `ColorBy`, highlight (uniform or per-atom `highlight_rgb`), chain filter, `update_coords`
for animation; `ColorBy.DISPLACEMENT` from a built transition; `ColorBy.LIGAND_SHELL` from `structure.shells`; `ColorBy.PLDDT` for fills only).

## `ip3r/ui/` (PyQt6)

| File | Purpose |
|---|---|
| `app.py` | `main()` — surface format, theme, window, initial load (or `--session FILE`) |
| `main_window.py` | layout and wiring; menus (File: open/save session); loads on workers; `CHECK_SITES` maps a check to what "Show on structure" highlights, `CHECK_COLOURS` to a colouring (the shell checks), `CHECK_TREE` to the Tree tab, `CHECK_GENOMES` to the Genomes tab and its layer, `CHECK_RANGE` to the Range tab |
| `scene_controller.py` | what the viewport draws: `MolecularView`, pore spheres, site/variant highlights, `show_fill` (the `FillOverlay`), `show_variants` (sphere batch, refused in another numbering), `move_overlays` (spheres follow morph/mode frames), side/top views, mode animation |
| `fill_overlay.py` | `FillOverlay`: the AlphaFold fill as its own `MolecularView` (`fill:*`, always pLDDT) plus seam bonds (`seams`, red when broken); `move` re-places it on a morph/mode frame |
| `fill_controller.py` | `FillController`: the Completeness selector — builds on a worker, latest request wins, rebuilt on every load, refusal shown in the panel; `fill_html` |
| `gl_widget.py` | `ViewportWidget` (ported; viewport sized from the bound FBO every frame) |
| `structure_panel.py` | deposition list, style, colour, layer, Completeness (+ fill summary and pLDDT/seam legend), subunits, measured sites, legend |
| `channel_panel.py` | `ChannelSummary` text, pore profile vs S0's, state comparison of the loaded family (`set_paralog`: ITPR3 or RYR1), `show_mutants` (RyR1 only), `show_unitary` (conductance per state vs the measured values; `unitary_rows` for the smoke test) |
| `modes_panel.py` | ANM table with irreps and κ (local modes named by residue), animation controls |
| `transition_panel.py` | Transition tab: end state, fit, method, per-family preset (`preset_for`: 8TKG→8TKF, or RyR1's `morph_start`→`morph_end`), frame slider/play (label shows the frame's gate), displacement colouring, element, overlap and gate-along-the-path plots |
| `transition_controller.py` | `build_transition` (worker; also the `AtomPath` and `GatePath`), `TransitionController` (install, `coords_at` — matched atoms interpolated, the rest ride their Cα, `show_frame`, `play`, `reset` — path built from the displayed structure) |
| `dynamics_panel.py` | Gating (bell; model: DYK, Mak 1998 with the flank comparison, or RyR1: Stern scheme vs Murayama bell, plus the scheme fitted to it and that fit under the fibre's Mg²⁺ at both K_Mg,A readings), Oscillations (+ window scan), and the Puffs sub-tab |
| `puffs_panel.py` | `PuffsPanel`: receptor (DYK / park-drive / RyR1 sparks, mean-field, in the cleft, or in the cleft fitted to Murayama's bell; RyR ignores IP3), coupled vs uncoupled traces, event-size histogram ("no events" when silent), "Scan coupling (both receptors)"; RyR1: free Mg²⁺ and the K_Mg,A reading, "Triggered sparks vs Mg²⁺" (cleft receptors, `spark_mg.mg_scan`); `result` for the smoke test |
| `findings_panel.py` | checks by paper, run on a worker, claim/method/verdict, exhibit, show on structure (`showable`: residue-keyed checks drawn on the displayed structure) |
| `tree_panel.py` | Tree tab: `draw_tree` on a toolbar canvas, tip labels / support toggles, "Vertebrates" zoom, click names a tip; loaded on a worker when first shown |
| `genomes_panel.py` | Genomes tab: `draw_grid` with layer / sort / class / above-bar controls, click names a genome; loaded on a worker when first shown; `show_layer` (from a check's "Show") |
| `range_panel.py` | Range tab: `draw_range` with min-proteomes / collapse-prokaryotes / genome-absence controls, click lists a clade's genome absences; loaded on a worker when first shown |
| `variants_panel.py` | S17 variants per paralog/class; "Draw on structure" (spheres, `draw_requested`), "VUS by layer" (stratum column, `draw_strip` plot, stratum colours); `follow(paralog)` on load; drawn only in matching numbering |
| `params_dialog.py` | `ParametersDialog`: the registry editor (filter, modified-only, edit with clamp report, reset selected/all, import/export in the `IP3R_PARAMETERS` format; `edit(key, text)` for the smoke test) |
| `params_banner.py` | `ParametersBanner`: the amber "parameters modified" strip, driven by the registry's listeners; hosted in a full-width toolbar (`MainWindow.params_strip`) |
| `session_controller.py` | `SessionController`: File → Save/Open session, `--session`; `capture`, `save_to`, `apply(session, parameters=)` (asks when the saved parameter set differs; applying re-measures), finished from `loaded` and `transition_built` because both the load and the morph run on workers; a restore is dropped if another deposit arrives first |
| `plot_canvas.py`, `workers.py`, `theme.py` | helpers |

## `ip3r/resources/` (committed)

`parameters.json`, `references.json` (built by `scripts/build_parameters.py`);
`sequences.json`, `domains.json`, `sites.json`, `constraint.json`,
`variants.json`, `structures.json` (built by `scripts/sync_genes.py` from
ip3r_genes, each with the source paths, SHA-256 and ip3r_genes commit).

## `scripts/`

| File | Purpose |
|---|---|
| `parameter_table.py` (+ `parameter_table_pd.py`, the park/drive constants; `parameter_table_perm.py`, permeation and wall charge; `parameter_table_graft.py`, AlphaFold fills; `parameter_table_ryr.py`, RyR1 bath, measured conductances, gating and spark constants; `parameter_table_ec.py`, the V channel and the voltage protocol), `param_entry.py` (the shared entry constructor), `reference_table.py`, `build_parameters.py` | the registry and its provenance gate (duplicate reference keys fail the build) |
| `sync_genes.py` | import resources from ip3r_genes; `--check` reports drift |
| `curate_ryr.py` | the RyR1 resource from UniProt / InterPro / RCSB: `sequence`, `domains`, `panel` (rules 1–6: `verdict`, `state_of`), provenance hashes; `make ryr` |
| `screenshot_app.py` | scripted GUI smoke test + README screenshots (ends by saving a session, loading elsewhere, restoring, and comparing every field) |
| `screenshot_sparks.py` | its RyR1 spark steps (`spark_step`, `SPARK_STEPS`): mean-field, cleft, fitted, fitted under 1 mM Mg²⁺ (silent), the triggered Mg²⁺ scan |
| `create_env.sh` | the `ip3r_sim` conda env |

## `tests/`

Morph gate (`test_morph_pore` — both ends atom-for-atom, endpoint frames measure as the deposits for ITPR3 and RyR1, a flipped side chain's chord, rigid side chains caught missing 8TKF's gate), fills (`test_graft` — a rigidly moved chain fills exactly and follows the deposit, a one-off numbering / stubbed flanks / a bent loop / a clashing neighbour each caught, 8TKG's 56 stretches, 9YKK/7LHF refused, fills beat a straight line on 16 hidden stretches), transition (`test_transition` synthetic calibrations; `test_transition_real`
— the drawn end must be 8TKF as a shape, with a case that must fail), physics (`test_anm`, `test_network_checks` — a planted flap after a gap is located and a solid slab has no local modes, 8TKG's #11–15 are the residue-86 flap, the single A mode moves with the cutoff while the collective A modes hold, stride 2 agrees with stride 1 and stride 4 does not, `test_gating`, `test_gating_mak` — the paper's constants, the plateau, a planted K_act dependence the flank test must catch, `test_calcium`, `test_puffs`, `test_park_drive` — stationary = generator null space, printed IP3 functions reproduced; `test_puffs_pd` — the split step reproduces the clamped stationary P_open to 3 %, and a 5 ms step must fail; `test_puff_compare` — recruitment counted by hand, only park/drive recruits the cluster; `test_permeation` — cylinder and Donnan closed forms, solver = series sum, charge conserved, stubs counted, only 8TKF conducts and falls short of both measurements; `test_spline` — frames orthonormal and rotation-minimising (a plane curve keeps its normal); `test_cleft` — release leaves across the edges exactly, grid converged, Eq. 13 by hand, couplings beside Stern's Fig. 9, chessboard neighbours, the Gillespie array samples P_open, cleft sparks end by inactivation; `test_spark_termination` — the fit lands on both flanks and recovers Stern from Stern's bell, rate scaling leaves the steady state alone, an unended spark is flagged, fitted sparks never end and none start at 37 °C, duration grows with Ki; `test_ryr_two_site` — exits follow the state code, null space = closed form with and without Mg²⁺, the slope ruler on a Hill function, K2/K1 spans slope 1–2, the fit meets both points and the slope (the one-site fit misses the slope), the steeper gate inactivates less at cleft Ca²⁺; `test_mg_competition` — Meissner's Table IV competition predicts Table II's Ka in three salts within 1.6× and no competition misses NaCl 3×, Na⁺ blunts the Mg²⁺ shift, the equivalent reproduces it; `test_spark_mg` — Mg²⁺ 0 leaves the scheme unchanged, closed form = null space under Mg²⁺, Ka moves by the competitive factor, Mg²⁺ inactivates as Ca²⁺ at rest, the trigger opens exactly the available channels, activation-site Mg²⁺ alone shuts every triggered spark uninactivated, the GUI's route and readings build the same scheme; `test_allosteric_v` — null space = Eqs. 6–7, detailed balance, the paper's 0.72 ceiling, fiber 827 gives Stern's plateaus, Stern's Fig. 11 follows the printed rates and Fig. 12 twice them; `test_lumen` — the conversion by hand, the content path's closed forms, a path of ones is the same couplon event for event, an empty store releases nothing, the path is self-consistent, Stern's constants reproduce Fig. 20 (content and corrected plateau), the record grid is fine enough; `test_couplon` — Monte Carlo = master equation through the step and after it, V→C coupling reciprocal, no V current no trigger, Stern's constants reproduce his Fig. 12 with control, the fitted scheme loses control, both-site Mg²⁺ keeps it with no peak; `test_ryr_gating` — null space = product of gates, the printed k_i typo would remove inhibition, coupling by hand, clamped cluster = stationary P_open, step converges and 5 ms is caught, coupling makes sparks; `test_ryr` — the resource's literature anchors (GGGIGD 4894, I4937), each curation rule rejects its own violation, the panel gates at I4937, 9HEO short of 801 pS, D4899Q refutes pairing; `test_salt_bridges` — one Arg cancels one Asp, the cutoff is N–O, a bridged lining group is dropped without adding its partner, 8TKF's D2478–R2471′ bridges), geometry
(`test_symmetry`, `test_pore`, `test_structures_real`), statistics
(`test_stats` — Fisher vs scipy and the tea-tasting value, logistic slope = log OR for a binary predictor, Wilson vs Newcombe; `test_newick`), grid (`test_genome_grid` — the channel rule, the bar, ordering, the real grid's counts), alignment (`test_pairwise` — score equals a cell-by-cell reference DP), shells (`test_shells`), VUS strata (`test_vus_strata` — the rule by hand, ties, a residue in two classes, the resource route = S17's table, one sphere per subunit), tree (`test_tree` — toy trees with known clades, root twin edge, a misplaced tip is foreign), modules (`test_modules` — spans hold their sites, column map lands on the residue, a tampered reference is refused), provenance (`test_parameters` — listeners fire once per effective change, import is replace-not-merge and refuses a bad file untouched, a memoised measurement made under an edit is dropped on reset,
`test_resources`), sessions (`test_session` — round trip, the field set pinned so a result cannot ride along, malformed files refused by name, `replace` semantics), rules (`test_sizes`), CLI, and
**`test_checks_calibration.py`** — every check flipped by a planted input
(`PLANTS`), sources proven complete, `not_run` without data, refusal under
modified parameters.
