# SCIENCE — the models, and what each check establishes

Every constant below is a registered parameter (`python -m ip3r params`);
this page says what the models are, not what the numbers are.

## The receptor

A C4 homotetramer of ~2,700-residue subunits (human ITPR1 2,758, ITPR2 2,701,
ITPR3 2,671). Each subunit carries an IP3-binding core at its N-terminus (the
β-trefoil PF08709 and the armadillo/MIR region), large regulatory RIH
domains, and a six-helix pore domain (PF00520) whose selectivity filter
(GGGVGD) faces the ER lumen and whose gate (Phe/Ile at the cytosolic
bundle crossing) sits ~20 Å above it. IP3 binds ~70 Å off the pore axis, one
site per subunit, and is contacted by ten residues of one subunit (6DQN;
Paknejad & Hite 2018). Residue numbers throughout are canonical human
UniProt numbers of the named paralog.

## Structure measurement

- **Axis.** Each subunit is superposed on the next (Kabsch 1976); the axis of
  that rotation is the four-fold axis, and its angle must be 90°. `ip3r_genes`
  S0 used the normal of the subunit-centroid plane instead; the two agree on
  6DQN to 0.02°.
- **Pore.** In 1.5 Å slabs along the axis, the minimum distance from the
  axis to a heavy-atom centre (`r_min`, S0's quantity) and the same less the
  atom's van der Waals radius (`r_free`, what the viewer draws). The filter
  is the narrowest luminal point, the gate the narrowest cytosolic point of
  the PF00520 span.
- **Numbering.** A deposit is in a paralog's numbering when ≥ 95 % of its
  side-chain-modelled residues match that sequence at the same number
  (S24's rule); stubbed residues are excluded and mismatching segments
  reported.

## Elastic-network modes (ANM)

Cα beads, springs within 15 Å weighted (d0/d)² (Atilgan 2001; Yang 2009),
every residue all four subunits resolve, strided. The character of a mode
under the 90° generator labels it A (+1), B (−1) or E (0). Symmetric
stimuli — IP3 at all four sites, a four-fold pore opening — couple at first
order only to A modes. Modes give directions and relative stiffness, not
amplitudes or time scales; the animation is illustrative.

**Collectivity.** A mode's κ (Brüschweiler 1995) is the fraction of sites
that effectively move. In coarse networks, weakly attached fragments give
near-zero-eigenvalue modes with κ ≤ 0.11; collective modes have κ ≥ 0.27.
Modes below `anm.min_collectivity` (0.2) are reported as artefacts and are
never called "the lowest A mode". At stride 4 on 8TKG the naive lowest A
mode was such a fragment, with κ 0.01 and overlap 0.008.

**Where the local modes come from** (`physics/network_checks.py`,
`local_modes`: the residue carrying most of the mode, summed over subunits,
and any unresolved stretch between its sampled neighbours). On 8TKG at
stride 3 all five (#11–15, B E E A B) are one piece: residue 86, the first
after the unresolved loop 77–85, held by 5 springs where the network's
median site has 17. Four copies of one flap give an A, a B and an E pair.
At stride 4 the same flap takes modes 2–5; at strides 1 and 2 it keeps
enough springs and there are no local modes. RyR1 9R8O has its own
(residue 896, even at stride 1; residue 1988 at stride 2), so no stride
removes them everywhere and the κ guard stays. Bridging the flap with
sequence-neighbour springs was tried and does nothing at any strength (1,
10, 100 × γ): the flap swings as a body rather than stretching along its
chain. What does remove it is keeping each site's coordination (a cutoff of
15 Å × stride^⅓), but that changes the model, not only its sampling (below).

**The stride is not free.** Against the stride-1 network (the reference;
7.8 s), the RMSIP of the first 20 modes (Amadei et al. 1999) is 0.97 at
stride 2, 0.89 at stride 3 and 0.77 at stride 4: the flap's modes displace
collective ones from the computed window. The first 10 modes agree to 0.99
at strides 2–3.

## Two states: transition, morph, overlap

`structure/transition.py` reduces two deposits of one paralog to a common
basis. Both deposits must be in that paralog's human numbering. A residue
enters the basis only if it is resolved on all eight chains, is not stubbed,
and has the amino acid the reference gives at that number. Subunits are
ordered right-handed about each deposit's cytosol-up axis. All four cyclic
correspondences are tried; for a C4 pair they fit identically. The end is
superposed onto the start *as deposited*, so the path can be drawn over what
is on screen. `structure/morph.py` interpolates linearly, then restores
peptide Cα–Cα distances (the PIEZO1 method).

`structure/morph_pore.py` then interpolates every heavy atom both deposits
resolve (matched by chain, residue and atom name; 70,952 for 8TKG → 8TKF,
none unmatched). The atom's offset from its own Cα is interpolated linearly
while the Cα follows the morph, so frame 0 is the start deposit and the last
frame the superposed end deposit, atom for atom. The axis is re-found on
every frame and the gate located by the deposit rule, so the endpoint frames
measure exactly as the deposits do (8TKG 2.73 Å at F2513, 8TKF 5.85 Å at
N2510; 9R8O 3.32 Å at I4937, 9HEO 5.05 Å at Q4933). Side chains riding their
Cα rigidly — the viewer's shortcut until then — end 8TKF's path at 5.16 Å:
8TKF's backbone lined by 8TKG's rotamers, 0.69 Å short of 8TKF. On RyR1 the
shortcut errs the other way (5.24 Å against 5.05 Å).

Along the path the gate opens monotonically (no overshoot) and is half-way
at fraction 0.42 (ITPR3) and 0.44 (RyR1). In ITPR3 the constriction moves:
I2517 lines it from t ≈ 0.07 to 0.69, and N2510, 7 Å further down, from
0.72. Everything between the ends is an interpolation, so the half-way point
is a property of the path, not a statement about the gating order. A
swinging side chain takes a chord; `lining_error` reports the worst
shortening of an atom–Cα distance among the gate's lining atoms per frame.

`physics/transition_modes.py` scores the elastic network of one state against
the observed move (overlap = |cos|, Tama & Sanejouand 2001). It removes the
rigid-body part exactly, so the result does not depend on the fit. It splits
the move into C4 isotypic components, and compares each overlap with a
random direction **of the same irrep make-up**.

8TKG → 8TKF (stride 2, the default since 2026-09-24; pore fit):

| quantity | value |
|---|---|
| basis | 2,194 residues × 4 subunits |
| RMSD, pore fit / overall | 3.14 / 16.11 Å |
| mean displacement: RIH_N, MIR, RIH_C | 22.0, 19.0, 16.9 Å |
| mean displacement: channel, gate, filter | 2.5, 3.4, 0.8 Å |
| A-symmetric share of the move | 100.0 % (C4 imposed in the maps) |
| collective A modes together (5) | 0.669 |
| lowest collective A mode (#5); best single (#10) | 0.390 (null 0.017); 0.498 |
| cumulative, 20 modes | 0.670 (null 0.039) |
| same, network of 8TKF scoring the reverse | lowest A 0.195; 20 modes 0.585 |

At stride 3 (the default until 2026-09-24) the lowest A mode gave 0.415 and
20 modes 0.638. Across strides 1–4 the lowest collective A mode overlaps at
0.39–0.49 and 20 modes at 0.61–0.67. At stride 5 the network falls apart (0.03, at the
null). The 8TKG network therefore points towards activation. From the other
end the picture is not symmetric: the lowest A mode of 8TKF points back only
weakly (0.20). Other pairs overlap less: 6DQJ → 8TKF reaches 0.36 over 20
modes, and 8TKH → 8TKF 0.14 (stride 2).

**The single mode is a cutoff choice; the A subspace is the result**
(`network_checks.cutoff_scan`, `transition --cutoff-scan`, measured
2026-09-24). With no sampling artefact at all (stride 1), the lowest A
mode's overlap falls steadily as the cutoff rises, while the three lowest A
modes (#5, 9, 10) together hold:

| cutoff (Å) | 12 | 13.5 | 15 | 16.5 | 18 | 21.6 |
|---|---|---|---|---|---|---|
| lowest A mode | 0.48 | 0.45 | 0.43 | 0.38 | 0.33 | 0.24 |
| A modes in the first 10, together | 0.51 | 0.66 | 0.67 | 0.67 | 0.67 | 0.68 |

A longer cutoff redistributes the move among near-degenerate A modes. The
earlier sentence here, "one A mode carries two-thirds of what 20 modes
capture", held only at 15 Å. At stride 3 the scan is noisier (collective A
together 0.59–0.67 over 13.5–21 Å) because the flap comes and goes, and at
12 Å the network falls apart (every mode local). RyR1 9R8O → 9HEO is
different: its lowest A mode overlaps only 0.12, and that is flat in the
cutoff (0.120–0.127 at stride 1).

## Unresolved stretches, filled from AlphaFold

`structure/graft.py` places an AlphaFold DB model's residues where a deposit
has none. This is the PIEZO1 simulator's local-anchor method, applied to every
stretch separately. Each stretch is placed by a Kabsch fit of the prediction's
Cα onto the deposit's on up to `graft.anchor_window` residues each side
(`graft.min_anchor` per side, or twice that on a terminus's one side). An
anchor must be resolved, must not be a backbone+CB stub, and must be the same
amino acid in both models at that number. On 6DQN this skips 8 stretches per
subunit inside or beside its unregistered 1434–1546 segment, rather than
fitting them to a segment of unknown register.

**Which model.** The model is chosen by identity by number against the
variant-painting bar (`numbering.min_identity`). AlphaFold DB (queried
2026-09-24) holds canonical ITPR3 (`AF-Q14573-F1`), ITPR1 isoform 4 only, rat
ITPR1 isoform 8 only, and a 181-residue ITPR2 isoform. The seven ITPR3
deposits match at 99.1–100 %. 9YKK (17.7 % at best) and 7LHF (22.1 %) are
refused.

**Per fill.** The anchor RMSD. The Cα distance across each seam, which is
3.80 Å for a real peptide. The mean pLDDT and the fraction ≥ 70. The residues
with a heavy atom within `graft.clash_distance` of any deposited heavy atom,
own seam neighbours excepted and other subunits included. On 8TKG (gaps
only): 56 stretches and 1,592 residues, mean pLDDT 38. 12 of 224 seams are
broken, mostly at 1558–1585, where the anchors themselves disagree (5.3 Å).
108 filled residues clash with the deposit.

**Calibration** (`structure/graft_calibration.py`). A stretch that another
ITPR3 deposit leaves unresolved, but the host resolves entirely, is cut out
of the host, filled by the same route, and scored against the host's own Cα.
Two baselines use no prediction: a straight line between the flanks, and the
prediction superposed on the whole chain. Over 16 stretches in 8TKG and 15
in 8TKF, the fill lands at a median 1.46 Å, against 5.59 Å for the line and
10.5 Å for the global fit. It beats the line in 28 of 31. pLDDT predicts the
error (Spearman −0.79). The seams of these true fills reach 5.42 Å (90th
percentile 4.74 Å). The first tolerance, 4.5 Å, failed 9 of those 62 seams,
so `graft.join_tolerance` is 5.5 Å.

**The calibration's limit.** The tested stretches have pLDDT 53–77 and are
4–10 residues long. The real gaps are mostly longer and below pLDDT 50, and
there the calibration says nothing. A fill there is a picture of where a
chain of that length could run, not a structure.

**Nothing measures on it.** The fill is a separate structure drawn beside the
deposit. The pore, the modes, the transition and the checks all see the
deposit only. `FilledModel.place` re-fits every stretch on its anchors in
any coordinates, so a fill follows a morph frame without being part of it.

## Gating (De Young & Keizer 1992; Li & Rinzel 1994)

Per subunit: an IP3 site, a fast activating Ca²⁺ site and a slow inhibitory
Ca²⁺ site. With the fast sites at equilibrium,

    m∞ = p/(p+d1),  n∞ = c/(c+d5),  dh/dt = a2 (Q2 (1−h) − c h),
    Q2 = d2 (p+d1)/(p+d3),          P_open = (m∞ n∞ h)³.

The steady state is bell-shaped in Ca²⁺ (Bezprozvanny et al. 1991). More
IP3 raises the bell and moves its inhibitory flank out (IP3 relieves Ca²⁺
inhibition). In this model the activating flank moves too, by less (0.1 →
10 µM IP3: 2.4× vs 1.8×); Mak et al. (1998) measured IP3 tuning inhibition
alone, which the DYK scheme does not reproduce.

## Gating (Mak, McBride & Foskett 1998)

Single IP3R-1 channels in *Xenopus* oocyte nuclear patches, fitted by a
biphasic Hill equation with one denominator (their Eqs. 1–2):

    P_open = P_max / [1 + (K_act/c)^H_act + (c/K_inh)^H_inh],
    K_inh(p) = K_∞ / [1 + (K_IP3/p)^H_IP3].

P_max 0.81, K_act 0.21 µM, H_act 1.9, H_inh 3.9, K_∞ 52 µM, K_IP3 50 nM and
H_IP3 4 were each read from the paper (PMC28128) before they were
registered. "Kinh being the only IP3-concentration-sensitive parameter" is
the model's content. It describes steady-state data only, with no
inhibition kinetics, so it cannot drive the cell or puff models (the
park/drive receptor, below, is the kinetic model used for puffs).

Both models are measured with one ruler (`physics.bell`): the Ca²⁺ at half
the bell's own peak on each flank. From 33 nM (the lowest IP3 at which the
paper says activation was unaffected) to 10 µM:

| model | half-activation | half-inhibition |
|---|---|---|
| De Young–Keizer | 2.01× | 2.76× |
| Mak 1998 | 1.016× | 6.22× |

The test can fail. With IP3 dependence planted into K_act, the same ruler
sees the activating flank move (a test asserts it). Below K_IP3 the Mak bell
collapses (10 nM: peak 0.11, K_inh 0.08 µM < K_act), as the paper reports
for 10–20 nM. The Hill curve gives K_inh(33 nM) = 8.3 µM against the 9.5 µM
the paper measured at that point. The two models' IP3 sensitivities differ
by roughly tenfold: Mak's K_inh saturates by 0.1 µM IP3, while DYK's bell
is still rising at 10 µM.

## Cell Ca²⁺ (closed-cell Li–Rinzel)

Release through the receptor, an ER leak and SERCA uptake, with total Ca²⁺
conserved. Measured here by simulation: sustained oscillations for IP3 in
0.36–0.63 µM (period ~11–13 s); damped spirals are not counted.

## Puffs

**The DYK cluster** (`physics/puffs.py`). N receptors × 4 subunits × 3
two-state sites, flipping with the DYK rates on a fixed step (Shuai & Jung
2002). The receptors are coupled through a mean-field cluster Ca²⁺
(Swillens et al. 1999): every receptor sees `ca_rest + ca_per_open ×
(number open)`.

**The park/drive cluster** (`physics/park_drive.py`, `physics/puffs_pd.py`).
The same cluster with the receptor swapped for Siekmann et al.'s (2012)
six-state IP3R-1 model, with Cao et al.'s (2013) gating variables. A drive
mode (C1, C2, C3, O6) is open 70 % of the time; a park mode (C4, O5) is
almost never open. Ca²⁺ and IP3 act only on the switch between modes:

    q24 = a24 + V24 (1 − m24 h24),   q42 = a42 + V42 m42 h42,
    dG/dt = λ_G (G∞(c) − G)

`h42` recovers at 0.5 s⁻¹ while the channel is closed and falls at 20 s⁻¹
while it is open. An open receptor sees its own mouth (`c + 120 µM`), which
is Cao's two-concentration scheme. All 41 constants were read from the
authors' code (Cao et al. 2014, Text S1). The 2013 paper's Table S1 was not
reachable. Its printed IP3 dependences are registered as `base + amp × Hill`;
a test proves the two forms equal. Integration splits each 0.1 ms step: the
constant-rate transitions are taken exactly (`expm(Q dt)`), then the mode
switch and the gating variables are updated. Against the stationary
distribution at clamped Ca²⁺, the error is 0.2 % at 0.1 ms and 10.5 % at 5
ms. The test uses 3 % and has a case (5 ms) that must fail.

**One ruler** (`physics/puff_compare.py`). Both simulators record, per 1 ms
bin, the number open at the bin's start (for the Fano factor) and the most
open at once in the bin (for events). An event is a run of bins with any
channel open, and its size is its peak. The *recruitment* statistic counts
the events that reach half the cluster (10 of 20).

Measured over 30 s at 0.2 µM IP3, seed 0, both receptors at the same couplings:

| coupling (µM per open) | DYK Fano | DYK ≥10 | PD Fano | PD ≥10 |
|---|---|---|---|---|
| 0 | 0.93 | 0 | 0.98 | 0 |
| 0.09 | 1.32 | 1 | 2.79 | 10 |
| 0.17 | 1.21 | 0 | 2.78 | 17 |
| 0.32 | 1.14 | 0 | 2.23 | 16 |
| 1.08 | 1.25 | 1 | 1.60 | 6 |

The result holds on seeds 1–3: at 0.1 µM, park/drive has 12–14 half-cluster
events and DYK 0–1. Uncoupled, 5.2 % of DYK channels are open at a time,
against 1.7 % for park/drive. Park/drive's many "blips" are mostly brief
(0.3 ms) park-mode flickers.

The park/drive coupling (`puff.pd_ca_per_open` = 0.1 µM) was chosen by the
rule that chose DYK's 1 µM: the value that most raises the Fano factor.
Cao's own microdomain gives about 0.11 µM per open channel.

**Limitations.** The cluster Ca²⁺ is mean-field and instantaneous, and the
store is never depleted. Above about 0.5 µM coupling, the park/drive cluster
settles into sustained partial activity (9 % open) rather than discrete
puffs. `docs/SCIENCE_PUFF_DOMAIN.md` (Round 7.2) replaces the mean field with
Cao's microdomain, fluo-4 and a store that can deplete. It reproduces Cao's
IPI trend and the fluorescence bend at N ≈ 12, and shows the sustained state
belongs to the receptor model, not to the missing store.

## Permeation: conductance, selectivity and the Ca²⁺ current

The drift-diffusion pore model, its calibrations, the K⁺ conductance of
every ITPR3 state, and the selectivity and unitary Ca²⁺ current under
Vais 2010's protocols are in [`SCIENCE_PERM.md`](SCIENCE_PERM.md).

## Ryanodine receptors

RyR1 structures, conductance, charge mutants, gating and sparks are in
[`SCIENCE_RYR.md`](SCIENCE_RYR.md).

## Paper 6: the module contrast

The **ligand core** is the smallest span that holds all ten IP3 contacts. The
**pore module** is PF00520 less the luminal loop. Both are rebuilt in
`core.modules` from the imported sites and domain map, and validated against
what they must contain. For each orthologue in S17's deep alignment, the
identity to the human reference is computed in each module. Only columns the
tip covers are counted, and a tip must cover at least 50 % of each module
(`ligand.module_min_coverage`). The test is on the paired difference core −
pore: an exact sign test, and a signed-rank test that drops zeros, corrects
for ties and applies no continuity correction. The residue → column map comes
from walking the alignment's reference row, and is refused unless that row's
ungapped sequence is the UniProt sequence. S17's `deep_col` column is not
used.

What agreement establishes: from the alignments, the published spans, tip
counts, means, sign counts and p-values follow by independent code, and the
reversal with the loop counted in is real. It says nothing about the
alignments themselves.

## Paper 6: ligand shells

For each residue, the distance is the minimum all-atom distance to the IP3 on
its own subunit, with hydrogens included as in S22's reader. Per deposit it
is the best over the four subunits, and the consensus is the median over the
six IP3-bound ITPR3 depositions. Shells are half-open: [0, 4.5), [4.5, 8),
[8, 11.5) and [11.5, 15) Å. Residues beyond 15 Å belong to no shell.

S22 carried the pocket to ITPR1/ITPR2 with a MAFFT pairwise alignment. Here
it is carried with a separate alignment: Gotoh's affine-gap global
alignment, BLOSUM62, gap open 10 and extend 0.5 (EMBOSS needle's defaults),
end gaps free (`core/pairwise.py`). All 250 transferred positions agree with
S22's. Conservation is S17's per-residue deep JSD. Each shell is tested
against every scored residue with a one-sided Mann–Whitney (normal
approximation, tie and continuity corrections). The trend is Spearman's ρ
with a t-distribution p.

"No step at 4.5 Å" is prose, not a table, so it is tested two ways:

1. The contact shell does not beat the second (one-sided Mann–Whitney).
   With 12 and 14 residues this test is weak on its own.
2. The contact boundary does not carry the largest drop in mean JSD
   between adjacent shells.

It holds both ways in all three paralogs. The largest drop is at the 11.5 Å
boundary (third → fourth) everywhere. Note that the FEL purifying fraction
(the paper's §8) falls after the *second* shell instead, so the two
instruments put the drop in different places.

## Paper 2: the tree

The P2 tree checks ask their questions of the committed `rooted.nwk` with this
project's Newick reader (`analysis/newick.py`, `analysis/tree.py`). No
ip3r_genes clade table goes into the answer. A tip's group comes from the
census prefix on its label. Cyclostomes are recognised by genus (*Myxine*,
*Petromyzon* and the other hagfish and lamprey genera).

- **A paralog's whole clade** is the MRCA of every tip whose record names it.
  The check requires that clade to contain nothing foreign, meaning no other
  paralog, no cyclostome and no invertebrate. A second rule, "expand the
  named core while the clade stays pure", gives the same 19/13/19.
- **A cyclostome clade** is a maximal clade made only of cyclostome tips. A
  clade "branches first among the vertebrates" when its parent holds exactly
  the vertebrate tips.
- **Support is counted per bipartition.** In a rooted tree the root's two
  edges are a single bipartition of the unrooted tree, and IQ-TREE labels
  both of them. Counting nodes gave 92 of 132; counting bipartitions gives
  the published 91 of 131. The bars are the registered parameters
  `tree.alrt_min` (80, Guindon 2010) and `tree.ufboot_min` (95, Hoang 2018).

## Papers 3 and 4: the genome × paralog grid

The retention sweep asked 309 vertebrate assemblies for four cells each:
ITPR1, ITPR2, ITPR3 and a ryanodine-receptor control. `analysis/genome_grid.py`
joins three per-cell tables into one grid: what the search found
(`contiguity_cells.tsv`), the S15a evidence state (`character_matrix.tsv`)
and protein-record recovery (`gene_recovery.tsv`).

- **The contiguity bar** is the registered `genomes.contiguity_bar_bp`
  (142,212 bp, D4's median measured ITPR span). This project splits genomes
  on contig N50 ≥ bar with its own comparison, and the check tests that the
  split agrees with the table's `spans_gene` column in every cell.
- **The method control** is recomputed with `analysis/stats.py`: Wilson
  intervals, a two-sided Fisher exact test (hypergeometric sum, R's 1e-7
  tie tolerance), and a logistic regression of found on log10 N50 (IRLS,
  Wald p). The Mann–Whitney test of found against missed N50 is
  **one-sided**. The first version doubled it and was off by exactly a
  factor of 2, which was the checker's error.
- **The recovery channel** is rebuilt from the table's count columns in
  this order: gene present? reference proteome? any family record? any
  full-length record resolving to a cell? one resolving to *this* cell?
  It agrees with the table's own label in all 1,236 cells.

## Paper 1: presence and absence across eukaryotes

`analysis/range_table.py` joins the S20 proteome sweep (per-proteome
presence, the taxonomy, the per-record assignments, the relaxed search) and
the S23 genome sweep (manifest, control ledger, copies, copy-number ledger).
A clade is the taxonomy table's own `clade`. Seven checks re-derive the
paper's numbers from those rows:

- **Presence two ways.** Presence is counted from the presence table's call
  count, and again record by record from the six assignment tables. The two
  agree in all 6,854 taxa. They can only be compared per taxon, because 72
  taxa have more than one proteome.
- **Absence targets from the rule.** S23's rule G3 is re-implemented from its
  one-sentence statement: a eukaryotic phylum or class with ≥
  `range.absence_min_proteomes` (10) swept and no call. It gives exactly the
  35 clades `absence_at_genome.tsv` reports.
- **Genome absences from the ledgers.** Each target's genomes, controlled
  genomes (`controlled_cross_kingdom` or `controlled_by_target` only),
  complete-gene genomes and trace-only genomes are counted from the
  per-genome ledgers. All 35 rows match the published table.
- **Substantial matches** are full E ≤ `range.substantial_evalue` (1e-5)
  and model coverage ≥ `range.substantial_coverage` (0.5), with each target
  counted once per lineage. All 24 cells of S20's four-lineage table
  reproduce.
- **The record chase.** Rule R5 calls a record a fragment if it is below
  `range.family_floor_aa` (2,000) *or* UniProt flags it as a fragment. The
  first checker applied only the floor and got 50 real genes instead of 47:
  three records of 2,366–2,858 aa are UniProt-flagged fragments. That was
  the checker's error.

## Paper 5 §8: the VUS stratification

`analysis/vus_strata.py` rebuilds S17's rule from its description. For one
gene on one layer, each class (P/LP, B/LB, VUS) is scored once per
**position**. A residue scores only at occupancy ≥
`constraint.min_occupancy` (0.5). On the deep layer that is the same set as
the `deep_reliable` flag: the two agree at every residue of all three
paralogs. A residue carrying alleles of two classes counts in both. The
thresholds are the labelled sets' own medians. A VUS at or above the P/LP
median is pathogenic-like, and one at or below the B/LB median is
benign-like; ties count. `P5.vus_stratification` rebuilds all 12 rows of
`vus_stratification.tsv` from `variants.tsv` and the per-residue tables, and
every field agrees.

The stratification takes every source, so the 13 curated UniProt P/LP
records are in its P/LP median. The classifier test (`P5.variant_auc`) is
ClinVar only. Run ClinVar-only, the table is unchanged: every curated record
sits at a position ClinVar already labels P/LP, or at one no layer scores. So
the difference in rules changes no number here.

The viewer computes the same stratification from the committed resources
(`constraint.json`, `variants.json`), and a test proves that it reproduces
S17's table.

## What the findings checks establish

| kind | what agreement means |
|---|---|
| `recomputed` | the published number follows from the deposited coordinates by a route sharing no code with the publication |
| `rederived` | the published output follows from the published inputs by independent arithmetic; says nothing about the inputs |
| `read` | the prose states what the table says |

A check is only trusted after it has been shown to flip on a planted change
(`tests/test_checks_calibration.py`).
