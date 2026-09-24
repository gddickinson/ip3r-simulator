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

## Two states: transition, morph, overlap

`structure/transition.py` reduces two deposits of one paralog to a common
basis. Both deposits must be in that paralog's human numbering. A residue
enters the basis only if it is resolved on all eight chains, is not stubbed,
and has the amino acid the reference gives at that number. Subunits are
ordered right-handed about each deposit's cytosol-up axis. All four cyclic
correspondences are tried; for a C4 pair they fit identically. The end is
superposed onto the start *as deposited*, so the path can be drawn over what
is on screen. `structure/morph.py` interpolates linearly, then restores
peptide Cα–Cα distances (the PIEZO1 method). Side chains ride their Cα.

`physics/transition_modes.py` scores the elastic network of one state against
the observed move (overlap = |cos|, Tama & Sanejouand 2001). It removes the
rigid-body part exactly, so the result does not depend on the fit. It splits
the move into C4 isotypic components, and compares each overlap with a
random direction **of the same irrep make-up**.

8TKG → 8TKF (stride 3, pore fit), measured 2026-09-23:

| quantity | value |
|---|---|
| basis | 2,194 residues × 4 subunits |
| RMSD, pore fit / overall | 3.14 / 16.11 Å |
| mean displacement: RIH_N, MIR, RIH_C | 22.0, 19.0, 16.9 Å |
| mean displacement: channel, gate, filter | 2.5, 3.4, 0.8 Å |
| A-symmetric share of the move | 100.0 % (C4 imposed in the maps) |
| lowest collective A mode (#5), overlap | 0.415 (null 0.021) |
| cumulative, 20 modes | 0.638 (null 0.048) |
| same, network of 8TKF scoring the reverse | lowest A 0.169; 20 modes 0.708 |

Across strides 1–4 the lowest collective A mode overlaps at 0.39–0.49 and
20 modes at 0.61–0.67. At stride 5 the network falls apart (0.03, at the
null). The 8TKG network therefore points towards activation, and one A mode
carries two-thirds of what 20 modes capture. From the other end the picture
is not symmetric: the lowest A mode of 8TKF points back only weakly (0.17).
Other pairs overlap less: 6DQJ → 8TKF reaches 0.32 over 20 modes, and
8TKH → 8TKF 0.16.

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

**Limitations.** The cluster Ca²⁺ is mean-field and instantaneous (Cao
integrates a microdomain ODE with fluo-4), and the store is never depleted.
Above about 0.5 µM coupling, the park/drive cluster settles into sustained
partial activity (9 % open) rather than discrete puffs.

## Unitary conductance (drift-diffusion over the pore)

Ported from PIEZO1's `physics/permeation.py`. For each species *i*, the
steady Nernst–Planck flux down a pore of accessible area
A_i(z) = π (r_free(z) − a_i)² is

  J_i = −D_i A_i [dc_i/dz + (z_i F / RT) c_i dφ/dz],  dJ_i/dz = 0,

Scharfetter–Gummel discretised. The potential is closed in the
electroneutral limit, because full Poisson coupling diverges when the Debye
length (5.8 Å here) exceeds the pore radius. Without wall charge, it is
closed by ohmic current continuity ∇·(σA∇φ) = 0. With charge, by local
electroneutrality Σ z_i c_i + X = 0 (a local Donnan partition). Hall's
access resistance 1/(4σa) is added at each mouth. The closed-form sum
R = ∫dz/(σA) + 2R_access is derived without the solver as its check.

**Inputs.** The profile is `r_free` of protein heavy atoms (not S0's
HETATM-inclusive one) over the pore-domain span ± 12 Å. The bath is
symmetric 140 mM KCl at room temperature, the condition both measurements
were made in. The fixed charge comes from the deposit's own side chains:
the carboxylate midpoint (Asp/Glu), NZ (Lys) or CZ (Arg), admitted when
within `pore_charge.lining_margin` (3 Å) of the atom-centre radius at its
height. It is spread by a 3 Å Gaussian and divided by the lumen area.
Stubbed side chains are counted as unplaced, never guessed; there are none
in the ITPR3 panel.

**Calibrations** (`tests/test_permeation.py`): a cylinder equals the exact
ohmic conductance to 0.1 %; the solver equals the series sum on an
hourglass to 2 %; a uniformly charged cylinder equals the exact Donnan
partition to 1 %; the charge kernel conserves charge; a pore below the K+
radius is shut. On 8TKF the result is grid-converged to 4 % (step 1 → 0.25 Å).

**Measured** (`python -m ip3r unitary`). All six non-activated ITPR3 states
are sterically shut (r_free 0.25–1.03 Å, below K+'s 1.38 Å). Activated 8TKF
(r_free 3.08 Å at the filter) gives 64 pS (series), 65 pS (neutral solver)
and 33 pS with its lining charges. Measured, rat ITPR3 in symmetric 140 mM
KCl: 358 ± 8 pS (Mak et al. 2000, oocyte nuclei; Vais et al. 2010 cites it
as 370) and 545 ± 7 pS (Vais et al. 2010, DT40 nuclei).

- Over the unmeasured constants (in-pore diffusivity 0.25–1× bulk, ion
  radius 1–2 Å), the neutral model spans 25–150 pS and the charged
  13–75 pS. No corner reaches 358 pS; the model is 2.4× short at best.
  (PIEZO1's same model was 1.5× *high*.)
- The resistance is spread along the ~50 Å pore; the filter slice holds
  about a third of it.
- The wall charges are seven rings: E2398, K2482, D2478, D2518, D2522,
  K2529 (and R2524/E2532 in other states). Alternating-sign rings act as
  junctions in series, since each carrier must cross a zone where it is the
  excluded co-ion. So the charges **lower** the conductance. Acidic rings
  alone give 174 pS, and without K2482 it is 64 pS.
- The charged number is not robust. The lining margin moves it 7–58 pS
  (0–8 Å) and the smoothing moves it 16–151 pS (1.5–6 Å). The partition
  density peaks at 18.6 M, above the 10 M packing ceiling, and the model
  rectifies (26 pS at −20 mV against 32 at +20) where the measured I–V is
  linear. It is a point-ion continuum at its limit. D2478 is also
  salt-bridged (2.5 Å) to R2471 of the neighbouring subunit, a charge the
  lining rule does not see.

What survives: the gate is the only state change that opens a conducting
pathway. Taken as a neutral continuum, the only open deposit's pore is
too narrow or too long to carry the measured conductance.

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

## What the findings checks establish

| kind | what agreement means |
|---|---|
| `recomputed` | the published number follows from the deposited coordinates by a route sharing no code with the publication |
| `rederived` | the published output follows from the published inputs by independent arithmetic; says nothing about the inputs |
| `read` | the prose states what the table says |

A check is only trusted after it has been shown to flip on a planted change
(`tests/test_checks_calibration.py`).
