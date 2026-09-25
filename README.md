# IP3R Structural Simulator

An interactive 3-D model of the **IP3 receptor** — the endoplasmic
reticulum's IP3- and Ca²⁺-gated Ca²⁺-release channel — driven by models with
a citable basis rather than by animation, and an instrument that
**illustrates, demonstrates and independently re-derives the results of the
[`ip3r_genes`](https://github.com/gddickinson/ip3r_genes) publication
project**. It is the IP3R counterpart of the
[PIEZO1 simulator](https://github.com/gddickinson/piezo1-simulator), whose
renderer, reader and conventions it reuses.

![6DQN coloured by functional element, IP3 contacts highlighted](docs/img/gui_element.png)

*Human ITPR3 with IP3 bound (PDB 6DQN), cytosolic cap up. Blue/cyan/green:
the IP3-binding β-trefoil, MIR and RIH domains; violet: the pore domain;
gold spheres: the ten residues within 4.5 Å of IP3 at each of the four
sites.*

## What it does

**See the channel.** Nine curated depositions — every ITPR reference and
gating-state structure `ip3r_genes` S11 selected, plus 6DQN — drawn as
cartoon, tube, spheres or sticks, coloured by functional element, subunit,
secondary structure, B-factor or **per-residue conservation** (the S17
tables, four layers, on a fixed scale with unscored residues grey). Residue
annotation is painted only when the deposit is verified to be in that
paralog's human numbering (rat 7LHF is not, and is left grey).

**Measure it.** The four-fold axis is found by superposing each subunit on
its neighbour; the pore profile, the selectivity filter (GGGVGD) and the gate
are located; every bound IP3 is found and its contacts listed.
`python -m ip3r states` measures the whole ITPR3 state panel the same way —
and shows the gating transition at the pore:

| PDB | state | gate radius (Å) |
|---|---|---|
| 8TKH | labile resting | 1.95 |
| 8TLA | higher-order inhibited | 2.44 |
| 7T3P | preactivated | 2.53 |
| 6DQN | IP3-bound | 2.55 |
| 6DQJ | apo | 2.69 |
| 8TKG | resting | 2.73 |
| **8TKF** | **activated** | **5.85** |

**Move it.** An elastic-network model of the tetramer gives its collective
modes, each labelled by its C4 irreducible representation: **A** (all four
subunits alike — the only kind that can couple to IP3 binding at all four
sites and to a symmetric pore opening), **B**, or the degenerate **E** pair.

**Watch it open.** The Transition tab (and `python -m ip3r transition`)
puts two states of one paralog on a residue-matched basis — 2,194 residues
on each of the four subunits of resting 8TKG and activated 8TKF — superposes
them on the pore domain, and morphs between them with peptide Cα–Cα
distances restrained (an interpolation, labelled as one). Side chains are
interpolated too, atom by atom, so the last frame is 8TKF itself and the gate
plotted along the path (2.73 → 5.85 Å, half-way at 0.42) is the gate on
screen; the rigid-side-chain shortcut, drawn dashed, ends 0.7 Å short. Each residue is
coloured by how far it moves, on a fixed 0–25 Å scale: the cytosolic RIH
and MIR domains move 17–22 Å on average, while the pore domain moves 2.5 Å
and the filter 0.8 Å. The tab then asks whether the resting
state's elastic network points towards the activated one. It does: the
collective A modes together overlap the observed displacement at 0.67,
where a random direction of the same symmetry scores 0.04 over all 20 modes.
No single mode is the answer. The lowest A mode (#5) gives 0.39 and #10
gives 0.50, and how the move splits among them depends on the 15 Å cutoff
(from 12 to 21.6 Å the lowest falls 0.48 → 0.24 while the A modes together
hold at 0.66–0.68; `transition --cutoff-scan`). So the tab's headline plot
is the A subspace: the collective A modes added lowest first, beside the
same-symmetry null and the √(A share) ceiling. The displacement is 100 %
A-symmetric, but that is inherited from C4-imposed reconstruction and is not
a finding. The network keeps every second Cα: at every third, the strand
after the unresolved 77–85 loop (residue 86) hangs on too few springs and
produces five spurious local modes. Any low-collectivity mode that remains is
named where it sits (the Modes tab and the report say so).

![The transition tab](docs/img/gui_transition.png)

**Gate it.** The De Young–Keizer receptor as reduced by Li & Rinzel: the
bell-shaped Ca²⁺ dependence of open probability, IP3 relieving Ca²⁺
inhibition, whole-cell Ca²⁺ oscillations (measured window 0.36–0.63 µM
IP3), and stochastic clusters in which Ca²⁺ coupling turns independent blips
into cooperative openings. The Gating panel also offers the Hill-type model
Mak, McBride & Foskett (1998) fitted to single IP3R-1 channels, in which IP3
tunes Ca²⁺ inhibition alone. From 33 nM to 10 µM IP3 it moves the
half-inhibition point 6.2× and half-activation 1.016×; De Young–Keizer moves
them 2.8× and 2.0× (`python -m ip3r gating --model mak`). The park/drive
receptor's stationary bell is there too (`--model pd`). It moves them 43×
and 1.09×, so it passes the same test. "The three models side by side"
draws all three bells at 33 nM and at 10 µM IP3.

![The three IP3R gating models](docs/img/gui_gating_compare.png)

![The Mak 1998 gating model](docs/img/gui_gating_mak.png)

The Puffs panel can fill the same cluster with park/drive receptors
(Siekmann et al. 2012, with the gating variables of Cao et al. 2013; every
constant read from the authors' code). Nearly all of these are parked at
rest, so the cluster stays quiet until one receptor enters drive mode. Then
its Ca²⁺ pulls the others in. Both receptors are measured with one ruler.
Over 30 s at 0.2 µM IP3, the De Young–Keizer cluster reaches half its
channels in at most 1 event at any coupling from 0 to 2 µM per open channel.
The park/drive cluster does so 10–17 times at 0.09–0.32 µM (Fano 2.8, against
≤ 1.32), and its event sizes split into blips and puffs with a valley
between (`python -m ip3r puffs --scan`).

![Park/drive puffs](docs/img/gui_puffs_pd.png)

`python -m ip3r microdomain` puts the park/drive cluster in Cao et al.'s
microdomain: Ca²⁺ pools that fill and drain, fluo-4, and a store that can
deplete. Puffs are then read from F/F0. As the h42 recovery rate rises from 0.1 to
5 s⁻¹, the puff rate rises 5.5× and the inter-puff-interval CV climbs from
0.79 toward 1, as in Cao 2013. With the release set so that the mean blip
matches Cao's, fluorescence amplitude bends at about 12 receptors while
Ca²⁺ bends less, because the dye saturates. The mean-field cluster's
sustained 9 % open state survives the microdomain's kinetics, and depleting
the store makes it more active, not less. That state belongs to the receptor
model (`docs/SCIENCE_PUFF_DOMAIN.md`). In the Puffs panel, the park/drive
receptor has a "Simulate in the microdomain" box. It draws F/F0 with the
puffs marked, the number open, and the inter-puff intervals against
Thurley's refractory density and an exponential with the same mean.

![The microdomain cluster](docs/img/gui_puffs_domain.png)

The Channel tab turns each ITPR3 deposit's pore into a K+ conductance by
drift-diffusion (ported from PIEZO1), with and without the charges of the
side chains that line it. Only activated 8TKF conducts. It gives 65 pS
uncharged (25–150 pS over the unmeasured diffusivity and ion radius)
against 358–545 pS measured, so the continuum model falls 2.4× short even
at its most generous. Its own lining charges lower the conductance rather
than raising it. Cancelling the salt-bridged ones (the filter's D2478 is
paired with R2471 of the next subunit) lowers it further, to 23 pS
(`python -m ip3r unitary`).

Most of that 2.4× is the model's geometry (`python -m ip3r shortfall`).
The 1-D model reduces each slice to its inscribed circle. The same
electrolyte in the voxelised lumen, solved in 3-D, conducts 1.3–2.0× more,
because the circle leaves out the lumen's corners. It is not a substate:
Schmitz et al. 2022's independent active-state deposit 7T3T reads the same
as 8TKF (85 against 106 pS). Nor is it the exit window, which moves the
answer by 2 %. At bulk diffusivity and a 1 Å K+ exclusion, RyR1's open
deposit gives 787 pS against 801 measured, and ITPR3 gives 261–278 pS
against 358–545.

The corners change how much current flows, but not where the voltage
falls (`python -m ip3r lumen`; Channel panel → "Draw the lumen"). In 8TKF,
7T3T and RyR1's 9HEO, the 3-D potential puts the half-drop point within
1.3 Å of the 1-D model's, and the filter's share within 4 points (8TKF: 32 %
against 35 %). The viewer draws the lumen coloured by that potential, and
the panel plots both curves.

![8TKF's lumen coloured by the potential, and where the voltage falls](docs/img/gui_lumen.png)

Selectivity is a ratio, so it does not depend on the unmeasured diffusivity
and tests the wall charge directly. Vais et al. 2010's own solutions were
run through the same pore (`python -m ip3r selectivity`). No reading of
8TKF's wall comes near the measured P_Ca:P_K of 15.2: 0.17 uncharged, 0.00
with its lining charges, 0.69 with only the acidic rings. The model's i_Ca
is at most 0.046 against 0.30 pA/mM. Under a local Donnan partition, the
lining lysines act as Ca²⁺ barriers. The calibration shows that charge in
discrete rings cannot make a continuum pore Ca²⁺-selective; a charged
tract can (`docs/SCIENCE_PERM.md`). Solving Poisson–Boltzmann across each
slice instead of assuming local Donnan (`selectivity --closure radial`)
screens the vestibule's K2529 ring, but lifts the charged reading only to
0.04 (0.01–0.16 over ε 80–10), so the closure is not the Ca²⁺ barrier.

Nor is protonation (`python -m ip3r protonation [8TKF|9HEO] [--corners]`).
Two independent pKa routes, a Tanford–Kirkwood network titrated by Monte
Carlo and PROPKA 3, both keep every lining group of 8TKF charged at Vais's
pH 7.3. The lysine rings stay charged even at a protein-like permittivity
of 4, where the acid rings lose up to half their charge. No combination of
rings formal or neutral exceeds P_Ca:P_K 0.69, and with the lysines charged
none exceeds 0.05. RyR1's open deposit, run through Xu et al. 2006's own
protocol as the control, fails the same way: 0.46 against 7.0 under every
pKa reading. The model also misses the order of Xu's charge mutants. The
missing piece is therefore the continuum model, not the IP3R wall. Xu's
E4900N also shows that PROPKA wrongly buries E4900 (pKa 8).

![Unitary conductance](docs/img/gui_unitary.png)

**Ryanodine receptors.** Rabbit RyR1 loads beside the IP3Rs: six deposits
chosen from the PDB by stated rules (`scripts/curate_ryr.py`: full-length,
wild type, activators only, ≤ 4 Å, best per stated state), one each closed,
primed, open, inactivated and closed-inactivated, plus a primed deposit from
the open state's own paper for the morph (9R8O → 9HEO). They go through the
same measurement: every one is in P11716 numbering, the shut states gate at
I4937, and only open 9HEO widens (5.05 Å). Its modelled K+ conductance is
136 pS neutral and 180 pS charged, against 801 pS measured (Xu et al. 2006),
so the continuum model is short on both receptors. RyR1 also allows a
stronger test than one number. Five charge-neutralising mutants were
measured, and the model makes each one on the structure
(`python -m ip3r mutants`). It gets the direction right for all four lining
residues and no effect for E4955Q, which does not line the pore. It misses
the largest effect: D4899Q cuts the conductance to 0.20× but the model
predicts 0.90×. And D4899, like ITPR3's homologous D2478, is salt-bridged,
so cancelling ion pairs predicts no effect at all. The bridge does not
neutralise that charge.

![RyR1 charge mutants](docs/img/gui_ryr_mutants.png)

**RyR1 gating and sparks.** Dynamics → Gating → "RyR1" draws Stern et al.
1997's two-gate scheme against Murayama et al. 2015's measured rabbit-RyR1
bell. The scheme activates where RyR1 does (3.9 vs 4.4 µM) but inactivates
6.7× too readily (48 vs 320 µM). Puffs → "RyR1 sparks" runs a 30-channel
cluster with the Ca²⁺ coupling derived from one channel's current at 30 nm.
Uncoupled it gives only blips. Coupled it gives 1.5 sparks per second that
recruit nearly the whole cluster, read with the same ruler as IP3R puffs.
With a single cluster Ca²⁺ the sparks last ~120 ms, against ~6 ms
measured. "RyR1 sparks in the cleft" puts the same channels in Stern's
junctional cleft (two rows, 60 × 15 nm, edges leaking), where each channel
sees its own Ca²⁺ from a steady diffusion solve. Sparks there end by local
inactivation after ~20 ms: six times shorter, still 3× the measurement, and
no geometric uncertainty closes the gap. The Gating tab also draws the
scheme fitted to both measured flanks (dashed: Ka 4.9, Ki 249 µM), and
"RyR1 in the cleft, gating fitted to Murayama 2015" runs it. Once a spark
starts, it never ends. With inactivation made as weak as the measured bell
says, the tens of µM in the cleft cannot shut the array at any
inactivation rate. So what ends a real spark is missing from a Ca²⁺-only
scheme.
Adding the fibre's 1 mM free Mg²⁺ (`spark-mg`; in the GUI, Puffs → Mg²⁺
and "Triggered sparks vs Mg²⁺", and the dash-dot bell in Gating), from open
sources, supplies it. Competing at the activation site, Mg²⁺ shuts
triggered sparks with no channel inactivated: the array's feedback falls
below one. How fast depends on the site's Mg²⁺ affinity. Meissner et al.
1997 measured it in the assay behind the fitted bell. Carried into that
assay's 0.17 M NaCl, where Na⁺ already holds the same site, it gives
K_Mg,A ≈ 770 µM (range 345–1,600 µM), and the activation site alone ends a
triggered spark in 32 ms at best (435 ms at the central value), against
6.3 ms measured. With Mg²⁺ also at the inhibitory site, which holds ~80 %
of channels shut before any trigger, sparks end in 6 ms under every
reading.
The V channels, the half of the couplon opposite the voltage sensors, are
now simulated (`ec`; Rios 1993's allosteric model, Stern's rates). With
Stern's constants the couplon reproduces his release under voltage clamp:
a peak, a plateau, and a stop when the membrane repolarises. With the
C channels fitted to the measured bell, no Mg²⁺ arrangement gives a peak.
Either release goes on after repolarisation, or the C channels hardly
open. Emptying the SR (`ec --depletion`) restores a peak only by releasing
more of the store than a fibre loses. Nor does a two-site inactivation
gate fitted to the bell's slope (`ec --two-site`) help. It inactivates
less at cleft Ca²⁺ than the one-site gate, and control is lost further.

![RyR1 sparks](docs/img/gui_sparks.png)
![RyR1 sparks in the cleft](docs/img/gui_sparks_cleft.png)
![RyR1 in the cleft, fitted to the measured bell](docs/img/gui_sparks_fitted.png)
![Triggered sparks against free Mg²⁺](docs/img/gui_sparks_mg.png)

**Check the publication.** The Findings tab re-derives 45 results from the
six `ip3r_genes` papers and its structural baseline, by three routes of
different strength: `recomputed` from coordinates with code the two projects
do not share, `rederived` from the publication's input tables with this
project's arithmetic, and `read` (the table read, the prose tested).

![The findings panel](docs/img/gui_findings.png)

Residue-keyed findings can be drawn on any structure in human numbering.
Paper 6's two modules, for example, are shown as Cα traces: the ligand core
in green and the pore module less the luminal loop in magenta.

![Paper 6's modules on 6DQN](docs/img/gui_modules.png)

**Colour by distance to IP3.** "Distance to IP3 (S22 shells)" paints every
residue by its all-atom distance to the IP3 bound on its own subunit, in
S22's four shells (contact < 4.5 Å, then 8, 11.5 and 15 Å). Residues beyond
15 Å, and subunits with no IP3, are grey. The shell checks' exhibit plots
conservation against that distance for all three paralogs.

![Ligand shells and conservation against distance](docs/img/gui_shells.png)

**Read the tree.** The Tree tab draws Paper 2's RyR-rooted maximum-likelihood
tree (134 proteins) from the committed `rooted.nwk`, parsed by this
project's reader. Each paralog's whole clade is boxed with its size and
SH-aLRT/UFBoot support, and so is the RyR outgroup. Hagfish and lamprey tips
are yellow diamonds, labelled with the support of their clade and of the
node where it joins. A white dot marks each node that clears both support
bars (80/95). "Vertebrates" zooms to the 57 vertebrate tips with their
labels. "Show" on any tree check opens this tab. "Beside --bnni" draws the
model-violation re-search next to it, both rooted on RyR, and lists what the
guard did to each clade claim: 9 of 10 held, and the ITPR1 core (never well
supported) fell from 47.8/95 to 47.5/73.

![Paper 2's tree](docs/img/gui_tree.png)

![The reported tree beside the --bnni re-search](docs/img/gui_tree_bnni.png)

**See every genome.** The Genomes tab draws Papers 3 and 4 as a grid: a row
for each of the 309 assemblies in the retention sweep and a column for each
cell (ITPR1–3 and the RyR control). Each row sits beside a contig-N50 strip
on a fixed log scale. You can colour the cells by what the sweep found, by
the known genes it missed, by the S15a evidence state, by how a protein
search could reach the gene, or by lesions: whether the cell carries more
frameshifts and stops than its own genome's identity-matched siblings
(S15b §8, rebuilt from the per-locus table). Rows sort by N50 (the contiguity bar is drawn),
by class, or by name, and can be filtered to a class or to assemblies above
the bar. Clicking a row names the genome and lists its four cells. "Show" on
a P3 or P4 check opens this tab on the right layer. `P3.lesion_strata`
opens it on the birds, where ITPR3's excess (25 genomes to 2) sits mostly
below the bar.

![Misses against contiguity](docs/img/gui_genomes.png)

![ITPR3's lesion excess in the birds](docs/img/gui_genomes_lesion.png)

**See where the receptor is.** The Range tab draws Paper 1 as one bar per
clade of the 6,928-proteome sweep. Each bar shows the fraction of swept
proteomes carrying an IP3 receptor call, on a fixed 0–1 scale, coloured by
supergroup; archaea and bacteria are collapsed to one row each. A red cross
marks a clade whose absence held in controlled genome assemblies (a small
one where it held for a class inside the clade). You can hide small clades
or expand the prokaryotes. Clicking a row lists its genome-level absences.
Double-clicking a clade (or View → "Genomes (S23)") draws its S23 genome
assemblies one by one: each genome's control verdict, what the copy ledger
found, whether its contig N50 reaches its own contiguity bar, and its
complete gene models on a fixed 0–20 scale. "Show" on any P1 check opens
this tab.

![Paper 1's range](docs/img/gui_range.png)

![One clade's S23 genomes](docs/img/gui_range_genomes.png)

**See the variants, and where the uncertain ones sit.** The Variants tab lists
the S17 harvest for one paralog and class. "Draw on structure" puts a sphere
on every variant residue of that class on all visible subunits: P/LP red,
B/LB blue, conflicting violet, VUS amber. It draws only on a deposit in that
paralog's human numbering; on any other it says why and draws nothing. Choose
a layer under "VUS by layer" and each VUS is placed against its own gene's
labelled medians on that layer, as in Paper 5 §8. A VUS at or above the P/LP
median is *pathogenic-like* (pale red), one at or below the B/LB median is
*benign-like* (pale blue), and an unscored one is grey. The table gains a
stratum column and the plot shows the three classes with both medians. This
is a stratification, not a call.

Each median is itself uncertain, so the plot shades its exact 95 % interval
(a distribution-free interval from order statistics). A VUS that a
threshold inside an interval could move is drawn hollow and marked "near
the P/LP (or B/LB) median" in the table. A class with five positions or
fewer cannot bound its median: ITPR2's P/LP median is one position and
ITPR3's is five. The panel says so, and no VUS of theirs counts as firmly
pathogenic-like.

![Variants and the VUS stratification](docs/img/gui_variants.png)

![ITPR2's thresholds, one of them unbounded](docs/img/gui_variants_bands.png)

**See what the map leaves out, as a prediction.** Representation →
Completeness adds the AlphaFold model's residues where a deposit has none.
"+ AlphaFold gaps" fills the internal stretches, each fitted on resolved
residues on both sides. "+ gaps and ends" also adds the termini, which have
one side to fit on. The fill is drawn in AlphaFold's own pLDDT colours
whatever the deposit is coloured by. Each seam is a bond from the deposit to
the fill: pale where it closes, red where it does not. The fill follows a
morph or mode frame. It is only drawn: every measurement in the application
still runs on the deposit. The model is chosen by measuring which downloaded
prediction is in the deposit's numbering, or failing that, which one is the
same protein through an alignment of the deposit's construct. So rat 7LHF is
filled from rat ITPR1 isoform 8. The two splice segments the isoform lacks
are left empty, and the panel says so. 9YKK (ITPR2) has no model and is
refused, and the panel shows why.
On 8TKG, 1,592 residues are filled over the four subunits, at a mean pLDDT
of 38. AlphaFold is least sure exactly where the map is empty. Filling 16
stretches that 8TKG does resolve (but another deposit does not) lands at a
median 1.1 Å, against 5.5 Å for a straight line (`python -m ip3r graft 8TKG --calibrate`).
Length does not break a fill, but low confidence does
(`python -m ip3r graft 8TKH --long`). Hidden 60-residue windows fill to
1.5 Å where AlphaFold is confident. The one resolved stretch below pLDDT 50,
8TKH 926–943, was placed 58 Å off, further than a straight line. So the
summary counts the residues below pLDDT 50 (1,308 of 8TKG's 1,592) and says
they are not positions.

![AlphaFold fills with their seams](docs/img/gui_alphafold.png)

## The checks, as of Round 7.7

52 checks: 50 confirmed, 2 discrepancies. Both discrepancies are genuine, and neither
touches a published paper's headline:

- **`P6.contacts_heavy_atom`.** S22's positive control — S0's ten IP3 contacts
  recovered in all six IP3-bound depositions — is reproduced exactly with
  S22's rule, which counts hydrogens. Under S0's own heavy-atom definition,
  Arg503 is 4.78 Å (8TKG) and 4.83 Å (8TKH) from IP3, so the control holds
  in 4 of 6 depositions, and in the other two only through a hydrogen.
- **`P5.report_both_metrics`.** The S17 report says the gate and filter are
  the most constrained elements "on both metrics and in all three
  paralogs". On the JSD the ITPR1 top two are RIH-associated and gate, and
  the ITPR2 top two are gate and the β-trefoil. The paper's own Results
  state the narrower, correct version, which `P5.gate_filter_most_conserved`
  confirms.

The table of every check, its kind and what it re-derived is in
`docs/SCIENCE_CHECKS.md`.

Every check is calibrated: `tests/test_checks_calibration.py` runs each one on
a copy of only the tables it declares, then plants a change and requires the
verdict to flip.

## Installing

```
bash scripts/create_env.sh          # or: conda create -n ip3r_sim --clone piezo1
conda activate ip3r_sim
make fetch                          # the registry's 16 mmCIF files, ~46 MB, into ref/
```

Clone `ip3r_genes` beside this repository (or set `IP3R_GENES_DIR`) for the
findings checks; everything else runs without it.

## Running

```
./run_app.command                   # the GUI (activates ip3r_sim; double-click in Finder)
python -m ip3r                      # the GUI, from an activated environment
python -m ip3r --session view.json  # the GUI, reopened on a saved session
python -m ip3r checks [--paper constraint] [--figures out/]
python -m ip3r states               # the ITPR3 gating states at the pore
python -m ip3r unitary              # their K+ conductance, vs 358/545 pS
python -m ip3r shortfall [--scan]   # every open deposit in 1-D and 3-D vs the measurement
python -m ip3r lumen [8TKF ...]     # where the voltage falls: 3-D potential vs the 1-D model
python -m ip3r selectivity          # 8TKF's P_Cl:P_K, P_Ca:P_K, i_Ca vs Vais 2010
python -m ip3r protonation [9HEO] [--corners]  # lining pKas (network, PROPKA) and selectivity under each
python -m ip3r states --paralog RYR1   # the curated RyR1 states (also unitary)
python -m ip3r mutants              # RyR1 charge mutants: model vs Xu 2006
python -m ip3r transition 9R8O 9HEO # RyR1 primed -> open
python -m ip3r ryr-gating | sparks --scan   # RyR1 bells; sparks over the coupling band
python -m ip3r sparks --cleft      # sparks with each channel's own Ca2+ in the cleft
python -m ip3r sparks --cleft --fit  # the same with Ka, Ki fitted to Murayama's bell
python -m ip3r spark-termination --scan fit|ki|rate|use|ratio|fraction|low-activity|speed  # what ends a cleft spark
python -m ip3r spark-mg [--scan] [--reading measured|selectivity|meissner]  # Mg2+ ends a triggered spark
python -m ip3r ec [--scan] [--reading R] [--use]  # the couplon under voltage clamp (V + C channels)
python -m ip3r ec --two-site       # ... with a two-site inactivation gate fitted to the bell's slope
python -m ip3r ec --depletion [--pool-scan]  # ... with the SR emptying (Stern's Fig. 20 pool)
python -m ip3r info 8TKF            # one deposit, measured
python -m ip3r modes 6DQN           # normal modes with C4 irreps
python -m ip3r transition 8TKG 8TKF # morph, displacement, mode overlap (--gate: pore per frame;
                                    #  --cutoff-scan, --stride-check: the network against itself)
python -m ip3r graft 8TKG --calibrate   # AlphaFold fills, seams, and how good they are
python -m ip3r graft 8TKH --long        # ... on long windows and islands; 7LHF by alignment
python -m ip3r gating | oscillate --window | puffs --ip3 0.2
python -m ip3r puffs --model park-drive | puffs --scan   # the two receptors
make help
```

Mouse: left-drag rotate, shift-drag pan, wheel zoom, click to identify a
residue (element and conservation shown). `Ctrl+1` side view, `Ctrl+2` down
the pore, `Ctrl+3` one IP3 site (its 15 Å pocket, the tetramer in front of
it clipped away), `Ctrl+0` fit to view, `Space` spin. Until you move the
camera, the view keeps the molecule filling the viewport as the window or
the visible subunits change. "Show on structure" for the IP3-contact and
ligand-shell checks opens the site view. The deposition list is grouped by
family, with RyR1 collapsed.

**Changing a parameter.** Help → Parameters… (`Ctrl+Shift+P`) lists every
registered number with its default, bounds and source; double-click a value
to change it (out-of-range values are clamped, and the clamp is reported).
While anything differs from its default an amber banner runs across the top
of the window, and the findings checks report *not run* rather than confirm.
Nothing typed is remembered after you quit. "Export…" writes the set as
JSON, which `IP3R_PARAMETERS=file.json python -m ip3r …` reproduces headless.
A control whose starting value is a parameter (the Puffs cluster size and
coupling, the microdomain run length) follows an edit unless you have typed
your own value there. Drawings nobody asked for (the Gating plot, the
displacement and ligand-shell colours and their legend) are redrawn.
Results you ran (a puff simulation, a channel measurement, modes, a
transition) keep the values they were computed with until you run them again.

![Parameter editor](docs/img/gui_parameters.png)

**Saving where you were.** File → Save session… (`Ctrl+Shift+S`) writes the
view as JSON: the deposit, style, colouring, layer, subunits, marked sites,
pore, camera, open tab, any transition built (end, fit, method, frame), the
normal mode animating and its amplitude, the Dynamics settings (gating model,
oscillation, puffs, microdomain) and the Variants view (paralog, class, layer,
drawn). A mode animation is recomputed and restarted; a simulation is not
re-run, only its settings come back. File → Open session… (`Ctrl+O`) or `--session` puts it back. A session holds
no coordinates and no results: it is re-derived from the same inputs on
opening. It does record the parameter overrides in force when it was saved.
If they differ from the current set, you are asked whether to apply them
(the deposit is then re-measured) or keep your own.

## How it is built, and how it checks itself

See [`INTERFACE.md`](INTERFACE.md) for the module map and
[`docs/SCIENCE.md`](docs/SCIENCE.md) for the models (the ryanodine
receptor in [`docs/SCIENCE_RYR.md`](docs/SCIENCE_RYR.md), and the couplon in
[`docs/SCIENCE_EC.md`](docs/SCIENCE_EC.md)). Every number a
calculation uses is a registered parameter with a unit, bounds and a source
(`python -m ip3r params`); curated data imported from `ip3r_genes` records
the SHA-256 of its source tables (`make sync-check`). Files stay under 500
lines. What is next is in [`ROADMAP.md`](ROADMAP.md).

## References

De Young & Keizer 1992 (PNAS 89:9895); Li & Rinzel 1994 (J Theor Biol
166:461); Bezprozvanny, Watras & Ehrlich 1991 (Nature 351:751); Mak, McBride
& Foskett 1998 (PNAS 95:15821); Swillens et al. 1999 (PNAS 96:13750); Shuai &
Jung 2002 (Biophys J 83:87); Siekmann et al. 2012 (Biophys J 103:658); Cao et
al. 2013 (Biophys J 105:1133); Cao et al. 2014 (PLoS Comput Biol
10:e1003783); Smith & Parker 2009 (PNAS 106:6404); Paknejad &
Hite 2018 (NSMB 25:660); Atilgan et al. 2001 (Biophys J 80:505); Yang, Song &
Jernigan 2009 (PNAS 106:12347); Tama & Sanejouand 2001 (Protein Eng 14:1);
Brüschweiler 1995 (J Chem Phys 102:3396); Kabsch 1976 (Acta Cryst A32:922); Hanley &
McNeil 1982 (Radiology 143:29); Barlow & Thornton 1983 (J Mol Biol
168:867); Xu et al. 2006 (Biophys J 90:443); Stern, Pizarro & Ríos 1997
(J Gen Physiol 110:415); Murayama et al. 2015 (PLoS One 10:e0130606); Ríos
et al. 1999 (J Gen Physiol 114:31). Full entries: `ip3r/resources/references.json`.
