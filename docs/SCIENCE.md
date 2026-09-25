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

**Which model, and which residue.** AlphaFold DB (queried 2026-09-24) holds
canonical ITPR3 (`AF-Q14573-F1`), ITPR1 isoform 4 only, rat ITPR1 isoform 8
only, and a 181-residue ITPR2 isoform. A `NumberMap`
(`structure/graft_numbering.py`) says which prediction residue fills which
deposit residue. There are two routes, tried in order:

- **By number**, against the variant-painting bar (`numbering.min_identity`).
  The seven ITPR3 deposits match `AF-Q14573-F1` at 99.1–100 %.
- **By alignment** (Round 7.9). The deposit's whole construct, including the
  residues no atom was built for, is read from the mmCIF's
  `_pdbx_poly_seq_scheme` (`io/poly_seq.py`). It is aligned to each model
  with `core.pairwise`, and the aligned pairs are the map. The model must
  agree with the construct over the pairs to `graft.align_min_identity`
  (0.99). Rat 7LHF (numbered in canonical P29994 6–2741) against rat isoform
  8 gives 1.000 over 2,681 pairs, with two unpaired segments: 318–332 (SI)
  and 1693–1732 (SII), the splice segments the isoform lacks. Human ITPR1
  isoform 4 gives 0.988 (the ortholog), ITPR3 0.658. 9YKK (ITPR2) reaches
  0.713 at best and is still refused.

A stretch that touches an unpaired segment is skipped and named ("the model
lacks residues 318–332"), and so is one where the model has extra residues.
The fill's atoms carry deposit numbers on either route. On 7LHF (gaps only):
40 stretches and 1,264 residues, mean pLDDT 32. Three stretches per subunit
are not filled: two touch the splice segments and one has no anchors
before it. Every filled residue's name equals the construct's.

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

**Long stretches** (Round 7.9; `python -m ip3r graft PDB --long`). The
stretches above have pLDDT 53–77 and are 4–10 residues long. The real gaps
run to 70 residues, with a median pLDDT of about 40. On 8TKG, 1,308 of the
1,592 filled residues are below 50. Two tests reach further:

- **Windows.** Resolved stretches of 10–60 residues, clean anchors either
  side, are hidden and filled. On 8TKG the median fill is 0.49 Å at 10
  residues and 1.52 Å at 60, against 4.3 and 16.1 Å for the line. 8TKH and
  7LHF (through the alignment) are alike. **Length does not break a
  fill.** But these residues are ordered (pLDDT ~82), so they are the easy
  case. The seams of right fills (< 2 Å) stay ≤ 5.22 Å in 802 of them, so
  `graft.join_tolerance` holds.
- **Islands.** No deposit resolves a run below pLDDT 50 with its own anchors
  either side. What the maps do resolve are islands: short resolved runs
  (5–60 residues) with a gap on each side. Hiding one and filling the whole
  span from the outer anchors is a real 27–137-residue fill, scored where
  the experiment saw residues. Over 28 islands in 9 deposits, the 26 whose
  residues have pLDDT ≥ 70 land at a median 1.26 Å, against 17.8 Å for the
  line, and all beat it. 7LHF 1025–1045 (pLDDT 63) is 11.0 Å against 13.1.
  **8TKH 926–943 is the only island below 50 (pLDDT 31)**, in a 62-residue
  span. Its fill lies 57.9 Å from the deposit, and the line 24.4 Å.

So pLDDT, not length, decides whether a fill is a position. The one test in
the real gaps' regime failed, and by more than a line would. That is one
stretch, not a rate, but nothing supports the opposite. The fill summary
counts the residues below `display.plddt_low` and warns that they are not
positions. They are still drawn, in AlphaFold's very-low colour: they show
how much chain the map leaves out, not where it is.

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

All three IP3R models (with park/drive's stationary bell, Round 7.3) are
measured with one ruler (`physics.bell`): the Ca²⁺ at half the bell's own
peak on each flank. From 33 nM (the lowest IP3 at which the paper says
activation was unaffected) to 10 µM:

| model | half-activation | half-inhibition |
|---|---|---|
| De Young–Keizer | 2.01× | 2.76× |
| Mak 1998 | 1.016× | 6.22× |
| Park/drive (stationary, gates at equilibrium) | 1.09× | 42.8× |

Park/drive passes Mak's test almost as well as Mak's own fit, although
Siekmann's rates were fitted to stationary records and not to this
comparison. Its half-activation sits at 0.28–0.39 µM from 10 nM to 10 µM
IP3. The inhibitory flank moves 40× because IP3 lifts the m42·h42 switch
out of park mode at high Ca²⁺ (1.4 µM at 0.1 µM IP3, 56 µM at 10 µM). Its
bell never reaches zero: park mode's O5 leaves a floor of about 10 % of
the peak at both ends. The Gating panel's "three models side by side" draws
the bells at 33 nM and at 10 µM, each relative to its own peak.

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

## The findings checks, paper by paper

How each check re-derives its paper's numbers (the module contrast, the
ligand shells and their substitution rates, the tree and its `--bnni`
re-search, the genome grid, presence and absence, the family-call
benchmark, the VUS strata) and what agreement establishes:
`docs/SCIENCE_CHECKS.md`.
