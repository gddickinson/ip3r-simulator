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
- **Salt bridges cancelled** (`physics/salt_bridges.py`, the *paired*
  reading). An ion pair is an acid and a base with a charged O and a charged
  N within 4 Å (Barlow & Thornton 1983, `pore_charge.salt_bridge_cutoff`).
  Pairs are matched one-to-one, closest first, across the whole deposit, so
  each bridge removes exactly +1 and −1. A lining group that is half of a
  pair is dropped. On 8TKF only the four D2478–R2471′ bridges (2.45–2.58 Å)
  qualify: the wall goes from −8 e to −4 e and the conductance from 33 to
  **23 pS** (8–54 pS over the diffusivity × radius sweep). The next pairs
  are close to the line. D2518–R2524′ sits at 4.23–4.34 Å and K2482–D2400
  at 4.31–4.40 Å, so a 4.5 Å cutoff drops those too and gives 38 pS. Over
  cutoffs of 3–6 Å the paired reading is 23–38 pS, always below the neutral
  65 pS. The peak partition density (18.6 M, above the ceiling) does not
  move, because it sits on the D2518/D2522 rings at z ≈ −63 Å in a 4.4 Å
  lumen, not at the filter. So cancelling ion pairs does not rescue the
  charged model. The gap to 358 pS is not a charge-counting artefact.
- **The paired reading is refuted by RyR1** (next section). RyR1 D4899 is
  the homologue of D2478 (own alignment of the PF00520 domains: R4892 ≡
  R2471 as well). It is bridged the same way (2.9 Å), yet neutralising it
  cuts RyR1's conductance to 0.20×. So a bridged filter carboxylate still
  acts as a charge. The *charged* reading is the better one. *Paired* is
  kept as a reported bound, not a correction.

What survives: the gate is the only state change that opens a conducting
pathway. Taken as a neutral continuum, the only open deposit's pore is
too narrow or too long to carry the measured conductance.

## Ryanodine receptor 1

**Why.** RyR1 is IP3R's closest relative (both share the trefoil, MIR, RIH,
RIH-associated and PF00520 domains), with many more deposits, a larger
conductance and measured pore mutants. It is the natural control for
anything measured on the ITPR3 pore.

**Resource** (`scripts/curate_ryr.py` → `resources/ryr1.json`, SHA-256 of
every response). Reference rabbit P11716 (5,037 residues), because the
deposits are rabbit; the GGGIGD filter is at 4894–4899 and the gate
isoleucine at 4937, as in the literature. Domains are InterPro's Pfam
matches; PF00520 (4789–4947) is named `channel`, so the pore-domain span is
found the same way as for the ITPRs. RyR1 has no conservation, sites or
variants here, so those colourings are grey on it.

**Panel, by rule**, from the 158 PDB entries mapped to P11716 (on 2026-09-24):
single-particle EM of four full-length chains (not domain crystals or local
refinements); wild type; only Ca²⁺/Zn²⁺/Mg²⁺, ATP/ACP/ADP, caffeine or
lipid as ligands; ≤ 4.0 Å; the state the title names, best resolution per
state. That gives 9OL6 closed (3.11 Å), 8RRX primed (3.10), 9HEO open
(3.40), 7TDG inactivated (3.80) and 7TDI closed-inactivated (3.30). The morph
needs one preparation, so the best primed deposit of the open state's own
paper, 9R8O (3.30, same DOI as 9HEO), is added. 140 entries fail a rule,
each recorded in the resource with the first rule it fails; 18 pass, and
the best per state is taken from those.

**Measured.** Every deposit is 100 % in P11716 numbering. The four shut
states gate at I4937 (2.6–3.3 Å r_min) and 7TDI at Q4933 (2.21 Å). Only
9HEO opens (5.05 Å, lined by Q4933). Primed 9R8O → open 9HEO morphs over
4,235 residues × 4 subunits (RMSD 2.72 Å after the pore fit). The start's
elastic network explains little of it: the cumulative overlap of 20 modes is
0.174 against a null of 0.031, where ITPR3's is 0.64. Most of the raw
displacement (71.5 %) is rigid-body, removed before the overlap.

**Conductance.** Symmetric 250 mM KCl, the bath of Xu et al. 2006's
recombinant RyR1 (801 ± 7 pS, n = 17, planar bilayer). 9HEO: 134 pS series,
136 pS neutral, 180 pS charged, 161 pS paired; the sweep reaches at most
308 (neutral) and 449 pS (charged). The neutral model is 5.9× short. So
ITPR3's 2.4× shortfall is not particular to ITPR3. The same continuum
under-predicts both receptors, and RyR1's acidic wall (−32 e) *raises* its
conductance where ITPR3's alternating rings lowered it. Primed 9R8O
(r_free 1.62 Å) conducts 4.5 pS neutral. Its charged solve does not
converge (Cl⁻ is excluded, and K⁺ alone cannot neutralise a basic zone), and
it is reported as n.c., never as 0.

**Charge mutants** (`physics/ryr_mutants.py`, `python -m ip3r mutants`). Each
mutant drops one residue's charge on all four subunits of 9HEO. The ratio
mutant / wild type cancels the transport constants, so it tests the wall
charge alone:

| mutant | measured | model, charged | model, paired | in 9HEO |
|---|---|---|---|---|
| D4899Q | 0.20 | 0.90 | 1.00 | lining, bridged to R4892′ |
| E4900N | 0.63 | 0.76 | 0.41 | lining |
| D4938N | 0.65 | 0.78 | 0.80 | lining |
| D4945N | 0.92 | 0.90 | 0.91 | lining |
| E4955Q | 1.01 | 1.00 | 1.00 | not lining |

The direction is right for all four lining residues, the size right for
D4945N, and the null right for E4955Q. The filter aspartate is 4× too weak
in the charged model and absent in the paired one. A continuum Donnan
partition spread over 3 Å does not carry what a ring of four carboxylates
does at a 5 Å filter.

## RyR1 gating and sparks

**Sources, and why these.** A search for RyR1 models whose constants can be
read from the source's own text (many papers are behind a script
challenge) found one measured bell and one kinetic scheme.

- *Murayama et al. 2015* (PLoS One, S1 Table): recombinant rabbit RyR1,
  Ca²⁺-dependent [³H]ryanodine binding, A = Amax fA (1 − fI) with
  KA 5.5 µM, nA 1.2, KI 0.27 mM, nI 1.5 (wild type, 25 °C). Binding is an
  activity index, not P_open, so it is compared by shape (flanks), never by
  height.
- *Stern, Pizarro & Ríos 1997* (J Gen Physiol, Table I): the skeletal "C
  channel", two gates in series. Activation opens on two Ca²⁺
  (k_o = 10 µM⁻² s⁻¹, k_o− = 500 s⁻¹) and inactivation closes on one
  (k_i, k_i− = 20 s⁻¹). k_i is printed "2 × 10⁻⁶ M⁻¹ s⁻¹", a sign typo: the
  text gives the inactivation Kd as 10 µM, so k_i = 2 µM⁻¹ s⁻¹. The printed
  value would remove inhibition entirely (tested). The authors say the
  model "has not been objectively 'fitted' to data".

The best-constrained readable kinetic scheme, Zahradníková et al. 1999, is
cardiac RyR2 and has no Ca²⁺ inhibition, so it was not used.

**Bells on one ruler** (`python -m ip3r ryr-gating`). Half-activation: 3.9 µM
(scheme) vs 4.4 µM (measured). Half-inhibition: 48 vs 320 µM. The scheme
activates where RyR1 does, but inactivates 6.7× too readily. Stern et al.
said as much of their inactivation site ("one or two orders of magnitude
lower than in bilayers").

**Sparks** (`physics/sparks.py`, `python -m ip3r sparks [--scan]`). The
cluster is the 30 Ca²⁺-gated C channels of Stern's 60-channel couplon, with
the mean-field coupling the IP3R puffs use. The coupling is derived:
free-diffusion Ca²⁺ one channel spacing (30 nm) from one open channel
(0.3 pA, D 5 × 10⁻⁶ cm² s⁻¹), which is 8.25 µM, an unbuffered upper bound.
Each step is exact (expm of the generator, tabulated per number open). The
step was measured: 1e-5 and 2.5e-5 s agree within noise, 1e-4 s lengthens
sparks by 15–20 %, and a 5 ms step is caught.

Read with the puff ruler (10 s):

- Uncoupled: only single-channel blips.
- Coupled: 1.5 sparks per second, most reaching 25–30 of 30 channels,
  with a size gap between blips and sparks; Fano 4.
- Coupling scan (× derived value): sparks switch on between 0.05× and 0.1×
  (0.44–0.80 µM per open channel) and saturate by ~0.2×. So buffering
  could cut the coupling five-fold without changing the answer.

**Where it fails.** Sparks last ~120 ms (median), against a measured
release of 6.3 ms (Ríos et al. 1999, frog). The trace shows why. After the
first near-whole-cluster peak, the mean-field cluster settles at the point
where 30 · P_open(0.1 + 8.25 n µM) = n (between 5 and 6 open). It stays
there until the number open falls to zero by chance. Stern et al.'s
termination relied on local geometry and voltage-sensor coupling, which a
single cluster Ca²⁺ cannot represent. Measured spark durations need a
spatial Ca²⁺ field, as they had, not a better constant.

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
