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
alone, which the DYK scheme does not reproduce (Round 4).

## Cell Ca²⁺ (closed-cell Li–Rinzel)

Release through the receptor, an ER leak and SERCA uptake, with total Ca²⁺
conserved. Measured here by simulation: sustained oscillations for IP3 in
0.36–0.63 µM (period ~11–13 s); damped spirals are not counted.

## Puffs

N receptors × 4 subunits × 3 two-state sites flipping with the DYK rates on
a fixed step (Shuai & Jung 2002), coupled through a mean-field cluster Ca²⁺
(Swillens et al. 1999). The measured signature of coupling is the Fano
factor of simultaneous openings (~1 independent; 1.43 coupled at 0.2 µM).
DYK's resting activity is high (n∞(0.1 µM) = 0.55), so puffs are modest.

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

## What the findings checks establish

| kind | what agreement means |
|---|---|
| `recomputed` | the published number follows from the deposited coordinates by a route sharing no code with the publication |
| `rederived` | the published output follows from the published inputs by independent arithmetic; says nothing about the inputs |
| `read` | the prose states what the table says |

A check is only trusted after it has been shown to flip on a planted change
(`tests/test_checks_calibration.py`).
