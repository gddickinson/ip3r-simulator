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

## What the findings checks establish

| kind | what agreement means |
|---|---|
| `recomputed` | the published number follows from the deposited coordinates by a route sharing no code with the publication |
| `rederived` | the published output follows from the published inputs by independent arithmetic; says nothing about the inputs |
| `read` | the prose states what the table says |

A check is only trusted after it has been shown to flip on a planted change
(`tests/test_checks_calibration.py`).
