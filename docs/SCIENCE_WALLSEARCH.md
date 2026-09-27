# What a wall can do to Cl⁻ and Ca²⁺ at once (Round 7.25)

`python -m ip3r wallsearch [PDB ...] [--only FAMILY] [--no-ceiling]
[--no-rings] [--no-reversal]`; `ip3r/physics/selectivity_bound.py`,
`ip3r/physics/wall_search.py`; `tests/test_selectivity_bound.py`,
`tests/test_wall_search.py`. Continues `SCIENCE_REVERSAL.md` (Round 7.23).

## The question

Round 7.23 read Vais 2010's two ratios at their own reversal protocol
through the 3-D charged lumen. The deposits' walls shut Cl⁻ out (8TKF
P_Cl:P_K 0.003, measured 0.27) and still select Ca²⁺ only about 1:1
(measured 15.2). Scaling 8TKF's charge from 0 to 2× moved the pair along a
path that never came near the measurement. The emergent item: search wall
models scored on both ratios at once, for charge placements or missing
groups that pass Cl⁻ at 0.27 and still select Ca²⁺.

## The series bound

Take linear response along a single-file path: cross-section A(z), point
ions, energy zψ(z) in kT. A species of valence z then has permeability
P_z ∝ D_z / ∫ e^{zψ}/A dz. Divide by the same path uncharged and write
X = e^{ψ} and w for the normalised weight 1/A. Then
ρ_z = 1 / E_w[X^z]. The function log E[X^t] is convex in t (Hölder), and
t = 1 lies at ⅓(−1) + ⅔(2), so

    E[X]³ ≤ E[X^{−1}] · E[X²]²   ⇔   ρ_Cl ρ_Ca² ≤ ρ_K³.

In the measured ratios, with a = (P_Ca:P_K)/(P_Ca:P_K)₀ and
b = (P_Cl:P_K)/(P_Cl:P_K)₀ each over the uncharged pore's:

    B = b a² ≤ 1.

Equality holds only for a uniform potential, a tract of constant charge
density. Any other profile, however its charge is placed, spread, screened
or signed, pays for a gain in P_Ca:P_K with at least its square in
P_Cl:P_K. Vais's pair over the model's uncharged 8TKF at reversal (0.29,
0.26) is B = 0.93 × 58² ≈ 3,100.

Calibrated (`test_selectivity_bound`): a uniform potential reads e^{−zc}
by hand, on the bound; 2,000 random series profiles (2–60 segments, areas
1–200, ψ spread up to 4 kT) never exceed it.

The bound has three exits, and the round measures each:

1. **Parallel paths.** Two paths of opposite sign side by side have
   log-convex conductance: with ψ = ∓2 kT, B ≫ 1 (tested). In 3-D a lumen
   could in principle carry cations on one side and anions on the other.
2. **Non-linearity.** The measurement is a reversal, not linear
   response, and the wall's screening differs between Vais's Cl⁻ and Ca²⁺
   experiments. Round 7.23's scan already reached B 1.2 at 0.05× charge.
3. **Species-specific energies.** A Ca²⁺-only energy is binding that is
   not the mean potential. It leaves the Cl⁻ experiment, which holds no
   Ca²⁺, untouched, so b stays 1.

## The instruments

- `selectivity_bound.linear_ratios(pore, fixed, extra)`: Round 7.19's 3-D
  linear response (`pb`, point ions) for any charge map, plus any
  species-specific energy on the grid.
- `well_ceiling(pore, z_lo, z_hi)`: a Ca²⁺-only well of
  `wallsearch.ceiling_depth` (15 kT) over one band. Its resistance is then
  gone, so the reading is the most any well there can lift P_Ca:P_K in
  linear response. 1/ceiling is the share of Ca²⁺'s resistance left
  outside the band. Calibrated on a tube: the ceiling equals 1/(1 − the
  band's share of the neutral drop), where the band counts half a voxel to
  one voxel past its last planes because its edge faces conduct as the well
  (Scharfetter–Gummel weights). It holds within 1 % from 15 to 20 kT.
  Deeper wells (30 kT) stall the weighted Laplace solve, and those
  readings are refused.
- `wall_search.ring_search`: C4 rings of +q and −q at one height, 45°
  apart, at the wall's radius, every `wallsearch.ring_step` (5 Å) through
  the span, at q = 1 and `wallsearch.ring_charge_max` (2 e per site). This
  is the most direct test of the parallel-path exit that the tetramer's
  symmetry allows. Read in linear response, with B.
- `wall_search.search`: candidates at Vais's reversal (Round 7.23's steady
  PNP, point ions), scored S = |ln(P_Ca:P_K/15.2)| + |ln(P_Cl:P_K/0.27)|.
  The families are `charge` (the deposit's lining charge × a geometric
  grid from `wallsearch.scale_min` to 1), `well` (no charge; a Ca²⁺-only
  well of 2–8 kT, step `wallsearch.depth_step`, over the span, the filter,
  the gate or filter to gate) and `well + charge` (the best well with the
  charge scaled in). The Cl⁻ experiment depends only on the charge, so it
  is solved once per scale.

## Results

All numbers are at 1 Å (`reversal3d.spacing`). 8TKF is the open ITPR3
deposit and 7T3T the open-state control; Vais's pair is ITPR3's. Each
deposit's search takes about 45 minutes (both ran together).

**Exit 1, parallel paths: shut where it would matter.** Opposite C4 rings,
B at each height (1 and 2 e per site; linear response):

| z (Å) | 8TKF | 7T3T |
|---|---|---|
| luminal vestibule (−103/−101) | 1.01 / 1.00 | 0.87 / 0.80 |
| filter (−88/−86) | 0.81 / 0.62 | 1.13 / 1.28 |
| gate (−73 / −71) | 0.94 / 0.89 | 1.71 / 3.08 |
| cytosolic end (−58/−61) | 0.98 / 1.00 | 2.71 / 7.89 |

In 8TKF's lumen no ring pair moves B past 1.01. Oppositely charged
patches a few ångström apart share one mixed potential. 7T3T's gate
region is 20 Å wider, and there the pairs do reach B 7.9. That gain is
almost all in Cl⁻ (×3.4); Ca²⁺ rises only ×1.5.

**Exit 2, non-linearity: small.** 8TKF's charge scaled, at reversal:

| scale | 0 | 0.01 | 0.03 | 0.1 | 0.32 | 1 | *measured* |
|---|---|---|---|---|---|---|---|
| P_Cl:P_K | 0.291 | 0.246 | 0.173 | 0.066 | 0.013 | 0.003 | *0.27* |
| P_Ca:P_K | 0.26 | 0.29 | 0.37 | 0.58 | 0.88 | 0.98 | *15.2* |
| B | 1 | 1.08 | 1.20 | 1.16 | 0.52 | 0.13 | *3,239* |

(7T3T: B ≤ 1 throughout, P_Ca:P_K ≤ 1.41.) The reversal's non-linearity
buys at most 20 %.

**Exit 3, a Ca²⁺-only well: the only large lever, and not large enough.**
In linear response, even an unlimited well lifts P_Ca:P_K only so far
(`well_ceiling`, × the uncharged pore):

| band | 8TKF | 7T3T |
|---|---|---|
| the membrane span | 46 | 28 |
| filter ± 3 Å | 1.5 | 1.3 |
| gate ± 3 Å | 1.1 | 1.2 |
| filter to gate | 1.9 | 3.0 |

Vais's value needs 58× over the uncharged pore's reversal reading on
8TKF. Even a well covering the whole span leaves 2 % (8TKF) to 4 % (7T3T)
of Ca²⁺'s resistance outside it, in the access regions. A well at one
constriction covers too little of the path to matter.

At reversal, with its Ca²⁺ counted by Poisson (P_Cl:P_K stays the
uncharged pore's, 0.291 / 0.367):

| well | 8TKF P_Ca:P_K | 7T3T P_Ca:P_K | 8TKF peak Ca²⁺ |
|---|---|---|---|
| span, 2 kT | 1.73 | 3.01 | 0.07 M |
| span, 4 kT | **4.48** | **6.86** | 0.33 M |
| span, 6 kT | 4.39 | 5.90 | 1.0 M |
| span, 8 kT | 1.85 | 2.32 | 2.2 M |
| filter, 2 / 4 / 6 / 8 kT | 0.36 / 0.27 / 0.03 / < 0 | 0.63 / 0.52 / 0.26 / < 0 | |
| gate, 2 / 4 / 6 / 8 kT | 0.28 / 0.20 / 0.01 / < 0 | 0.60 / 0.59 / 0.50 / 0.34 | |

The ratio peaks near 4–5 kT and then falls. A deeper well fills with
Ca²⁺, and that uncompensated positive charge repels further Ca²⁺ entry.
A deep, uncompensated local well is a block: past 6 kT the GHK ruler reads
the Ca²⁺ reversal on the wrong side of the Cl⁻-only value (< 0). Adding the
deposit's own charge to the best well compensates it. 8TKF then reaches
5.37 and 7T3T 8.89, the highest readings found, but P_Cl:P_K falls to
0.003 / 0.001.

**The best wall.** On 8TKF it is the 4 kT span well, uncharged: P_Cl:P_K
0.29, P_Ca:P_K 4.48, score 1.30, against the uncharged pore's 4.15 and the
deposit's own wall's 7.40. On 7T3T it is the same well with 0.01× charge:
0.29, 6.71, score 0.90 (uncharged 3.69). Both are 2–3× short of 15.2.

## What it means

The emergent item asked which charge placement passes Cl⁻ at 0.27 and
still selects Ca²⁺. None does, and none can. For point ions in series, the
bound B ≤ 1 holds for every placement, sign, spread and screening of
charge. The two exits open to a mean-field wall are worth ≤ 1.2
(non-linearity) and ≤ 8 (parallel paths, only in 7T3T's wide gate, and
then in Cl⁻), against 3,239 needed. The measured pair is not a mean-field
electrostatic property of any wall with this lumen.

A Ca²⁺-specific energy leaves Cl⁻ alone, so it is the only lever that
moves the right ratio. Even that lever reaches only a third to a half of
15.2, and only when it spans the whole membrane at ~4 kT. It is
limited by the access resistance (linear ceiling 28–46×) and by the
charge of the Ca²⁺ it gathers. The measured P_Ca:P_K therefore needs an
interaction that this continuum does not carry: Ca²⁺ occupancy that
blocks K⁺ in the bi-ionic condition (the anomalous-mole-fraction
mechanism), or a Ca²⁺ affinity whose site is compensated as it fills.

Left open:
- a compensated site (fixed charge −2e per bound Ca²⁺, so filling costs
  no field), and a well that extends into the vestibules' access regions;
- the same search on RyR1 against Xu's single ratio, which has no Cl⁻
  constraint.
