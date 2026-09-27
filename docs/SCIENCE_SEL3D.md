# Selectivity in 3-D (Round 7.19)

`python -m ip3r sel3d [PDB ...] [--mutants] [--spacing H] [--wall-volume V]`;
`ip3r/physics/csc3d.py`, `ip3r/physics/selectivity3d.py`;
`tests/test_csc3d.py`. Continues `SCIENCE_CSC.md` (Round 7.17) and
`SCIENCE_PERM.md` (the 3-D wall charge, Rounds 7.11–7.15).

## The question

Round 7.17 gave the 1-D pore charge–space competition (hard spheres + MSA).
RyR1's filter then bound Ca²⁺ as Gillespie 2008's does, yet P_Ca:P_K only
rose from 0.46 to 0.64 against Xu's 7.0. The reason: an uncharged gate
window is in series and holds 72 % of Ca²⁺'s resistance, and the 1-D
local-neutrality closure carries no field into an uncharged slice. The 3-D
wall field (Rounds 7.11–7.15) does reach past its charges. So:

1. With the wall's field reaching the gate in 3-D, and the csc excess in
   every voxel, does P_Ca:P_K approach the measured values (RyR1 7.0, ITPR3
   15.2)?
2. Does D4899Q then stand out among Xu's mutants?

## 1. The ruler: permeability from linear response

Between identical baths at small voltage, each species conducts as its
bulk σ_i = z_i² e² D_i c_i / kT times a Boltzmann-weighted Laplace solve on
its own volume (Round 7.11). GHK at 0 mV in the same bath gives
g_i = z_i² e² P_i c_i / kT, so **P_i = D_i × g_i/σ_i**, and

    P_Ca:P_K = D_Ca G_Ca / (D_K G_K),   G_i = the weighted Laplace conductance (m)

This is the independence regime's permeability. It is the integral Round
7.17 took in 1-D (`csc_readings.shares`, 1/P = ∫ e^E / (D A) dz), and
access resistance is included. **It is not a reversal potential.** In 1-D
the two differ (9HEO csc: 0.87 in linear response, 0.64 under Xu's
protocol), so every deposit prints its 1-D linear-response readings (Donnan,
csc) and its 1-D csc reversal reading beside the 3-D ones.

The bath is symmetric: the family's KCl (RyR1 250 mM, IP3R 140 mM) with 10
mM CaCl₂ (`selectivity.ryr1_cacl2_lumen`) on both sides.

## 2. The fluid on the voxels

- **The wall's groups** (`csc3d.wall_fluid`) are two half-charged oxygens
  per acid and one sphere per base. They are placed exactly as the fixed
  charge is placed: `slice` (the 1-D charge per length over each plane's
  lumen) or `local` (each group a 3-D Gaussian about its own centre). Their
  charge *is* the fixed-charge map (tested to 1e-12).
- **`local + csc` / `slice + csc`** (`csc3d.local_csc`): each voxel is in
  equilibrium with the bath under local neutrality, with Round 7.17's
  `csc.partition` applied voxel by voxel. Identical voxels (the slice
  placement fills a plane with one state) are solved once, and voxels no
  charge or group reaches are the bath.
- **`pb + csc`** (`csc3d.pb_csc`): Poisson's equation on the lumen as in
  Round 7.11's `pb` (no field into protein, u = 0 on the baths). Each
  species' excess over the bath, μ_i − μ_i,bath, is taken from the locally
  neutral fluid and held fixed in its Boltzmann factor
  (`charge3d.poisson_boltzmann(offset=)`).

**Why the excess comes from a neutral reference.** It was first tried
fully self-consistently: the fluid solved at the voxel's held potential,
with its stiffness in Newton's Jacobian. That solve finds no root below
u ≈ −5 kT/e. The MSA is a theory of a neutral mixture. With no neutrality
its screening term (−l_B Γ z², Γ growing as √ρ) outgrows ln c, and the
voxel's cations run away. Gillespie's density functional evaluates its
screening on a reference fluid for the same reason. Two consequences follow:
- Outside the groups' reach the reference is the bath, so the gate feels
  the mean field alone.
- Inside, the excess is the neutral fluid's, while PB's local composition
  departs from neutrality. This is the approximation the design accepts.

**Calibrated first** (`tests/test_csc3d.py`):
- the wall fluid carries the fixed map (both placements);
- with the terms off, `local + csc` = Donnan (1e-10) and `pb + csc` = PB
  (1e-9);
- in a long charged tube `pb + csc` = `local + csc` at the centre (1e-5,
  the Donnan limit);
- the fluid favours Ca²⁺ over K⁺ beyond the mean potential (Gillespie's
  point);
- past a charged band, `pb + csc` carries a field that local neutrality
  cannot (> 20×, the gate's question);
- the ruler on a tube = the 1-D Bernoulli series by hand (1e-7), and the
  neutral tube's ratio is D_Ca/D_K exactly.

## 3. Found

P_Ca:P_K in linear response (0.5 Å grid):

| | measured | 1-D donnan | 1-D csc | 1-D csc, reversal | 3-D neutral | slice | local | pb | slice + csc | local + csc | pb + csc |
|---|---|---|---|---|---|---|---|---|---|---|---|
| RyR1 9HEO | 7.0 | 0.50 | 0.87 | 0.64 | 0.63 | 0.72 | 0.75 | 1.01 | 0.95 | 0.85 | **1.08** |
| ITPR3 8TKF | 15.2 | 0.03 | 0.34 | 0.07 | 0.54 | 0.10 | 0.35 | 1.03 | 0.66 | 0.88 | **1.39** |
| ITPR3 7T3T | 15.2 | 0.06 | 0.63 | 0.40 | 0.84 | 0.33 | 0.41 | 1.41 | 1.19 | 1.05 | **1.79** |

Where Ca²⁺'s resistance lies under `pb + csc` (share of its drop across
± 3 Å of each constriction): 9HEO filter 3 %, gate 53 %; 8TKF 6 % / 31 %;
7T3T 9 % / 0 %.

- **The filter binds Ca²⁺ in 3-D too.** Under `pb + csc` the lumen's
  highest Ca²⁺ is 8.6 M in 9HEO against 1.2 M K⁺ (point ions under `pb`:
  5.6 against 5.9). In 8TKF it is 6.0 against 0.7.
- **The 3-D field does reach the gate, and it is worth about ×1.6 in RyR1**
  (neutral 0.63 → pb 1.01). The csc excess adds 7 % on top (1.08). In 1-D
  csc added ×1.75. In 3-D most of what csc did in 1-D, the field already
  does.
- **ITPR3 changes most**: the 1-D series crossing of the lysine rings
  (Donnan 0.03) disappears in 3-D (slice 0.10 → pb 1.03), as Round 7.11
  found for K⁺ g.
- **Every deposit stays 5–11× short**, and 9HEO's gate still holds over
  half of Ca²⁺'s resistance. The gate's field is present in 3-D. At this
  wall, placed this way, it is not strong enough.

**Xu's mutants** (`sel3d 9HEO --mutants`), P_Ca:P_K over the wild type
(each residue neutralised on all four subunits). The 1-D column is the
linear-response csc reading in the same bath (Round 7.17's reversal
reading gave 2.61):

| | measured | slice | local | pb | slice + csc | local + csc | pb + csc | 1-D csc |
|---|---|---|---|---|---|---|---|---|
| D4899Q | ×0.14 | ×0.95 | ×0.97 | ×0.89 | ×0.94 | ×0.92 | ×0.85 | ×0.91 |
| E4900N | ×0.64 | ×0.70 | ×1.02 | ×0.98 | ×0.93 | ×1.01 | ×1.00 | ×0.77 |
| D4938N | ×0.47 | ×0.96 | ×0.97 | ×0.84 | ×0.88 | ×0.93 | ×0.83 | ×0.84 |
| D4945N | ×0.93 | ×0.99 | ×1.00 | ×0.99 | ×0.98 | ×1.00 | ×0.98 | ×0.94 |
| E4955Q | ×1.19 | ×1.00 | ×1.00 | ×1.00 | ×1.00 | ×1.00 | ×1.00 | ×1.00 |
| Σ \|ln(model/measured)\| | | 2.93 | 3.34 | 3.05 | 3.10 | 3.24 | 3.02 | 2.80 |

No reading singles out D4899Q. In 3-D, D4899Q and D4938N cost about the
same (×0.85 and ×0.83 under `pb + csc`), while Xu measured ×0.14 and
×0.47. The field at the gate has no more to say about D4899 than the
filter's physics had.

## 4. Robust

`pb + csc` (and `pb`) P_Ca:P_K with one constant moved:

| constant | value | 9HEO pb | 9HEO pb + csc | 8TKF pb | 8TKF pb + csc |
|---|---|---|---|---|---|
| (registered) | | 1.01 | 1.08 | 1.03 | 1.39 |
| grid `pore3d.spacing` | 1.0 Å | 0.99 | 1.06 | 0.99 | 1.34 |
| `csc.structural_volume` | 0 | 1.01 | 1.04 | | |
| `csc.permittivity` | 40 | 1.01 | 1.12 | 1.03 | 1.64 |
| `permeation.permittivity_pore` | 78.4 | 1.19 | 1.21 | 1.20 | 1.56 |
| `pore_charge.smoothing` | 6 Å | 1.75 | 2.74 | 1.42 | 2.88 |
| `pore_charge.smoothing` | 1.5 Å, wall volume 0 | 0.84 | 0.75 | 0.75 | 0.62 |

At 1.5 Å with the wall's volume counted, `local_csc` refuses. The oxygens
alone pack 1.69 at their Gaussians' centres (9HEO 304 voxels, 8TKF 150),
where no fluid fits, and the refusal names the constants.

The ratio is converged in the grid and indifferent to the wall's volume.
Lowering the MSA's ε helps Cl⁻ as much as Ca²⁺ (P_Cl:P_K 0.14 → 0.66 in
9HEO, 1.41 in 8TKF, against Vais's 0.27 for ITPR3), so it is not
the missing selectivity. **The one constant that moves it is how far each
charge is spread**, and it moves it the way the gate question predicts.
At 6 Å the wall's charge reaches further toward the gate: 9HEO rises to
2.74, and the gate's share of Ca²⁺'s resistance falls from 53 % to 40 %.
At 1.5 Å it stays in the filter, and 9HEO falls to 0.75. Even at 6 Å the
ratio stays 2.6× short, and 6 Å is twice a side chain's positional
uncertainty.

## 5. What it means

The candidate Round 7.17 left was the gate's electrostatics. In 3-D that
field is present and reaches the gate. With the filter binding Ca²⁺ as
Gillespie's does, RyR1 reads 1.1 in linear response (about 0.8 at reversal,
by the 1-D translation) against 7.0. ITPR3 reads 1.4–1.8 against 15.2.
Neither the continuum's placement of the charge (7.11), the protein behind
the wall (7.13), the image cost (7.15), ion size (7.17), nor all of them at
the gate together (this round) closes the gap.

What remains is outside the continuum at this wall:
- **the gate's own geometry**: the deposits' gates are narrower than the
  channels that conduct the measured currents (Round 7.6's shortfall), and
  Ca²⁺'s resistance sits there;
- **charges the deposit does not place**: unresolved or mis-protonated
  groups nearer the gate;
- **the dehydration and binding** a mean field of spheres does not carry.
