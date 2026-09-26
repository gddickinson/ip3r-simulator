# The lining salt bridge: charged, and what its field does (Round 7.13)

`python -m ip3r bridge [PDB ...] [--scan] [--mutants]`;
`ip3r/physics/bridge_charge.py`, `ip3r/physics/dielectric3d.py`;
`tests/test_dielectric3d.py`. Continues `SCIENCE_PERM.md` (Rounds 7.4, 7.11).

## The question

8TKF's four D2478 carboxylates line the filter. Each is bridged
(N–O 2.45–2.58 Å) to R2471 of the next subunit. Every salt bridge that
touches the lining is this one (7T3T: 3.5 Å; RyR1 9HEO's homologous
D4899–R4892′: 2.9 Å). The guanidinium sits 3.6 Å further from the axis at
the same height (r 7.5 → 11.1 Å, z −90 Å). So the pair is a radial dipole
behind the filter wall. Round 7.11's 3-D readings turned on how the pair
was counted: PB ×3.1 with D2478 alone ("full"), ×0.5 with the pair removed
("paired"). Neither is a physical state.

## 1. Is D2478 protonated? No

Two independent routes (Round 7.4's), each run with the partner in and out.
Apparent pKa ranges over the four copies, at Vais's pH 7.3 (Xu's 7.4 for
RyR1):

| reading | 8TKF D2478 | 7T3T D2478 | 9HEO D4899 |
|---|---|---|---|
| network, sigmoidal ε | 2.0–2.2 | 2.5 | 4.1 |
| ... base left out | 4.0 | 4.0 | 5.1 |
| network ε 10 | 1.1–1.3 | 1.5–1.6 | 5.5 |
| ... base left out | 4.5 | 4.5 | 7.4 |
| network ε 4 | −3.6 to −3.1 | −2.3 | 6.7 |
| ... base left out | 4.7–4.8 | 4.4–4.5 | 10.3–10.4 |
| PROPKA 3 | 5.2 | 4.9–5.0 | 5.0 |
| base (Arg) pKa, all readings | 12.3–19.3 | 12.2–19.3 | 12.9–19.4 |

The ITPR3 acid is ionised (≤ −0.99 e) under every reading, with or without
its partner. The partner lowers its pKa by 2–8 units, and the arginine
never titrates. RyR1's D4899 is the weaker acid: it would be neutral at
ε 4 without R4892′, and with it is still −0.8 e or more. So **the pair is
two charges**, and the question is where their field falls.

## 2. The dielectric closure

The `local` and `pb` closures (Round 7.11) cannot represent a charge behind
the wall. They normalise each lining group's Gaussian over the lumen voxels
it reaches, and `pb` lets no field into the protein. The new `dielectric`
closure solves nonlinear Poisson–Boltzmann on the **whole box**:

    −∇·(ε ∇u) = (F² h² / ε0 R T) (Σ z c e^{−z u} [lumen] + X)

- ε is the pore water's (`permeation.permittivity_pore`, 40) on voxels
  open to the ion probe, with the membrane sealed as for the lumen. It is
  `dielectric.eps_protein` (4; Schutz & Warshel 2001) elsewhere. A face
  takes the harmonic mean of its two voxels.
- The ions live in the conducting lumen only. u = 0 on its bath faces, and
  every other box face is insulating.
- X is every modelled Asp/Glu/Lys/Arg in the box (`scope="all"`), or the
  lining alone (`"lining"`). Each is a Gaussian of `dielectric.charge_width`
  (1 Å) about its own charge centre, over every voxel it reaches: a buried
  charge stays buried.
- Newton, with the step capped by its change in the lumen (the protein's
  part is linear). The slope conductance is Round 7.11's linear response
  on this u.

**Calibrated first:**
- A charge in a uniform medium gives Coulomb's l_B/(εr) to 3 %.
- Across an ε 40 | ε 4 interface, ε·E is continuous to 1e-6.
- No charge gives no field. Charge is conserved in lumen and protein alike.
- With the protein's ε → 0 and the charge in the lumen, the solve is
  Round 7.11's `pb` to 0.02 kT/e.
- A base 3.6 Å behind an acid at a 3 Å pore's wall cancels 0.51 of the
  axis potential in uniform ε. The Coulomb ratio is 0.46, and the salt
  screens the farther charge a little. In ε 4 it cancels less (0.64).
- The grid: 8TKF's dipole reads ×3.18 / 3.11 / 3.11 at 1.0 / 0.75 /
  0.5 Å, so readings are at 1 Å (about 20 s each).

## 3. Found

K⁺ g over the same deposit's neutral 3-D g (h 1 Å):

| reading | 8TKF | 7T3T | 9HEO |
|---|---|---|---|
| `pb` (7.11, lining in the lumen) | ×3.17 | ×4.24 | ×2.00 |
| dielectric, lining only | ×3.85 | ×5.59 | ×2.79 |
| **dielectric, all: the pair as a dipole** | **×3.18** | **×3.47** | **×3.81** |
| ... base omitted ("full") | ×5.30 | ×5.77 | ×4.29 |
| ... pair omitted ("paired") | ×1.79 | ×2.01 | ×3.74 |

- **The pair counts about half.** In ITPR3 the dipole sits midway (on a
  log scale) between the two limits Round 7.11 had to choose between. In
  8TKF it equals 7.11's unpaired `pb` (×3.17), and it is far from paired
  PB's ×0.5.
- **The sign is settled.** Once the rest of the protein's charge is in the
  field, the wall raises g under every reading. The box holds −24 e in
  8TKF, −28 e in 7T3T and −64 e in 9HEO. Round 7.11's "paired PB lowers
  g" came from a wall with its lining pairs removed and nothing else
  charged.
- **Robust to the closure's own choices** (8TKF, dipole): ε_protein 2–20
  gives ×3.18–3.25, charge width 0.5–3 Å gives ×3.17–4.00, water ε 20–80
  gives ×2.67–3.81, a box of 30 → 40 Å gives ×3.18 → 3.19, and a bath margin of 25 → 35 Å gives ×3.18 → 3.17. The order
  pair omitted < dipole < base omitted holds in every row.
- **D4899Q is still missed.** Xu 2006's mutants under the new closure
  (9HEO, ε_protein 4) give D4899Q ×0.79 against ×0.20 measured (`pb`
  ×0.78), and E4900N ×0.88 against ×0.63. No placement of point charges,
  in the lumen or behind it, reaches Xu's D4899Q. Charge–space
  competition (finite ion size at the filter) is the one candidate left.

## What this does not settle

- Born (image) repulsion of an ion near the low-ε wall is left out here,
  as in every closure. It would lower every reading, and most in the
  narrowest planes.
- The partner's charge enters at its deposited position. A 3–4 Å rotamer
  change would move it, and the width scan is the only proxy for that.
- The absolute conductance is still the continuum's. This round moves
  ratios, not the 1.3–2.1× shortfall of Round 7.6.
