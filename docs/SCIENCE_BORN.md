# The image cost of the low-ε wall (Round 7.15)

`python -m ip3r born [PDB ...] [--scan] [--mutants]`;
`ip3r/physics/born3d.py`, `ip3r/physics/born_readings.py`;
`tests/test_born3d.py`. Continues `SCIENCE_BRIDGE.md` (Rounds 7.13–7.14).

## The question

Every closure so far treats an ion as a point in a mean field. A real ion
also polarises what surrounds it. Near a wall of lower permittivity than
the water, that polarisation pushes it away: the image force (Parsegian
1969). Round 7.13 named it as the one electrostatic term every closure
omits. It was also the cheapest step toward the charge–space candidate for
RyR1's D4899Q, which no placement of point charges reaches (×0.79 against
Xu's ×0.20). Two questions:

1. How large is the cost in these pores, and what does it do to K⁺
   conductance, neutral and charged?
2. Does it single out D4899Q among Xu 2006's mutants?

## 1. The self-energy

An ion's image cost at r is half its own reaction potential:

    W(r) = ½ z² [u_het(r) − u_bath(r)]          (kT)

Here `u_het` is the potential at r of a unit charge at r in the deposit's
permittivity map, linearly screened by the salt wherever an ion centre
reaches. `u_bath` is the same charge in bulk water with the bath's
screening. Both are solved on the same grid, so the grid's self-term
cancels. W depends on geometry and the bath only, never on the wall
charge. So one field serves every reading of a deposit, and it is cached
(`data/cache/born`, keyed by the inputs' hash).

- **Local solves.** The salt screens the reaction field, so each lumen
  voxel is solved on its own cube of ±`born.box_half_width` (12 Å), with
  the potential zero on the faces. A voxel farther than `born.reach`
  (10 Å) from every low-ε voxel takes W = 0 unsolved. That leaves about
  10⁵ solves per deposit at 1 Å: 13–14 min on ten processes, with each
  cube's operator poured into a fixed sparsity template and CG started
  from the bulk solution.
- **The map.** Water is the region an ion's own sphere sweeps: every voxel
  within the ion's radius of a voxel its centre reaches (membrane sealed,
  as for the lumen). No ion centre is then closer than its radius to the
  low-ε region (tested on 8TKF). Round 7.13's closure put the boundary at
  the ion *centres*, which would put every wall-contact ion inside the
  protein's dielectric. A first try with a 1.4 Å water probe did the same
  in crevices a K⁺ reaches but a water probe does not (W 20–30 kT there).
  The dielectric closure takes the same map when read with W
  (`surface="swept"`), so the mean field and the image see one wall.
- **How W enters.** Each species' Boltzmann factor gains e^{−z²W}, both
  in the Poisson–Boltzmann equilibrium (`dielectric_pb(self_energy=)`)
  and in the linear-response conduction (`WallField.energy`). A symmetric
  salt with W alone stays neutral (tested), so the neutral pore with W is
  one Boltzmann-weighted Laplace solve per species.
- **Screening** is the bath's. Screening by the mean field's own ions
  against a bath reference would also fold the Debye–Hückel activity of a
  dense counter-ion cloud into W. That is charge–space's term, not the
  wall's. A one-step try gave a number that mixed the two (8TKF ×1.24
  against ×3.42), so it is not reported.

**Calibrated first:**
- The unscreened planar image, l_B Δ/(4 ε_w d), by hand; the vacuum
  Bjerrum length by hand.
- The screened planar wall against Debye–Hückel (the half-space Green's
  function's reflected part, by quadrature): within 5 % at the registered
  box and 1 Å (+0.4 % at 2.5 Å, −3.2 % at 3.5 Å), and within 2 % at
  0.5 Å in a 16 Å box.
- Zero in a uniform medium; negative (attraction) when the wall's ε is
  higher; the box cut 12 vs 20 Å within 3 %; the planar wall at the reach
  is < 0.02 kT; the solver tolerance 1e-7 vs 1e-10 < 1e-4 kT; parallel =
  serial; the template = the operator.
- W alone leaves a symmetric salt neutral, and it weakens the screening of
  a fixed charge (the potential deepens).

## 2. Found

On the axis (unit charge; ×4 for Ca²⁺), and K⁺ g against each deposit's
neutral 3-D g without W (1 Å; IP3R in 140 mM, RyR1 in 250 mM KCl):

| | 8TKF | 7T3T | 9HEO |
|---|---|---|---|
| W at the filter, axis | 1.19 kT | 1.54 kT | 0.77 kT |
| W at the gate, axis | 0.76 kT | 0.86 kT | 0.66 kT |
| neutral (pS) | 90 | 79 | 220 |
| **neutral + image** | **×0.26** | **×0.20** | **×0.28** |
| dipole on 7.13's map | ×3.18 | ×3.47 | ×3.81 |
| dipole on the swept map | ×3.38 | ×3.58 | ×3.86 |
| **dipole + image** | **×3.42** | **×3.58** | **×2.18** |
| base omitted + image | ×6.48 | ×6.58 | ×2.77 |
| pair omitted + image | ×1.99 | ×1.49 | ×2.07 |
| dipole + image (pS) | 309 | 284 | 481 |
| measured (pS) | 545 / 358 | 545 / 358 | 801 |

- **The image costs the neutral pore 3.5–5×.** About 1 kT on the axis at
  the filter, and more off it, is enough to cut an uncharged pore's
  conductance to a fifth or a quarter. So the absolute shortfall of Round
  7.6 would be worse, not better, for an uncharged wall.
- **It costs a charged ITPR3 pore nothing.** Under the dipole the lumen is
  a K⁺ tract (Cl⁻ carries 0–1 %). Where the fixed charge dominates, the
  counter-ion density is pinned by neutrality, and the image cost is paid
  by the potential, not the concentration. The local Donnan form says it
  directly: c₊c₋ = c²e^{−2W} and c₊ − c₋ = |X|, so c₊ ≈ |X| whatever W is,
  and only the co-ion is lost. 8TKF moves ×3.38 → ×3.42, and 7T3T not at
  all. Moving the boundary from the ion centres out by an ion radius (the
  swept map) moves the dipole 1–6 %.
- **RyR1 is not pinned.** 9HEO's dipole falls ×3.86 → ×2.18 with the
  image. Its wall charge does not hold the counter-ions along the whole
  path, so where it does not, both ions pay e^{−W}. The model's 481 pS is
  then 0.60 of Xu's 801 pS, against 1.06 without W.

### Xu 2006's mutants (9HEO, each residue neutralised on all four subunits)

| mutant | measured | dipole (7.13) | swept | + image |
|---|---|---|---|---|
| D4899Q | ×0.20 | ×0.79 | ×0.79 | **×0.48** |
| D4938N | ×0.65 | ×0.74 | ×0.74 | **×0.28** |
| D4945N | ×0.92 | ×0.98 | ×0.98 | ×0.96 |
| E4900N | ×0.63 | ×0.88 | ×0.88 | ×0.73 |
| E4955Q | ×1.01 | ×1.00 | ×1.00 | ×1.00 |
| summed \|ln error\| | | 1.91 | 1.91 | 1.92 |

- **The image makes every charge matter more, but not the right one.**
  D4899Q moves about a third of the way toward Xu's ×0.20 on a log scale, and
  E4900N comes within 16 %. But D4938N, the cytosolic-side ring, falls to
  ×0.28 against ×0.65, and becomes the model's largest effect. The summed
  error is unchanged (1.91 → 1.92). Xu's order puts D4899Q first by 3×.
  The image term only scales the mean field's dependence on each charge by
  how narrow the pore is where that charge sits. D4899's effect needs
  something that acts at the filter alone.
- So the image cost is not what the continuum lacks at the filter. What is
  left is ion size (charge–space competition), which is local to the
  crowded filter by construction.

## 3. Robust to the method's own choices (8TKF)

`python -m ip3r born 8TKF --scan`. Each row is a fresh W field:

| setting | filter W (axis) | neutral + image | dipole, swept | dipole + image |
|---|---|---|---|---|
| box ±8 Å | 1.01 kT | ×0.32 | ×3.38 | ×3.55 |
| **box ±12 Å (registered)** | **1.19 kT** | **×0.26** | **×3.38** | **×3.42** |
| box ±16 Å | 1.25 kT | ×0.24 | ×3.38 | ×3.39 |
| reach 6 / 10 / 14 Å | 1.19 kT | ×0.26 | ×3.38 | ×3.42 |
| ε protein 2 | 1.35 kT | ×0.22 | ×3.39 | ×3.60 |
| ε protein 10 | 0.85 kT | ×0.38 | ×3.38 | ×3.19 |

- The box is not quite converged in the narrow pore: 12 → 16 Å adds 5 %
  to the filter's W and takes the neutral pore from ×0.26 to ×0.24. So the
  neutral cost is if anything under-read. The charged pore moves < 1 %.
- The reach changes nothing at two decimals, so 6 Å would do (a GUI
  could use it to solve faster).
- The protein's permittivity sets W's size (0.85–1.35 kT) and the neutral
  cost (×0.22–0.38). The charged pore stays within 7 % of its no-image
  reading at every ε. The pinning result does not depend on it.

## What this does not settle

- W exceeds 10 kT on 0.4 % of 8TKF's lumen voxels (631 of 152,534; 2.7 %
  exceed 5 kT). They sit in dead-end pockets at the box's side, 27–30 Å
  off the axis, where a voxel's cube runs past the grid and is filled by
  repeating its edge. They carry almost no current; the axis values above
  are 1–2 kT.
- The image is screened by the bath's ionic strength. In the ITPR3 filter
  the counter-ion density is several molar, which would screen it further.
  That only strengthens "costs the charged pore nothing".
- The in-pore water's ε (`permeation.permittivity_pore`, 40) holds in the
  bath too. A Born cost of moving from ε 80 to 40 would be uniform along
  the lumen and is not included.
- Ca²⁺ pays 4W: about 5 kT on 8TKF's axis at the filter against K⁺'s 1.2.
  A neutral pore would lose Ca²⁺ selectivity. A charged one pins Ca²⁺ as
  it pins K⁺. The selectivity readings (1-D, Round 7.4) do not include W.
