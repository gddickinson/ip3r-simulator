# Selectivity at bi-ionic reversal in 3-D (Round 7.23)

`python -m ip3r reversal [PDB ...] [--mutants] [--scale F ...] [--spacing H]`;
`ip3r/physics/pnp3d.py`, `ip3r/physics/reversal3d.py`;
`tests/test_reversal3d.py`. Continues `SCIENCE_SEL3D.md` (Round 7.19) and
`SCIENCE_GATE.md` (Round 7.21).

## The question

Rounds 7.19 and 7.21 read P_Ca:P_K from linear response: identical baths,
a small voltage, each species' conductance in its equilibrium field. Every
open wall read at most 1.8, gate or no gate. Xu 2006 (RyR1, 7.0) and Vais
2010 (ITPR3, 15.2) measured something else: the voltage at which no current
flows when the two baths differ, converted to a ratio by GHK. The two
readings agree only for ions that move independently in a constant field.
A filter that binds Ca²⁺ can set the reversal potential without setting the
conductance ratio. In 1-D the two already differ (9HEO csc 0.87 in linear
response, 0.64 at reversal). The protocol was the last difference between
model and measurement left untested in 3-D.

## The solve: steady Poisson–Nernst–Planck on the voxels

At reversal the net current is zero but each species still carries
current, so neither the concentrations nor the potential are at
equilibrium. `pnp3d.steady_state` solves that state on Round 7.19's voxels:

- **Nernst–Planck in Slotboom form.** With E_i = z_i ψ + μ_i (kT; μ_i the
  held charge–space excess), J_i = −D_i e^{−E_i} ∇n_i with
  n_i = c_i e^{E_i}. Steady state is ∇·(e^{−E_i}∇n_i) = 0, which is Round
  7.11's weighted Laplace problem with Scharfetter–Gummel faces. With n_i
  held at each bath's value (lumen c_lum, cytosol c_cyt e^{z v}), the
  solution is n_lum + (n_cyt − n_lum) φ_i, and the particle current is
  D_i g_i (n_cyt − n_lum). Each species moves in its own volume (its own
  radius). In electrostatic voxels it cannot reach, it takes n from its
  nearest voxel.
- **Poisson.** ψ = v φ_0 + u, where φ_0 is the uncharged lumen's Laplace
  potential and u = 0 on both baths. Laplace's operator annihilates v φ_0,
  so u obeys Round 7.11's Poisson–Boltzmann with each species' charge
  n_i e^{−z_i ψ − μ_i}: the applied field, the excess and ln n_i all enter
  its `offset`. The protein carries no field (as `pb`).
- **Gummel.** The two are alternated, the other held, until u moves
  < `pnp3d.gummel_tolerance` (10⁻⁴ kT/e). Undamped, every solve here
  converged.
- **Reversal.** `reversal` finds the root of the net current by Brent's
  method. Each solve starts from the nearest one already found. It takes
  6–7 solves per experiment. The grid is `reversal3d.spacing` = 1 Å.

**The protocols** are the 1-D model's (tested equal):
- `xu` (RyR1): 250 mM KCl on both sides, 10 mM CaCl₂ in the lumen, read
  with Xu's Eq. 1 (GHK with P_Cl taken as 0);
- `vais` (IP3R): P_Cl:P_K from 30 mM KCl + 110 mM NMDG-Cl in the lumen
  against 140 mM KCl, then P_Ca:P_K from 10 mM luminal CaCl₂ in 140 mM KCl,
  read with that P_Cl:P_K. NMDG⁺ is confined to the luminal bath below the
  membrane span (the 1-D model excludes it at the luminal mouth).

**Readings**: `neutral` (no wall charge; Poisson is still on, because the
junction needs it), `pb` (Round 7.11's local placement, point ions) and
`pb + csc` (Round 7.19's: each species' excess from the locally neutral
fluid, held). The csc reference bath is the side that holds every permeant
species and is neutral without NMDG⁺: the lumen for Xu and for Vais's Ca²⁺,
the cytosol for Vais's KCl. Outside the wall groups' reach the excess is
zero on both sides, which leaves out the two baths' difference in activity
coefficient. That difference is ≤ 0.012 kT (`bath_gamma_difference`,
tested < 0.02). The local-neutrality closures have no non-equilibrium form
here.

**Calibrated first** (`tests/test_reversal3d.py`):
- at v = 0 between identical baths the solve is Round 7.11's PB
  equilibrium to 10⁻¹⁰ kT/e, and carries no flux;
- at 0.1 mV each species' flux is Round 7.19's linear-response
  conductance (2 × 10⁻⁵);
- an uncharged tube reaches Planck's liquid junction as the Debye length
  shrinks against the gradient: 0.6 mV off at 0.03/0.14 M over 60 Å,
  0.02 mV at 0.3/1.4 M over 160 Å;
- a charged tube reaches Teorell–Meyer–Sievers: at 300 Å it is 0.1 mV off
  for a cation-selective wall and 1.0 mV for an anion-selective one, and
  doubling the length halves the miss (the ends' double layers);
- an excluded NMDG⁺ sets the Donnan jump and Planck from there: 0.2 mV off
  with its region 5 Å deep (1.5 mV at 20 Å, its own polarisation);
- with no root inside the bracket, the search refuses by name.

On a real deposit the uncharged pore obeys GHK. 9HEO's neutral reversal,
read with its own linear-response P_Cl:P_K, gives 0.608 against its
linear-response 0.607. Xu's ruler, which takes P_Cl as 0, reads the same
reversal as 0.064, because a neutral pore passes Cl⁻. The ruler assumes a
cation-selective pore.

## Found

P_Ca:P_K, 1 Å grid. "Linear" is Round 7.19's reading on the same grid.

| deposit | measured | 1-D csc, linear | 1-D csc, reversal | neutral (own P_Cl) | pb reversal | pb linear | pb + csc reversal | pb + csc linear |
|---|---|---|---|---|---|---|---|---|
| RyR1 9HEO | 7.0 | 0.87 | 0.64 | 0.61 (0.61) | **0.95** | 0.99 | **0.90** | 1.06 |
| ITPR3 8TKF | 15.2 | 0.34 | 0.07 | 0.26 | **0.98** | 0.99 | **1.15** | 1.34 |
| ITPR3 7T3T | 15.2 | 0.63 | 0.40 | 0.52 | **1.41** | 1.36 | **1.54** | 1.68 |

V_rev (mV, pb + csc): 9HEO +1.72 (Xu); 8TKF −38.98 (KCl/NMDG) / +3.59
(Ca²⁺); 7T3T −39.14 / +4.63.

- **The protocol does not rescue the ratio.** In 3-D the reversal reading
  lies within 15 % of linear response on every deposit and reading, and
  mostly below it (9HEO 0.90 against 1.06, 8TKF 1.15 against 1.34). The
  1-D model's large difference between the two (8TKF 0.34 → 0.07) came
  from its local-neutrality closure's barriers, not from the protocol.
  The filter still binds Ca²⁺ at reversal (9HEO 8.2 M under `pb + csc`,
  8TKF 5.8 M). Binding in the filter lowers the reversal reading slightly;
  it does not raise it.
- **Xu's ruler is not what is short.** Read with the model's own P_Cl:P_K
  instead of 0, 9HEO's charged readings rise only 10 % (pb 1.04, pb + csc
  1.18).
- **The charged wall over-excludes Cl⁻.** Vais's P_Cl:P_K is 0.27. The
  neutral 8TKF pore reads 0.29 and 7T3T 0.37, while every charged 3-D
  reading gives 0.001–0.003 (the 1-D model ≤ 0.07). A mean-field negative
  wall that shuts out Cl⁻ this completely still selects Ca²⁺ over K⁺ only
  1–1.5-fold. The measured channel does both: it passes Cl⁻ and
  selects Ca²⁺ 15-fold.

**Xu's mutants at reversal** (`reversal 9HEO --mutants`), P_Ca:P_K over the
wild type:

| | measured | pb | pb + csc | (7.19 linear, pb + csc) | (1-D csc reversal) |
|---|---|---|---|---|---|
| D4899Q | ×0.14 | ×0.85 | ×0.80 | ×0.85 | |
| E4900N | ×0.64 | ×0.96 | ×0.98 | ×1.00 | |
| D4938N | ×0.47 | ×0.73 | ×0.66 | ×0.83 | |
| D4945N | ×0.93 | ×0.97 | ×0.95 | ×0.98 | |
| E4955Q | ×1.19 | ×1.00 | ×1.00 | ×1.00 | |
| Σ \|ln(model/measured)\| | | 2.84 | 2.67 | 3.02 | 2.61 |

The reversal spreads the mutants a little further apart, but in the wrong
order. D4938N costs the most, D4899Q does not stand out, and Xu measured
the opposite.

## Robust

`pb` / `pb + csc` at reversal, with Round 7.19's linear reading beside:

| change | 9HEO pb | 9HEO pb + csc | 8TKF pb | 8TKF pb + csc |
|---|---|---|---|---|
| (registered, 1 Å) | 0.95 (0.99) | 0.90 (1.06) | 0.98 (0.99) | 1.15 (1.34) |
| `pore_charge.smoothing` 6 Å | 1.80 (1.77) | 2.70 (2.84) | 1.32 (1.38) | 2.36 (2.80) |
| grid 0.5 Å | 0.96 | 0.91 | | |

The 1 Å grid is within 1.5 % of 0.5 Å. With the charge spread to 6 Å the
ratio rises about as far as it did in linear response (Round 7.19), and
reversal again reads at or just below it.

**The wall's charge scaled** (`reversal 8TKF --scale …`, `pb`, point
ions, Vais's two experiments). This is the path a mean-field wall of 8TKF's
shape can take, from uncharged to twice the deposit's charge:

| scale | 0 | 0.05 | 0.1 | 0.25 | 0.5 | 1 | 2 | *measured* |
|---|---|---|---|---|---|---|---|---|
| P_Cl:P_K | 0.29 | 0.13 | 0.066 | 0.018 | 0.007 | 0.003 | 0.001 | *0.27* |
| P_Ca:P_K | 0.26 | 0.43 | 0.58 | 0.83 | 0.95 | 0.98 | 0.86 | *15.2* |

P_Cl:P_K falls as soon as the wall is charged. It halves by a twentieth of
the deposit's charge. P_Ca:P_K never reaches 1 and turns down past the
deposit's own charge. No scale comes near the measured pair.

## In the viewer: each ion at reversal (Round 7.24)

Channel panel → lumen box → **Steady state: at reversal** solves one
experiment of the family's protocol under one reading
(`physics/lumen_reversal.py`; headless `python -m ip3r reversal PDB
--lumen [READING ...] [--experiment Ca2+|Cl-]`). The surface is cut from
the reversal grid's electrostatic volume (1 Å, Ca²⁺'s radius), and is
coloured by one ion's

- **concentration** `c = n e^{−(zψ + μ)}`, the charge Poisson counts, on
  a fixed log scale (`display.lumen_conc_min`–`display.lumen_conc_max`,
  1 mM–10 M);
- **electrochemical drop** `(n − n_lumen)/(n_cytosol − n_lumen)`, 0 at the
  luminal bath and 1 at the cytosolic. With flux `−D e^{−E} ∇n`, this is
  the ion's own series resistance accumulated from the lumen.

Calibrated on tubes: an uncharged 0.3/1.4 M junction's concentrations are
Planck's straight lines and electroneutral (to 1 % of the range); across an
acidic band each ion's drop equals ∫e^{E}dz accumulated, normalised, to
2×10⁻⁴ (the other ion's resistance misses by 0.72). The counter-ion K⁺
drops 0.2 % across the band and the co-ion Cl⁻ 72 %.

The Ca²⁺ experiment (Vais's for 8TKF and 7T3T, Xu's for 9HEO), shares of
each ion's drop within ± 3 Å of each constriction:

| deposit, reading | V_rev | Ca²⁺ peak (z) | Ca²⁺ filter / gate | Ca²⁺ steepest | K⁺ filter / gate | Cl⁻ filter / gate |
|---|---|---|---|---|---|---|
| 8TKF neutral | −0.09 mV | 0.01 M | 28 / 11 % | −84.9 | 31 / 13 % | 49 / 10 % |
| 8TKF pb + csc | +3.59 mV | 5.8 M (−89.9) | 7 / 31 % | −76.9 | 24 / 20 % | 45 / 7 % |
| 9HEO neutral | +0.13 mV | 0.01 M | 25 / 20 % | −75.9 | 31 / 20 % | 27 / 30 % |
| 9HEO pb + csc | +1.72 mV | 8.2 M (−90.9) | 3 / 60 % | −75.9 | 22 / 35 % | 28 / 18 % |
| 7T3T neutral | +0.37 mV | 0.01 M | 26 / 13 % | −83.9 | 30 / 15 % | 37 / 14 % |
| 7T3T pb + csc | +4.63 mV | 6.2 M (−90.9) | 6 / 0 % | −82.9 | 25 / 4 % | 27 / 30 % |

(filter z −85.9 / −86.9 / −86.9 Å; gate −76.4 / −76.9 / −66.4 Å.)

**What it shows.** The wall gathers Ca²⁺ to 6–8 M a few ångström luminal of
the filter, 600–800× its bath, and that well carries almost none of Ca²⁺'s
drop (3–7 %, from about a quarter in the neutral pore). Ca²⁺'s resistance
sits in the uncharged stretch cytosolic of the well: at the gate in 9HEO
(60 %) and 8TKF (31 %), and just past the filter in 7T3T, whose gate is
20 Å away and wide. Round 7.19 found this in linear response, and Round
7.21 found that widening the gate only moves it. At reversal the picture
holds, and the viewer shows it directly: P_Ca:P_K is set where the
deposits place no charge.

## What it means

GHK's two readings of permeability, linear response and reversal, agree
to 15 % in the 3-D charged lumen. Neither approaches 7.0 or 15.2. The
protocol was the last difference between model and experiment that the
continuum left untested at this wall, and it is closed with a negative
result. Together with Rounds 7.11–7.21, no electrostatic or charge–space
treatment of the deposits' own wall charge gives the measured Ca²⁺
selectivity, however it is placed, screened, sized, imaged, gated or read.

Left outside the model:
- **charges the deposit does not place**: unresolved lining residues,
  groups in the cytosolic vestibule, the lipid headgroups at the mouths;
- **specific binding and dehydration** that a mean field of spheres does
  not carry. The Cl⁻ reading points here. A wall that selects Ca²⁺ by its
  mean charge excludes Cl⁻ in proportion. The measured channel passes Cl⁻
  at a quarter of K⁺ and still selects Ca²⁺ 15-fold, which a single mean
  field cannot do.

*Tested in Round 7.25 (`SCIENCE_WALLSEARCH.md`), negative:* no point-ion
wall can pass Cl⁻ and select Ca²⁺ together, because the series bound
(P_Cl:P_K)(P_Ca:P_K)² ≤ 1 holds over the uncharged pore. A Ca²⁺-only well
that spans the whole membrane reaches only 4.5 (8TKF) / 6.9 (7T3T) at
reversal.
