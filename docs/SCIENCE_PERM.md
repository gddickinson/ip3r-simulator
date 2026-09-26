# SCIENCE_PERM — ions through the pore

Split from `SCIENCE.md`: the permeation model and what it measures on the
ITPR3 panel. RyR1 uses the same model (`SCIENCE_RYR.md`).

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

## The shortfall, in three dimensions (Round 7.6)

`structure/pore_volume.py`, `physics/ohmic3d.py`, `physics/shortfall.py`;
`python -m ip3r shortfall [--scan]`.

**Why.** The 1-D model is short 2.4× even uncharged at its most generous
corner. Three candidates were named: 8TKF is a substate, the continuum
fails at 3 Å, or the exit is not where the profile window ends. Each is
tested here by one measurement. RyR1's 9HEO (5.9× short in 1-D) is the
control throughout.

**The instrument.** The same electrolyte (bath, diffusivities, the hard-
sphere exclusion of each ion's centre) now fills the pore's real shape. A
0.5 Å voxel is open when its centre clears every protein heavy atom by
vdW + ion radius. Inside the pore-domain span, only voxels within 25 Å of
the axis may conduct, because beyond the helices the atom model has empty
space where a cell has bilayer. Past the span, the box's side faces are bath.
Laplace's equation is solved over the component joining the two baths
(7-point finite volumes, preconditioned CG), and each species contributes
σ_s × g_s. Neutral only: this changes the geometry and nothing else.

**Calibrations** (`tests/test_ohmic3d.py`, 10 tests):
- A cylinder through a slab matches the series resistance plus Hall's two
  access terms (0.97 at 0.5 Å, 0.995 at 0.25 Å) and converges with the grid.
- Two pores conduct twice one.
- A blind hole changes nothing.
- A cap on the axis past the membrane, which a profile would call shut,
  conducts because current goes round it.
- A seal of r is a drawn r-Å pore.
- The hard-sphere test is exact for one atom.
- On 8TKF, seals of 20 and 25 Å agree, and 35 Å leaks through the lipid
  space. The leak is the failure the seal exists to prevent.

**Measured** (neutral K+, in-pore diffusivity 0.5 × bulk, each family's own
bath; 3-D extrapolated to h → 0 from 1.0 and 0.5 Å):

| deposit | 1-D (Round 4) | lumen-area 1-D | 3-D | bulk D, K+ 1 Å, 3-D | measured |
|---|---|---|---|---|---|
| 8TKF, ITPR3 activated | 65 | 135 | 106 | 261 | 358 / 545 |
| 7T3T, ITPR3 active (Schmitz 2022) | 64 | 110 | 85 | 278 | 358 / 545 |
| 9HEO, RyR1 open | 136 | 273 | 278 | 787 | 801 |

("Lumen-area 1-D" integrates dz/(σA) with A the real area of the in-plane
region around the axis, instead of the inscribed circle.)

- **A substate: not supported.** 7T3T comes from another laboratory and
  another preparation, and reads 85 pS against 8TKF's 106. Both are
  3.4–4.2× short of Mak's 358 pS, and RyR1's open deposit shows the same
  gap. One deposit's substate cannot explain a shortfall that two
  independent open ITPR3 structures and the RyR1 control share.
- **The exit window: no.** Moving the 1-D window from 6 to 24 Å past the
  span changes 8TKF by 2 % (64.7 → 63.2 pS); access is 2 % of the
  resistance. The 3-D solve has no window, and a box of 40 Å or a bath
  margin of 40 Å moves it < 0.5 %. No lateral exit bypasses the long,
  narrow stretch: its resistance is spread over ~50 Å at r_free 4–6 Å, and
  the filter slice holds about a third.
- **The continuum's geometry: yes, in part.** The inscribed circle
  undercounts the lumen, whose C4 corners the circle leaves out. The same
  electrolyte in the real shape conducts 1.3–2.0× more (8TKF 1.6×, 7T3T
  1.3×, 9HEO 2.0×). Most of that is area: the lumen-area 1-D reading is
  already 1.7–2.0× the circle. The 3-D solve then takes some back where
  current cannot use a corner fully.
- **What remains is the two unmeasured constants.** At the registered
  0.5 × bulk diffusivity, reaching 358 pS would need 1.7× bulk (8TKF) or
  2.1× (7T3T), which is impossible. At the sweep's favourable corner (bulk
  diffusivity, 1.0 Å K+ exclusion), **RyR1 reaches 787 pS against 801
  measured**, and ITPR3 reaches 261–278 pS: 1.3–1.4× short of Mak's
  358 pS and 2.0–2.1× short of Vais's 545. Round 4's 2.4× (1-D, same
  corner: 150 pS) was mostly the inscribed circle.

**What it means.** In the pore's real shape, a neutral continuum with
near-bulk mobility carries RyR1's measured conductance, and comes within
1.4× of the lower ITPR3 measurement. Two things are still missing for
ITPR3. First, the lower measurement and the higher one differ by 1.5×
between themselves (oocyte vs DT40 nuclei). Second, the wall charge has not
been tested in 3-D. The 1-D charged reading lowers the conductance, because
rings of alternating sign act as junctions in series. A cation-selective
pore that concentrates K+ should instead raise it. Whether the junctions
survive when the charges sit in the lumen's corners rather than in a
cylinder is the next question.

### Where the voltage falls (Round 7.10)

The 3-D solve also returns the potential φ in the lumen
(`physics/lumen_field.py`, `python -m ip3r lumen`; drawn in the viewer from
the Channel panel). It is read over S0's window for the bath's cation and
set beside the 1-D model's φ, which is the cumulative `∫ dz / A` of the
inscribed circle. Each curve is normalised to its own drop across the
window, so the comparison is of *where*, not how much.

**Calibration** (`tests/test_lumen_field.py`):
- a uniform cylinder drops linearly on both routes, and the share left
  outside the window is Hall's two access terms;
- a 6 Å neck of radius 2 in a radius-5 pore holds the share the 1-D sum
  gives by hand, on both routes, and pulls the half-drop point into itself;
- a blocked pore has no field and no drawn surface.

**Measured** (K⁺, neutral wall, ± 3 Å about each constriction):

| Deposit | In window | Half-drop z, 3-D / 1-D (Å) | Filter, 3-D / 1-D | Gate, 3-D / 1-D |
|---|---|---|---|---|
| 8TKF | 99 % | −83.2 / −83.8 | 32 % / 35 % | 12 % / 11 % |
| 7T3T | 99 % | −83.6 / −83.8 | 32 % / 35 % | 14 % / 17 % |
| 9HEO (RyR1) | 99 % | −78.1 / −79.4 | 24 % / 25 % | 22 % / 24 % |

**What it means.** The lumen's corners raise the conductance 1.3–2.0×, but
they do so almost evenly along the pore. Where the field falls moves by at
most 1.3 Å, and the filter's share by at most 4 points. So the 1-D model's
*profile* of the field is sound, even where its magnitude is not. A charge
or a blocker placed on the 1-D axis feels roughly the voltage it would feel
in 3-D. This does not settle the wall charge in 3-D. There the question is
the charges' distance from the ions, not the field of a neutral pore.

### The wall charge in 3-D (Round 7.11)

**Question.** In 1-D the lining charges of 8TKF *lower* its conductance,
because rings of alternating sign act as junctions in series: each carrier
must cross a zone where it is the excluded co-ion. Does that survive when
the charges sit where they are, in the lumen's corners, rather than smeared
over a cylinder?

**Method** (`physics/charge3d.py`, `physics/charged3d.py`,
`python -m ip3r wall3d [--scan] [--mutants]`). No 3-D drift-diffusion solve
is needed. Between identical baths at small voltage every species' flux is
`J = −(D c_eq/kT)∇μ`, divergence-free, with μ fixed at the baths. So the
slope conductance is the ions' *equilibrium* distribution followed by one
Laplace solve per species, with conductivity σ_s e^{−z_s u} per voxel
(Scharfetter–Gummel face weights), the species in parallel. The same
statement in 1-D (`linear_response_1d`) reproduces the existing Gummel
solver at 1 mV to 0.6 % on 8TKF and 9HEO, neutral, charged and paired.
(At the registered 20 mV the charged 1-D readings are 8–10 % nonlinear:
8TKF 33.0 pS at 20 mV, 30.1 at zero voltage. The 3-D readings are
zero-voltage, so they are compared with the zero-voltage 1-D numbers.)
The equilibrium comes from one of three placements, in order of how far
each departs from the 1-D model:

- `slice`: the 1-D charge per length (the same axial Gaussians, width
  `pore_charge.smoothing`) spread over each plane's real lumen region, and
  local Donnan per voxel;
- `local`: each group a 3-D Gaussian of that width about its own charge
  centre, over the lumen voxels within `charge3d.gaussian_reach` widths
  (charge conserved), and local Donnan per voxel;
- `pb`: the `local` density with nonlinear Poisson–Boltzmann in the lumen
  (ε = `permeation.permittivity_pore`). u = 0 on the bath voxels and no
  field into the protein, which is the ε_protein → 0 bound. Solved by Newton
  with a capped step and CG.

The electrostatics use the smallest ion's volume, and each species conducts
on its own volume. Ions are points, as in 1-D.

**Calibrations** (`tests/test_charged3d.py`). The weighted Laplace on a
tube equals its 1-D series to 1e-8, and a uniform energy scales g by
e^{−E}. Both placements conserve charge, and a charge that reaches no voxel
is reported. PB reaches Donnan (asinh) in a long charged tube to 0.1 %,
decays at the Debye length to 2 %, and obeys Gauss's law to 0.1 %. In a 3 Å
cylinder the slice closure *is* the 1-D reading (×0.027 against ×0.0267 for
two opposite rings), and local and PB stay within 1.3× of it: in a narrow
cylinder the junctions survive every closure. Ratios converge with the grid
(8TKF 1 → 0.5 Å moves each ratio by < 3 %).

**Measured** (h 0.5 Å; ratio to the same deposit's neutral 3-D reading):

| | 1-D | slice | local | pb |
|---|---|---|---|---|
| 8TKF, all lining charges (−8 e) | ×0.46 | ×0.80 | ×1.24 | ×3.08 |
| 8TKF, salt bridges paired (−4 e) | ×0.31 | ×0.45 | ×0.59 | ×0.52 |
| 7T3T, all (−16 e) | ×0.96 | ×1.44 | ×1.19 | ×4.24 |
| 7T3T, paired (−12 e) | ×0.42 | ×0.46 | ×0.54 | ×0.56 |
| 9HEO (RyR1), all (−32 e) | ×1.44 | ×1.33 | ×1.28 | ×1.93 |
| 9HEO, paired (−28 e) | ×1.23 | ×1.20 | ×1.08 | ×1.54 |

- **For the full wall, no.** Once the real cross-section dilutes the
  charge, 8TKF's junctions weaken (×0.46 → ×0.80). Once each group sits at
  its own centre, they reverse (×1.24), because a co-ion excluded from a
  corner still passes along the axis. Screened, the net −8 e raises g
  threefold (302 pS against 98 neutral). The sign is set by the placement,
  and its size by the smoothing width: over 1.5–6 Å, 8TKF runs
  ×0.39–5.07 (slice), ×1.21–5.45 (local) and ×1.97–5.65 (PB). The lining
  margin (1–8 Å) moves local and PB by < 10 %, and ε 20–80 moves PB
  ×2.4–4.0.
- **With D2478's salt bridges paired, yes.** In 8TKF and 7T3T the paired
  wall lowers g under PB at every width (×0.53–0.83) and under slice at
  every width (×0.23–0.46). Local lowers it at every width but the
  narrowest, 1.5 Å (×1.13 and ×1.06). Which reading is right is whether
  D2478 is charged. That is a question about the deposit, not the
  closure.
- **RyR1's mutants do not choose.** Xu 2006's ratios (h 1 Å), measured /
  1-D / slice / local / PB: D4899Q 0.20 / 0.90 / 0.89 / 0.80 / 0.78,
  D4938N 0.65 / 0.78 / 0.79 / 0.85 / 0.72, D4945N 0.92 / 0.90 / 0.95 /
  0.97 / 0.96, E4900N 0.63 / 0.76 / 0.86 / 0.97 / 0.94, E4955Q 1.01 /
  1.00 across. Summed |log| error: 1.89 / 2.02 / 2.14 / 1.90. Every closure
  misses D4899Q by 4×. So where the charge sits in the lumen is not what the
  point-ion continuum lacks at RyR1's filter. Round 7.4 reached the same
  conclusion from selectivity: the missing physics is charge–space
  competition (finite ion size), which remains the open item.

**What it means.** The 1-D result that "the ITPR3 lining charges lower the
conductance" is an artefact of the cylinder for the full wall, and holds only
if the filter's D2478 is neutralised by its bridge to R2471′. The charged 3-D
reading of 8TKF is not a number to quote: it spans 0.4–5.7× the neutral
reading over one unmeasured width, and its largest values (PB, 302–510 pS)
would reach Mak's 358 pS for reasons the calibration cannot support.

### Where the K+ drop falls with the charge (Round 7.12)

`physics/lumen_charge.py` reads 7.11's wall field on the drawn lumen. It has
two quantities. The first is the equilibrium wall potential u. The second is
K+'s electrochemical drop in linear response: the Laplace solve with
conductivity σ e^{−u}, whose integral is 7.11's K+ conductance. Its 1-D
counterpart is ∫dz/(A e^{−u}). In linear response each species has its own
drop, and the electrical potential would need Poisson at first order, so
the panel calls this the K+ drop and nothing more.

A cation well carries little of K+'s drop, because the ion is abundant
there. At 8TKF's filter the D2478 ring is such a well (−4 to −4.5 kT/e).
With the full wall, the filter's ±3 Å holds 1 / 7 / 18 % of the drop
(slice / local / pb), against 32 % neutral. Under slice and local the drop
moves to the K2482 / E2398 zone, 12 Å luminal. With D2478's salt bridges
paired the share is 23 / 29 / 38 %. 7T3T behaves the same way; in RyR1 9HEO
the drop is steepest at −76.4 Å under every reading. A plateau near one half
makes the half-drop point ill-conditioned (10 Å between two placements), so
the steepest point and each constriction's share are the reported numbers.

## Selectivity and the unitary Ca²⁺ current (Vais 2010)

`physics/selectivity.py`, `python -m ip3r selectivity [8TKF]`.

**Why.** The conductance shortfall is an absolute scale, and the in-pore
diffusivity is one of the two constants nobody has measured. A
permeability ratio does not depend on it: scaling every diffusivity by one
factor scales every flux by that factor, so the reversal potential does not
move (tested). Selectivity therefore tests the **wall charge**, the part of
the model that the conductance could not test cleanly.

**The experiments, as Vais et al. 2010 ran them** (lum-out nuclear patches of
rat InsP3R-3; every concentration registered under `selectivity.*`):

| Quantity | Cytosol (pipette) | Lumen (bath) | Measured |
|---|---|---|---|
| P_Cl : P_K | 140 mM KCl | 30 mM KCl + 110 mM NMDG-Cl | 0.27 ± 0.01 |
| P_Ca : P_K | 140 mM KCl | 140 mM KCl + 10 mM CaCl₂ | 15.2 ± 0.6 |
| i_Ca at 0 mV | 140 mM KCl, 3 µM Ca²⁺ | 140 mM KCl + 0.16 / 0.55 / 1.1 mM Ca²⁺ | 0.30 ± 0.02 pA/mM |

The reversal potential is the root of the solver's pore current. It is
read with Vais's Eq. 1 (general GHK), with P_Cl : P_K from the first
experiment as they used it. The model has concentrations, not activities,
so the ruler is applied to the concentrations the model saw. i_Ca is the
pore current at 0 mV, before the access correction; the error is about
0.1 mV across ~1 GΩ of access. Current is positive from lumen to cytosol.

**The solver change this needed.** With different baths, each mouth is now
in Donnan equilibrium with *its own* bath (before, the mean of the two).
That is Teorell–Meyer–Sievers. An impermeant bath ion (NMDG⁺), left out of
the species list, then sets a jump at its mouth even with no wall charge.
Symmetric baths give the same numbers as before (tested; the conductance
suite is unchanged).

**Calibrations** (`tests/test_selectivity.py`):
- the GHK ratio inverts GHK's own reversal potential;
- Vais's arithmetic is reproduced: the 0.30 pA/mM slope gives P_Ca =
  1.5 × 10⁻¹⁸ m³/s, and GHK from 545 pS and the ratios gives ~1.5 × 10⁻¹⁷
  at the ~104 mM activity they used;
- an uncharged cylinder gives Planck's liquid junction;
- a uniformly charged cylinder gives TMS at four charges of both signs;
- a −20 M wall approaches K⁺'s Nernst potential;
- an excluded NMDG⁺ gives the Donnan jump at √(30 × 140) mM plus Planck;
- the instrument *can* report Ca²⁺ selectivity. A 3.5 Å pore charged at
  −30 M along its whole length gives P_Ca : P_K 60 (8 at −3 M). The same
  wall concentrated in one 5 Å ring gives < 0.6 at any magnitude, because
  Ca²⁺ must cross the uncharged stretches with no Donnan enrichment.

**Measured on 8TKF** (the only conducting ITPR3 state):

| Wall | P_Cl : P_K | P_Ca : P_K | i_Ca (pA/mM) | g (pS) |
|---|---|---|---|---|
| neutral | 0.33 (0.74 with NMDG⁺ inside, slow) | 0.17 | −0.003 | 65 |
| charged (−8 e) | 0.01 | 0.00 | 0.0000 | 33 |
| paired (−4 e) | 0.05 | −0.07 | −0.0005 | 23 |
| acidic only (−16 e) | 0.00 | 0.69 | +0.046 | 174 |
| *measured* | *0.27* | *15.2* | *0.30* | *545* |

- **No reading comes within 20× of P_Ca : P_K = 15.2.**
- **The lining bases are Ca²⁺ barriers.** In the local-Donnan closure, K2529
  (four Lys at 11.6 Å radius in the cytosolic vestibule) becomes a +2.9 M
  wall at a +77 mV barrier. K2482 does the same at the luminal side. A
  divalent ion is excluded as e^{−2ψ}, so the channel passes almost no Ca²⁺.
  Neutralising K2529 alone is not enough (0.03). Neutralising both bases
  gives 0.69.
- **Even acidic-only, the charge is in rings.** E2398, D2478, D2518 and
  D2522 are four discrete rings, not a tract. So, per the calibration, the
  electroneutral model cannot give Ca²⁺ selectivity from them however
  strongly they are charged.
- **Anions.** Every charged reading makes the pore nearly anion-tight
  (P_Cl : P_K ≤ 0.05), but the channel passes Cl⁻ (0.27). The neutral
  pore's value is not a clean number either: it depends on whether NMDG⁺ is
  excluded at the mouth (0.33) or enters the pore slowly (0.74; the truth,
  exclusion at the constriction, lies between).
- **No GHK excess.** Vais's i_Ca is 8× (nominal concentrations) to 10×
  (activities) smaller than GHK predicts from their g and ratios; they
  attribute this to ion–ion and ion–channel interactions. The model
  obeys GHK to 4 % (acidic: 0.046 against 0.044 pA/mM). It has no
  interaction of that kind, so it cannot show the excess.

**What it means.** The conductance shortfall already said this continuum was
at its limit. Selectivity says what is missing. IP3R's Ca²⁺ preference
(15×) needs either a continuous charged tract, which 8TKF's lining does not
have, or the physics a point-ion, local-Donnan continuum leaves out:
ion size and crowding (charge–space competition), dielectric exclusion, and
screening within a wide vestibule, where the Debye length (5.8 Å at the model's
ε = 40; 8 Å in bulk water) is comparable to the radius. The last of these is exactly why K2529's ring
counts as a +2.9 M wall. The same comparison on RyR1, also Ca²⁺-selective, is the natural
control and is not yet run.

## Screening across the slice: the radial Poisson–Boltzmann closure

**Question.** Local Donnan spreads each slice's counter-charge uniformly over
the cross-section. At K2529 the slice is 9.9 Å wide (r_free; the charge
centres are at 11.6 Å) and R/λ_D = 1.7, so the ring should be screened
near the wall and leave a core the ions pass through. How much of the Ca²⁺
barrier is that closure?

**Model** (`physics/radial_pb.py`, `solve_pnp(..., closure="radial")`,
`python -m ip3r selectivity --closure radial`). In each slice, the
cylindrical Poisson–Boltzmann equation

  (1/r) d/dr (r dψ/dr) = −(F/ε) Σᵢ zᵢ c̄ᵢ e^{−zᵢψ/φ_T},  ψ′(0) = 0,  ψ′(R) = F X R / 2ε

with X the same fixed-charge density the Donnan closure uses (charge per
πR², R floored at the K⁺ radius) and c̄ᵢ the reservoir the slice is in
radial equilibrium with. Integrating over the disc gives Σ zᵢ⟨cᵢ⟩ + X = 0
for any c̄: the slice is neutral as a whole, as in Donnan, but its potential
varies across it. Each species sees its cross-section average
Γᵢ = ⟨e^{−zᵢψ/φ_T}⟩, which enters the axial Nernst–Planck as its own
potential wᵢ = −φ_T ln Γᵢ / zᵢ. A flat ψ makes every wᵢ the Donnan
potential, and the Gummel step then reduces to the Donnan loop's. The
reduction to 1-D assumes radial equilibrium and a wall charge that varies
slowly along z against R. The second is the weaker assumption, since the
charge is smoothed over 3 Å. Ions are points, as in the closure it replaces,
so the comparison isolates the closure.

Finite volume on 64 equal-width cells in r/R. The discrete Gauss law is
exact, and every slice's Newton system is solved in one banded call.
The permittivity (ε = 40, `permeation.permittivity_pore`) now enters the
answer, not only the reported Debye length.

**Calibration** (`tests/test_radial_pb.py`):

- R ≪ λ_D gives local Donnan in every species.
- A weak wall matches the Debye–Hückel cylinder, ψ = A I₀(κr) with
  A κ I₁(κR) = F X R / 2ε, to 0.1 %.
- Gauss's law holds to 10⁻⁸ in every slice.
- Ca²⁺ sits deeper than K⁺ in a wide negative wall's potential, which
  Donnan cannot show.
- 64 cells agree with 256 to 10⁻⁴ V.
- Near the Donnan limit, the solver gives the Donnan solve's conductance.
  A neutral pore is unchanged. On a wide charged pore the radial answer lies
  between the Donnan and uncharged ones.

**Measured (8TKF, charged reading, Ca²⁺ protocol at 0 mV).** Offsets are
per unit charge. The Ca²⁺ energy is twice the value shown.

| Ring | R (Å) | R/λ_D | X (M) | Donnan ψ (mV) | radial w_K (mV) | radial w_Ca (mV) |
|---|---|---|---|---|---|---|
| E2398 | 5.7 | 0.99 | −4.5 | −63 | −61 | −64 |
| K2482 | 6.7 | 1.15 | +3.9 | +82 | +71 | +67 |
| D2478 | 5.3 | 0.91 | −10.0 | −88 | −86 | −93 |
| D2518 | 4.4 | 0.76 | −17.9 | −98 | −96 | −103 |
| D2522 | 6.4 | 1.11 | −8.4 | −85 | −84 | −91 |
| K2529 | 9.9 | 1.71 | +2.9 | +78 | +60 | +55 |

The closure lowers K2529's Ca²⁺ barrier from +155 to +111 mV, a factor of
about 6 in partition. It lowers K2482's by 30 mV. That is not enough:

| Reading | P_Ca:P_K Donnan | P_Ca:P_K radial (ε 80 / 40 / 20 / 10) | g Donnan → radial (ε 40) |
|---|---|---|---|
| charged | 0.00 | 0.01 / **0.04** / 0.08 / 0.16 | 33 → 46 pS |
| charged − K2529 | 0.03 | 0.04 / 0.06 / 0.10 / 0.18 | — |
| acidic only | 0.69 | 0.69 / 0.69 / 0.69 / 0.70 | 174 → 172 pS |
| measured | 15.2 | | 545 pS |

The paired reading gives −0.04 (Donnan −0.07). P_Cl:P_K stays ≤ 0.07 on
every charged reading, against 0.27 measured.

**What it means.** Vestibule screening is real, but it is not the Ca²⁺
barrier. Even at ε = 10, where the ring is screened hardest, the charged
wall gives 0.16, a hundredth of the measurement. K2482 still sits in a
narrow slice (R/λ_D ≈ 1.1), and removing K2529 alone changes almost
nothing. The acidic rings are all in slices narrower than the Debye length,
where Donnan is already the right limit, so the closure leaves the
acidic-only reading exactly where it was (0.69). What the previous section
listed as missing narrows to three things: whether the two lysine rings
are charged at all (their pKa in that environment), ion size and crowding,
and dielectric exclusion.

## Protonation of the lining groups (Round 7.4)

**Question.** Every charged reading above gives each lining Asp/Glu −1 and
each Lys +1. After the radial closure, the charge state of the two lysine
rings (K2482, K2529) was the whole Ca²⁺ barrier, and full ionisation of the
eight D2518/D2522 carboxylates in a 4.4 Å lumen was the least plausible
assumption left. Are they charged at Vais's pH 7.3? And if the answer were
anything, could it matter?

**Two routes to the pKa, sharing no code or constants.**

- *Network* (`physics/pka.py`). Every Asp, Glu, His, Lys and Arg side
  chain within 20 Å of a lining group is a two-state site, with free energy
  G(θ) = Σ θᵢ ln10 (pH − pKaᵢ) + ½ Σ qᵢ Wᵢⱼ qⱼ. The model pKas are
  Thurlkill 2006's pentapeptides (Asp 3.67, Glu 4.25, His 6.54, Lys 10.40)
  and Fitch 2015's Arg (13.8). W is Coulomb with Mehler & Solmajer's
  sigmoidal ε(r), screened at the bath's Debye length. The sites are
  titrated together by Monte Carlo: single flips, pair flips, and a
  heat-bath move of each ring's four copies. The estimator averages each
  site's conditional probability. Desolvation and hydrogen bonds are left
  out, so for an acid the network's charge is the most that charge–charge
  coupling alone allows.
- *PROPKA 3* (`physics/pka_propka.py`, Olsson 2011). This is empirical,
  with desolvation and hydrogen bonds, run on the deposit's own atoms within
  25 Å of the lining groups (35 Å gives the same pKa). PROPKA gives each
  group one pKa with its neighbours in their default states, so it cannot
  show a ring titrating together. That is the network's part.

**Calibration** (`tests/test_pka.py`, `tests/test_protonation.py`). An
isolated site is Henderson–Hasselbalch exactly. A fixed neighbouring charge
shifts the pKa by exactly W/ln10. Two coupled acids follow their four-state
partition function. The Monte Carlo matches exact enumeration of a
12-site double ring at ε = 4 to 0.01. That test caught a sampler that
could not work: at ε = 4 a ring of four acids has two degenerate
half-protonated states, and single and pair flips stuck in one of them
(error 0.008–0.07 over seeds). The ring heat-bath move brings the error
to ≤ 0.0005. Also tested: the Bjerrum and Debye lengths against textbook
values, the network radius and PROPKA's context wide enough, the ring
copies in agreement, and Xu's Eq. 1 identical to GHK without a Cl⁻ term.

**8TKF at pH 7.3, 140 mM** (`python -m ip3r protonation`): mean charge per
lining ring.

| Ring | network | network ε 10 | network ε 4 | PROPKA pKa |
|---|---|---|---|---|
| E2398 | −1.00 | −0.89 | −0.50 | 4.6 |
| K2482 | +1.00 (pKa 12.9) | +1.00 | +1.00 | 10.6 |
| D2478 | −1.00 | −1.00 | −1.00 | 5.2 |
| D2518 | −1.00 (pKa 4.1) | −0.98 | −0.79 | 5.3–5.4 |
| D2522 | −1.00 | −0.99 | −0.64 | 4.5 |
| K2529 | +1.00 (pKa 10.9) | +1.00 | +1.00 | 10.4 |

Both routes keep every lining group charged at pH 7.3. The lysines stay
charged at any permittivity: the acid rings around them raise their pKa.
The acids lose charge only at a protein-like ε ≤ 10. The D2518 ring's
four carboxylates are 8.8 Å apart and R2524′ sits 4.3 Å from each, so they
repel each other much less than a 4.4 Å lumen suggests.

**Selectivity under each reading** (Vais's protocols):

| Reading | wall | P_Ca:P_K | P_Cl:P_K | g (pS) |
|---|---|---|---|---|
| formal | −8.00 e | 0.00 | 0.01 | 33.0 |
| network | −7.99 e | 0.00 | 0.01 | 33.0 |
| network ε 10 | −7.46 e | 0.00 | 0.02 | 31.4 |
| network ε 4 | −3.69 e | −0.01 | 0.03 | 24.2 |
| PROPKA | −7.92 e | 0.00 | 0.01 | 32.9 |
| *measured* | | *15.2* | *0.27* | *545* |

**A bound that needs no pKa** (`--corners`). Each of the six rings was set
formal or neutral, in all 64 combinations (about 2.5 min). The largest
P_Ca:P_K is 0.69, with both lysine rings neutral and every acid ring
charged: the "acidic" reading. With both lysines charged, which both routes
require, the largest is 0.05. Eight random fractional states inside the box
give at most 0.21. **No protonation state of 8TKF's wall comes within 20×
of 15.2.**

**The control: RyR1's open deposit** (9HEO) through Xu et al. 2006's own
experiment: 250 mM KCl, 10 mM CaCl₂ luminal, pH 7.4, read with their Eq. 1.
Their wild type measures 7.0.

| Reading | wall | P_Ca:P_K | g (pS) |
|---|---|---|---|
| formal | −32.0 e | 0.46 | 180 |
| network | −31.5 e | 0.45 | 179 |
| network ε 10 / 4 | −22.3 / −13.5 e | 0.31 / 0.07 | 160 / 119 |
| PROPKA | −28.7 e | 0.16 | 147 |

Xu's charge mutants, formal wall (P_Ca:P_K relative to wild type):

| Mutant | measured | model |
|---|---|---|
| D4899Q | ×0.14 | ×0.74 |
| E4900N | ×0.64 | ×0.24 |
| D4938N | ×0.47 | ×0.93 |
| D4945N | ×0.93 | ×0.97 |
| E4955Q | ×1.19 | ×1.00 (not lining) |

**What it means.**

- **The IP3R wall is charged as assumed.** Protonation is not why the
  model misses Vais's selectivity, and no pKa could make it so.
- **The failure is the continuum model, not the IP3R wall.** RyR1 is also
  Ca²⁺-selective, and its conductance model meets its mutants. Its
  selectivity is missed 15× under every reading, and the model's ranking
  of the mutants is wrong too: the filter's D4899Q should matter most, and
  the model makes E4900N matter most. A point-ion drift-diffusion pore
  with any electroneutral closure cannot produce either channel's Ca²⁺
  preference. The standing explanation for RyR1 is charge–space
  competition: finite ion size in a crowded, charged filter, which is
  Gillespie's model for the same mutants. That physics is absent here, and
  it is the next thing to add.
- **The mutants also test the pKa routes.** Xu's E4900N changes both the
  conductance (×0.63) and P_Ca:P_K (×0.64), so E4900 carries charge.
  PROPKA buries it (pKa 8.0, −0.18 e), and the network keeps it at −1.00.
  D4945N barely moves anything (×0.93), which the network's ε ≤ 10 bound
  (−0.19 e) allows and does not require.
