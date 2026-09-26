# Charge–space competition (Round 7.17)

`python -m ip3r csc [PDB ...] [--scan]`; `ip3r/physics/csc.py`,
`ip3r/physics/pnp_closures.py`, `ip3r/physics/csc_readings.py`;
`tests/test_csc.py`. Continues `SCIENCE_PERM.md` (Round 7.4's
protonation) and `SCIENCE_BORN.md` (Round 7.15).

## The question

Round 7.4 found the continuum pore short of both families' Ca²⁺
selectivity. The IP3R control, RyR1 (9HEO), gave P_Ca:P_K 0.46 under Xu
2006's protocol, against 7.0 measured. It also got the mutants' order
wrong: D4899Q came out ×0.74, where Xu measured ×0.14. Rounds 7.13 and
7.15 then ruled out the placement of the charges and the image cost.
What remained was the physics that Nonner, Catacuzzeno & Eisenberg 2000
and Gillespie 2008 use to explain the L-type and RyR filters. Ions have
size, and a crowded, charged fluid screens a divalent better than a
monovalent. Two questions:

1. Does the pore model, with that physics, bind Ca²⁺ in RyR1's filter as
   Gillespie's model does?
2. Does it then reach Xu's six P_Ca:P_K values and single out D4899Q?

## 1. The model

Each species gains a local excess chemical potential. It is a function of
the densities of every species in the slice (kT):

    mu_i = mu_i^HS + mu_i^MSA

- **Hard spheres** use the Boublik–Mansoori–Carnahan–Starling–Leland
  mixture. μ is the analytic derivative of its free energy through the
  ξ_n = (π/6) Σ ρ_j σ_j^n.
- **Screening** uses the mean spherical approximation (Blum 1975), in
  Nonner 2000's form. Γ and η are solved by damped iteration in every
  slice, at a uniform ε of 78.4 (Gillespie's).
- **The fluid.** The fluid has the permeant ions, with diameters twice the
  registered radii (K⁺ 2.76, Ca²⁺ 2.00, Cl⁻ 3.62 Å, Gillespie's own). Each
  lining acid becomes two oxygens of 2.8 Å carrying half its charge (a mean
  protonation lowers their valence, not their number). Each base is one
  2.8 Å sphere. Water is an uncharged 2.8 Å sphere at 55.5 M in the bath.
  The wall's groups have the fixed-charge map's own Gaussians and floored
  area, so their charge *is* the map (tested).
- **Their volume** (`csc.structural_volume`). The free radius already
  excludes the deposited carboxylate atoms. So counting the oxygens' volume
  in the lumen (1, Gillespie's tethered groups) counts it twice, and 0
  omits it. Both are read.

**The closure.** Each slice is in equilibrium with its reservoir under
local electroneutrality:

    c_i = a_i exp(−z_i ψ − μ_i(c)),   Σ z_i c_i + X = 0

with a_i the activity. Water satisfies the same equation with z = 0. This
is solved by Newton on (ln c, ln c_water, ψ) in every slice at once. The
Jacobian of μ comes from finite differences, and each step is capped
(`csc.max_step`) and halved where it would over-pack or not lower the
residual. Picard iteration on μ ran away at the RyR1 filter (100 M K⁺, no
water), because μ changes by ~10 kT between steps.

**In the drift-diffusion.** Each species carries the offset
w_i = ψ + μ_i kT / z_i e. This is the radial closure's structure, so both
now share one Gummel loop (`pnp_closures.solve_offsets`; the radial
closure moved there unchanged and passes its old tests). The ends are in
equilibrium with their own baths' activities. So the bath's own activity
coefficients are in the model, as they were in Xu's cuvettes. Their Eq. 1
reads concentrations, and so does this model's ruler.

**Calibrated first** (`tests/test_csc.py`):
- one component = Carnahan–Starling (to 1e-12);
- both μ = the derivative of their free energies (BMCSL; Blum & Høye's
  MSA energy + Γ³/3π), to 1e-6, including a half-charged species;
- the MSA's dilute limit is Debye–Hückel, and the restricted primitive
  model's closed form holds to 1e-10;
- with both terms off, one slice gives the Donnan potential, and the full
  Gummel loop the Donnan closure's current (1e-6); a bath returns itself;
- Gillespie's Fig. 8: in a −20 M slice, Li⁺-sized monovalents (1.33 Å)
  crowd Ca²⁺ out, and Cs⁺-sized ones (3.40 Å) let it in.

## 2. Found

**The filter binds Ca²⁺ as Gillespie's does.** Take 9HEO's most charged
filter slice (z −90.9 Å, −21.4 M: D4899 and E4900) at 150 mM KCl + 1 mM
Ca²⁺. It holds 10.3 M Ca²⁺ against 0.86 M K⁺ (packing 0.41). Ca²⁺'s
advantage over K⁺, split as his Fig. 7, is:

| | mean potential | screening | excluded volume |
|---|---|---|---|
| 9HEO filter, csc | +2.44 kT | +4.18 kT | +0.86 kT |
| Gillespie 2008 (RyR, PNP/DFT) | — | ~4 kT | ~0.5–1 kT |

In 8TKF's filter slice (−15.9 M) the terms are +2.26, +3.84 and +0.46 kT
(7.2 M Ca²⁺ against 1.5 M K⁺).

**Yet P_Ca:P_K stays below 1.** Xu's protocol (250 mM KCl, 10 mM CaCl₂
luminal), with × over the wild type:

| | measured | donnan | csc | csc, wall volume 0 |
|---|---|---|---|---|
| wild type | 7.0 | 0.46 | 0.64 | 0.47 |
| D4899Q | ×0.14 | ×0.74 | ×0.85 | ×0.95 |
| E4900N | ×0.64 | ×0.24 | ×0.70 | ×0.86 |
| D4938N | ×0.47 | ×0.93 | ×0.75 | ×0.57 |
| D4945N | ×0.93 | ×0.97 | ×0.86 | ×0.72 |
| E4955Q | ×1.19 | ×1.00 | ×1.00 | ×1.00 |
| Σ \|ln(model/measured)\|, mutants | | 3.54 | 2.61 | 2.83 |

The K⁺ conductance is 172 pS under csc (180 donnan; 801 measured). The
mutant g ratios are ×0.79–0.94 against ×0.20–1.01, and D4899Q comes out
×0.94 against ×0.20. Charge–space brings E4900N's ratio in (×0.70
against ×0.64). It narrows the summed error but does not single out
D4899Q, which is among the mildest mutants in the model and by far the
worst in Xu's recordings.

**Why: the gate is in series.** Take the linear-response resistance of
each ion about equilibrium in the mixed bath,
1/P_i = ∫ exp(z_i w_i) / (D_i A_i) dz. Under csc the filter holds 12 %
of K⁺'s resistance (4 % under Donnan) and 0 % of Ca²⁺'s: the filter
itself is now Ca²⁺-selective. But the gate window (±3 Å about Q4933,
r_free ≈ 3.5 Å, uncharged) holds 45 % of K⁺'s and 72 % of Ca²⁺'s. In an
uncharged segment both ions sit at bath concentration. So the segment's
own ratio is (D_Ca/D_K) × the area ratio = 0.54, by hand from the
profile, and it caps the pore: the whole pore reads 0.64 at reversal and
0.87 in linear response. Under the local-neutrality closure no charge
reaches the gate. So it is geometry and diffusivity, not the filter's
physics, that sets RyR1's P_Ca:P_K here. In 8TKF, Ca²⁺'s resistance lies
elsewhere (the lysine rings, as Round 7.4 found), and csc moves P_Ca:P_K
only from 0.00 to 0.07 against 15.2.

**What would change it.** The 9HEO gate is narrower than the channel
that conducts 801 pS (Round 7.6's shortfall). A real gate also sits
within a Debye length of D4938 and the filter's rings, whose field the 1-D
local-neutrality closure cannot carry into an uncharged slice. Widening
the gate stretch (z −84 to −70 Å) by hand to 8 Å raises csc only to 0.75
(a first trial, with the wall inferred from the net charge), because the
widened stretch is still uncharged and in series. The candidate is
therefore the gate's electrostatics: the 3-D wall field (Rounds
7.11–7.15), which does reach into the gate, carrying the csc excess and
read as a permeability ratio.

## 3. Robust

`csc 9HEO --scan`: the wild type's P_Ca:P_K, and D4899Q over it, under csc:

| constant | values | wild type | D4899Q |
|---|---|---|---|
| `csc.structural_volume` | 0 / 0.5 / 1 | 0.47 / 0.54 / 0.64 | ×0.95 / ×0.92 / ×0.85 |
| `csc.permittivity` | 40 / 60 / 78.4 | 0.31 / 0.56 / 0.64 | ×0.54 / ×0.79 / ×0.85 |
| `csc.water_diameter` | 2.5 / 2.8 / 3.1 Å | 0.68 / 0.64 / 0.75 | ×0.78 / ×0.85 / ×0.91 |
| `csc.oxygen_diameter` | 2.4 / 2.8 / 3.2 Å | 0.59 / 0.64 / 0.91 | ×0.89 / ×0.85 / ×0.64 |

The wild type stays at 0.31–0.91 against 7.0, and D4899Q at ×0.54–0.95
against ×0.14, under every scanned constant.
