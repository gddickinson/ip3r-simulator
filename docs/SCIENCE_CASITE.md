# A Ca²⁺ site that blocks K⁺ (Round 7.27)

`python -m ip3r casite [PDB ...] [--region BAND] [--kind KIND]
[--no-required]`; `ip3r/physics/ca_site.py`, the `hidden` and `coupling`
energies of `ip3r/physics/pnp3d.py`; `tests/test_ca_site.py`. Continues
`SCIENCE_WALLSEARCH.md` (Round 7.25).

## The question

Round 7.25 showed that Vais 2010's pair (P_Cl:P_K 0.27, P_Ca:P_K 15.2) is
not a mean-field property of any wall in these lumens. For point ions in
series, (P_Cl:P_K)(P_Ca:P_K)² over the uncharged pore's is at most 1, and
the pair needs about 3,200. A Ca²⁺-only well was the one lever that moved
the right ratio: it leaves Cl⁻ alone. It peaked at 4.5 (8TKF) and 6.9
(7T3T), for two reasons. The Ca²⁺ it gathers is uncompensated charge that
repels further entry, and its linear-response ceiling is 46× / 28× the
uncharged pore, because the access regions keep the rest of Ca²⁺'s
resistance. The emergent item asked about two interactions the well
lacks: a site compensated as it fills, and K⁺ blocked by the Ca²⁺ it
holds.

## The model

**The site.** `casite.sites` sites (4, one per subunit) are spread evenly
over the band's lumen voxels (the membrane span; density s = 0.27 M on
8TKF, 0.30 M on 7T3T). Each binds one Ca²⁺. The affinity is given as the
empty site's pull d (kT), so K_d = s / (e^d − 1). Binding is fast and
bound Ca²⁺ moves with the free, as in 7.25's well. The total is then
c_f (1 + s / (c_f + K_d)), which is a Ca²⁺ energy
μ = −ln(1 + s / (c_f + K_d)) that saturates as the site fills. In Slotboom
form the free Ca²⁺ is c_f = n e^{−2ψ} whatever μ is. μ is therefore a
coupling: `pnp3d.steady_state` recomputes it from each Gummel iterate,
damped by `pnp3d.coupling_damping`, and requires it to settle with u. The
occupancy is θ = c_f / (c_f + K_d).

* **Compensated**: each bound Ca²⁺ comes with −2e fixed, so a site is
  neutral full or empty. μ is then a `hidden` energy: transport sees it,
  Poisson counts only c_f.
* **Uncompensated**: Poisson counts the bound Ca²⁺ as well. This is 7.25's
  well, made saturable.

**The block.** A K⁺ crossing the band finds a fraction θ of it held by
Ca²⁺ and is stopped by an occupied site with probability f
(`casite.block`, 1 = single file). In the mean field this is the energy
−ln(1 − fθ) on K⁺ over the band. Poisson counts the K⁺ it excludes. It is
the anomalous-mole-fraction mechanism in its simplest continuum form.

Neither term acts on Cl⁻. Vais's Cl⁻ experiment holds no Ca²⁺, so it
reads the uncharged pore's P_Cl:P_K exactly (8TKF 0.291, 7T3T 0.367;
measured 0.27). The wall is otherwise uncharged, and the ions are points,
as in 7.25's search.

## Calibration (`tests/test_ca_site.py`, on a tube)

- A hidden energy leaves Poisson's u identical to the no-site solution
  (to 1e-12). The Ca²⁺ in the band is e^d times the bath's. The same
  energy made visible moves u by 0.05 kT or more.
- At equilibrium, uncharged, the compensated site holds c(1 + s/(c + K_d))
  and θ = c/(c + K_d), both checked by hand.
- A dilute site (many sites, trace Ca²⁺) is 7.25's fixed well, to 1e-5.
- Away from equilibrium, uncompensated and blocking: the state's energies
  are what the coupling returns for that state (within 2e-4 kT).
- The Cl⁻ experiment's fluxes and u are bit-identical with and without a
  site.
- At reversal, the compensation and the block each raise P_Ca:P_K over the
  plain site, and the two together raise it further than either.
- 8TKF pinned: a compensated, blocking 8 kT site is nearly full
  (θ > 0.95) and passes 15.2, while P_Cl:P_K stays the uncharged pore's.

## Results

The site covers the membrane span at 1 Å. d runs from 2 to 12 kT in steps
of 2. The last column is the depth at which P_Ca:P_K crosses 15.2 (found by
root-finding to 0.02 kT). Raw output is in `data/casite/`.

P_Ca:P_K at Vais's reversal (the uncharged pore gives 0.257 on 8TKF and
0.516 on 7T3T):

| d (kT) | 2 | 4 | 6 | 8 | 10 | 12 | 15.2 at |
|---|---|---|---|---|---|---|---|
| **8TKF** uncompensated | 1.66 | 4.15 | 4.97 | 4.37 | 4.03 | 3.97 | — |
| compensated | 2.05 | 6.13 | 8.52 | 8.99 | 9.06 | 9.07 | — |
| uncompensated + block | 1.80 | 5.01 | 6.31 | 5.51 | 5.06 | 4.97 | — |
| compensated + block | 2.33 | 11.40 | 36.1 | 55.5 | 60.1 | 60.8 | **4.41 kT** |
| **7T3T** uncompensated | 2.91 | 6.62 | 7.04 | 5.62 | 4.92 | 4.77 | — |
| compensated | 3.49 | 9.56 | 12.60 | 13.18 | 13.27 | 13.28 | — |
| uncompensated + block | 3.10 | 7.67 | 8.50 | 6.75 | 5.88 | 5.70 | — |
| compensated + block | 3.85 | 15.31 | 44.2 | 68.0 | 73.8 | 74.7 | **3.98 kT** |

Occupancy θ (the band's mean at the Ca²⁺ reversal) rises with d in the
same way for every kind: 8TKF 0.05–0.07 at 2 kT, 0.23–0.41 at 4, 0.54–0.86
at 6, and above 0.96 from 10 kT. At the crossing, the compensated blocking
site on 8TKF has θ = 0.50 (2.0 of 4 sites held),
V_Ca = +18.17 mV. On 7T3T it has θ = 0.34 (1.3 held), V_Ca =
+17.47 mV.

**Uncompensated sites still fall with depth,** with or without the block.
They peak at 6 kT (8TKF 4.97 / 6.31, 7T3T 7.04 / 8.50) and fall as they
fill with charge that repels further Ca²⁺. This is 7.25's finding; a
saturable site changes the height of the peak only a little.

**Compensation removes the fall but not the ceiling.** The compensated
site rises with depth and plateaus once full, at 9.07 (8TKF) and 13.28
(7T3T). That is 35× / 26× the uncharged pore, close to 7.25's
linear-response ceilings for any Ca²⁺-only energy over the span (46× /
28×). Access holds the rest of Ca²⁺'s resistance. A Ca²⁺ affinity alone,
however well compensated, stops short of 15.2 in both lumens.

**The block is what crosses.** On a compensated site, blocking K⁺ passes
15.2 at d = 4.41 kT on 8TKF and 3.98 kT on 7T3T, and reaches 61 / 75 at
full occupancy. P_Cl:P_K stays at the uncharged pore's value, so the
score (|ln| distance to both ratios) at the grid's nearest point is 0.36 /
0.31. Round 7.25's best was 1.30 / 0.90, and the uncharged pore scores
4.15 / 3.69. On an uncompensated site the block adds only 20–30 %, because
that site never fills far at reversal.

**What the crossing means.** K_d = s / (e^d − 1) is 3.3 mM on 8TKF and
5.7 mM on 7T3T. That is a low-millimolar Ca²⁺ affinity in the pore. At
Vais's 10 mM luminal Ca²⁺ it leaves the site partly occupied at reversal,
and each occupied site then stops K⁺.

## What it means

Round 7.25 found that no mean-field wall gives the measured pair. This
round finds the smallest addition that does. It needs two things
together:

1. **Ca²⁺ occupancy that blocks K⁺.** An affinity alone is capped near
   the linear ceiling (9.1 / 13.3). Blocking K⁺ lowers P_K in the Ca²⁺
   experiment only, and nothing in the series bound limits that, because
   the bound assumes each ion moves independently.
2. **A site that stays neutral as it fills.** Without compensation the
   site's own charge keeps occupancy, and so the block, low.

This is the anomalous-mole-fraction picture of Ca²⁺-selective channels,
here as its minimal continuum form. It makes a testable prediction:
measured P_Ca:P_K should depend on the luminal Ca²⁺ concentration, and
K⁺ current through ITPR3 should fall at millimolar luminal Ca²⁺ more than
the drop in K⁺ activity explains. The model does not say where the site
is. Spread over the whole span it works; 7.25 showed that a site at one
constriction covers too little of the path. The 4 sites, the block f = 1
and the band are hypotheses, each a registered parameter. The finding is
that the pair is reachable once Ca²⁺ and K⁺ interact through occupancy,
and not before.
