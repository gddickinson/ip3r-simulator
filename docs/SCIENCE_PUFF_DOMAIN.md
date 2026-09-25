# SCIENCE — the puff microdomain (Round 7.2)

The park/drive cluster of `SCIENCE.md` § Puffs, with its instantaneous
mean-field Ca²⁺ replaced by Ca²⁺ pools that fill and drain, and puffs read
from fluo-4 fluorescence as the experiments read them.

## The model

**Pools** (`physics/microdomain.py`). These are Cao et al. 2014's three pools,
in cytosolic units, with every constant read from their model code (Text S1):

    cs  = γ2 (ct − (cb + b)/γ1 − c)                store (ER units)
    dc  = Jdiff + Jleak − Jserca + Jin − Jpm        cytosol
    dcb = γ1 (Jipr − Jdiff) − Jdye                  microdomain
    dct = Jin − Jpm                                 total
    db  = kon (B − b) cb − koff b                   fluo-4 bound (Cao 2013 Eq. 12)

    Jipr = k_ipr N_open (cs − cb),   Jdiff = k_diff (cb − c)

γ1 = 100 makes the microdomain relax at γ1 k_diff = 1000 s⁻¹. Fluo-4 lengthens
that to about 10 ms. The indicator sits in the microdomain, and its bound
Ca²⁺ is counted in the total. Without that, the dye would take Ca²⁺ from the
store: `test_store_balance_by_hand` fails if it is left out. Fluo-4's
constants are Shuai, Rose & Parker 2006's (k_on 150 µM⁻¹s⁻¹, k_off 300 s⁻¹,
K_d 2 µM, which is the K_d Cao 2013's text states; 25 µM total). Cao 2013's own
point-domain constants are in its Table S2, which PMC serves only behind a
browser challenge.

**What a receptor sees.** A closed receptor sees `cb`, and an open one sees its
mouth, `1.2 × cs`, the code's `120 (cs/100)`. That rule gives 539 µM at the
starting store. The mean-field model's `pd.ca_mouth = 120 µM` is the same rule
at a store of 100 µM. The store scan's "deep mouth" condition shows the
difference moves the open fraction by < 1 point.

**Store.** The code starts at ct = 45 µM, c = cb = 0.1 µM, so the store is
449 µM. That is not the pools' steady state, which sits at 630–1340 µM for
0.05–0.5 µM IP3 and takes minutes to reach. `rest_state` keeps the code's
store and puts the cytosol at its steady state against it (c = 0.095 µM at
0.2 µM IP3). Three clamps (`CLAMPS`): `store` holds cs (Cao 2013's
constant-flux assumption, and the default); `none` lets it deplete; `bath`
also holds the cytosol, so the microdomain exchanges with fixed surroundings
(the mean-field model's assumption, with the kinetics and dye kept).

**The cluster** (`physics/puffs_domain.py`). The receptors are
`puffs_pd.ReceptorCluster`, the same object the mean-field cluster now
drives. Each 0.1 ms step integrates the pools by RK4 at the step's number
open, then steps the receptors toward the gate equilibria at the new `cb`
and mouth. The refactor left the mean-field cluster bit-identical for a
given seed.

**Reading puffs** (`physics/puff_stats.py`). An event is a run of 1 ms bins
whose dF/F0 exceeds half the steady dF/F0 of one open receptor. Park-mode
flickers (0.3 ms) stay below that level; drive-mode openings cross it. A
blip is an event with one channel open. A puff is an event above 1.875 mean
blips, which is Cao's dF/F0 > 3 against a mean blip of 1.6, carried as a
ratio. IPIs run onset to onset. Thurley et al.'s density (Cao Eq. 14) is fitted by
maximum likelihood; `refractory_lr` is its likelihood ratio against the
exponential (χ²₁, so > 3.84 means a significant refractory period). The fit
recovers (λ, ξ) = (0.3, 0.8) from 2000 synthetic intervals, and an
exponential sample fits with ξ/λ > 20 and LR < 3.84 (tests).

## Calibrations (`tests/test_puffs_domain.py`)

- The rest state is steady under each clamp. With the store free, only the
  total drifts, at Jin − Jpm.
- Store balance by hand: d(cs)/dt from the bookkeeping equals
  γ2 (Jserca − Jleak − Jipr).
- One open receptor in a bath reaches exactly
  `cb = (k_diff c + k_ipr cs)/(k_diff + k_ipr)`, the dye reaches its
  equilibrium there, and the rise is the mean-field coupling (0.112 µM).
- Ca²⁺ is linear in the number open. dF/F0 is compressed more the further
  the microdomain passes K_d: 20 open compress it to 0.23 of linear at 4×
  release, and less at 1×.
- With no release and the mouth set to the microdomain, 400 receptors
  sample the stationary P_open of `park_drive` (within 15 %).

## Results (`python -m ip3r microdomain --scan …`)

**IPIs against the h42 recovery rate** (N = 10, 0.1 µM IP3, 900 s per point,
store clamped). Cao 2013 Fig. 4 reports λ 0.1–0.5 and ξ 0.5–2.2 s⁻¹ as a_h42
runs from 0.1 to 5 s⁻¹, with the shape turning exponential; its CV is 0.79
at a_h42 = 1 and 0.65–0.95 across the scan.

| a_h42 (s⁻¹) | puffs | rate (s⁻¹) | λ | ξ | LR | CV |
|---|---|---|---|---|---|---|
| 0.10 | 114 | 0.13 | 0.16 | 0.61 | 2.2 | 0.79 |
| 0.27 | 230 | 0.26 | 0.27 | 5.3 | 6.5 | 0.87 |
| 0.71 | 409 | 0.45 | 0.51 | 4.2 | 24.6 | 0.88 |
| 1.88 | 564 | 0.63 | 0.71 | 4.9 | 47.6 | 0.92 |
| 5.0 | 663 | 0.74 | 0.86 | 4.8 | 74.3 | 0.93 |

The puff rate rises 5.5× and the CV climbs toward 1, which matches Cao's
direction and magnitude. At the slow end (λ 0.16, ξ 0.61, CV 0.79) the
numbers land on Cao's, but the refractory fit is not significant on 114
intervals (LR 2.2). At the fast end ξ settles near 5 s⁻¹. That is a refractory
period of ~0.2 s, which a_h42 no longer sets: puffs take time to end, and
onset-to-onset intervals include them. The LR grows with the interval count,
so it measures confidence, not shape.

**Amplitude against cluster size** (0.2 µM IP3, 300 s per point). At Cao
2014's release (1×), puffs stay below the dye's K_d. The mean puff reaches
0.19–0.71 µM, the dye compresses F by 25 % at N = 25, and nothing bends.
Cao 2013's mean blip is dF/F0 1.6. Here it is 0.58–0.72 at 1× and
1.33–1.66 at 2.5× (`domain.blip_scale`, the scan's default). At 2.5×:

| N | 3 | 4 | 6 | 9 | 12 | 18 | 25 |
|---|---|---|---|---|---|---|---|
| dF/F0 | 3.86 | 4.36 | 5.33 | 6.58 | 8.14 | 9.22 | 10.43 |
| Ca²⁺ (µM) | 0.51 | 0.60 | 0.79 | 1.07 | 1.48 | 1.89 | 2.41 |
| peak open | 2.5 | 3.0 | 3.9 | 5.0 | 6.5 | 7.7 | 9.2 |

dF/F0 gains about 0.5 per receptor up to N = 12 and 0.17 beyond. Ca²⁺ gains
0.09–0.14 per receptor, then 0.07. So F bends at N ≈ 12, as in Cao's Fig. 8,
and bends harder than Ca²⁺. Ca²⁺ is not perfectly linear either, because a
bigger cluster recruits a smaller share (peak open 6.5 of 12, 9.2 of 25). At 4×
(puffs up to 3.6 µM) F/Ca falls from 6.6 to 3.4. The bend is a property of the
dye once puffs reach its K_d, as Cao and Solovey et al. argue. Whether it shows
depends on the release rate, which Table S2 fixes and we could not read.

**The sustained-open question** (Round 4; 0.2 µM IP3, N = 20, 60 s). Above
~0.5 µM coupling the mean-field cluster sat at a sustained 9 % open. The
release was raised 1–20× and run under five conditions, each adding one
thing to the one before:

| ×k_ipr | µM/open | mean-field | + deep mouth | + kinetics, dye (bath) | + cytosol (store held) | + store free | store min (µM) |
|---|---|---|---|---|---|---|---|
| 1 | 0.11 | 5.3 % | 4.5 % | 3.9 % | 5.5 % | 5.5 % | 443 |
| 2.1 | 0.24 | 8.5 | 7.5 | 7.2 | 7.7 | 7.5 | 371 |
| 4.5 | 0.50 | 9.3 | 8.5 | 8.1 | 6.5 | 7.3 | 265 |
| 9.5 | 1.06 | 9.1 | 8.4 | 7.5 | 4.0 | 5.8 | 174 |
| 20 | 2.24 | 9.0 | 8.4 | 6.5 | 2.1 | 4.1 | 139 |

In the bath, from 0.5 µM/open up, the sustained state survives: 6.5–8 % open,
a channel open 65–69 % of the time, Fano 1.3–1.8 against 2.7 at 1×. It is not an artefact
of instantaneous Ca²⁺ or of the shallow mouth. Freeing the store *raises*
activity against the held store. Depletion lowers the release and so the
Ca²⁺ that inhibits. What lowers activity is the cytosol filling: 27–60 µM
when one cluster releases into the whole cell, which is the model standing in
for a cell it is not. Even then, with the store free, a channel is open 51–65 % of the
time: sustained activity, not puffs. **Answer: the model, not the missing store.** At
microdomain Ca²⁺ of a few µM the park/drive receptor holds a steady
partial activity. The clamped stationary bell does not predict this level
(it gives 32 % at 0.11 µM/open), because an open receptor's own mouth
inhibits it, and the bell leaves that out.

## Limits

- Cao 2013's Table S2 was not read. Its release rate and its J_decrease
  (V_d c/(c + K_d)) are replaced by Cao 2014's linear exchange, which the
  2013 text says gives the same statistics. `domain.blip_scale` is matched
  to its blip, not read from it.
- λ_h42 follows Cao 2014's code (0.5 closed / 20 open s⁻¹), not Cao 2013's
  Eq. 10 (a_h42 + 100 c⁷/(c⁷ + 20⁷)). The a_h42 scan varies the closed rate.
- One seed per point. The a_h42 scan's ξ at 0.27 s⁻¹ (5.3) is out of line with
  its neighbours, which reads as noise on 230 intervals.
- The Puffs panel does not show the microdomain yet (Round 7.3).
