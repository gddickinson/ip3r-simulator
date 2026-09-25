# SCIENCE_EC — excitation–contraction coupling at the couplon

The couplon of Stern, Pizarro & Ríos 1997: two rows of RyR1s in the
junctional cleft, half opened by voltage sensors (V channels) and half by
their own Ca²⁺ (C channels). The C-channel gating, the cleft and the Mg²⁺
readings it uses are in [`SCIENCE_RYR.md`](SCIENCE_RYR.md).

## The V channels: the couplon under voltage clamp (Round 6.6)

**Sources.** Rios, Karhanek, Ma & Gonzalez 1993 (J Gen Physiol 102:449,
supplied as a PDF): the release channel opposite a tetrad of voltage
sensors, written as an MWC allosteric protein (`physics/allosteric_v.py`).
The channel is closed or open as a whole. Each of four sensors moves
independently (k_c, k_−c = α/2·exp(±(V − V̄)/8K), Eqs. 2–3). Each active
sensor multiplies the opening rate by 1/f and the closing rate by f. That
gives ten states and the closed forms of Eqs. 6–7 for Po(V) and Q(V). Stern
et al. 1997 used this model for their V channels, "with time-dimensioned
rate constants increased by a factor of 2", at 0.1 pA against the C
channels' 0.3 pA. They made the V channels "neither activated nor
inactivated by Ca²⁺" (`physics/couplon.py`).

**Which row.** Rios's Table I has four reference fibres. Fiber 827 is the
one their kinetic figures use, and it is the only one that gives Stern's V
plateaus (Figs. 11–12: 0.457 / 0.052 / 0.003 at 0 / −30 / −50 mV; the
others give 0.24–0.71 at 0 mV). It also reproduces Rios's own statement that
Po saturates at 0.72 with every sensor active.

**The factor of 2, read from the figures.** Stern's Fig. 11, their
stand-alone check of the V model, follows the printed rates (rise rms 0.020
at ×1, 0.087 at ×2, and the decay after repolarisation agrees). Their
Fig. 12, the couplon runs, follows ×2 (V at 10 ms: 0.43, against 0.42 at ×2
and 0.26 at ×1). So Fig. 11 was made before the factor was applied, and the
couplon uses ×2 (`ec.stern_rate_scale`). Both figures print "0.0" where
0.1 belongs on the y axis, so the values were read from the tick positions.

**Exact simulation.** The V channels do not depend on Ca²⁺, so each one is
simulated first by Gillespie's method through the voltage steps (a step
boundary is exact because the process is memoryless). The C array is then
simulated by Gillespie's method with each V opening and closing as a
scheduled change of Ca²⁺. `cleft.v_coupling_matrix` gives the V → C
coupling from the same finite-volume solve. At 0.1 pA a V channel gives
6.5 µM at the C channel across the row and 5.6 µM at each neighbour along
it. The Monte Carlo mean of the V channels equals the master equation
through the step and after it (Stern's own check). The V → C coupling is
reciprocal to the C → V coupling within 5 %.

**Calibration: Stern's couplon reproduced** (`python -m ip3r ec`, 200
couplons, 100-ms steps from −90 mV). Stern's figure values are in brackets.

| step | V plateau | C peak | C plateau | flux peak/plateau |
|---|---|---|---|---|
| 0 mV | 0.459 [0.457] | 0.644 at 5 ms [0.64, ~4 ms] | 0.117 [0.105] | 2.7 |
| −30 mV | 0.053 [0.052] | 0.381 at 11 ms [0.39] | 0.088 [0.078] | 3.8 |
| −50 mV | 0.003 [0.003] | 0.112 [0.123] | 0.042 [0.039] | 2.6 |

After repolarisation the C channels shut (Po ≤ 0.005), which is Stern's
"control". At −50 mV single C events last a median 11 ms (peak 7 open), and
those reaching half the array last 20 ms, the same as the triggered cleft
sparks of Round 6.3. The couplon is Stern's.

**The question.** Round 6.5 left two candidate rows. Mg²⁺ at the
activation site alone was physiological if the voltage sensor lifts the
inactivation-site block (Laver 2018), but it was 5–70× too slow. Mg²⁺ at
both sites ended sparks in 6 ms, but only 6 of 30 channels could open. In
Stern's couplon the V channels carry no Mg²⁺ or Ca²⁺ action at all, which
is exactly the lifted block. So the C channels were given the scheme fitted
to Murayama's bell under the fibre's 1 mM Mg²⁺, and the couplon was
stepped as above (200 couplons; K_Mg,A 769 µM, Meissner's reading):

| C channels | C peak / plateau at 0 mV | flux peak/plateau (0 mV) | C after repolarisation (−30 mV) | −50 mV events |
|---|---|---|---|---|
| Stern 1997 | 0.64 / 0.12 | 2.7 | 0.004 | 11 ms; half-array ones 20 ms |
| fitted, no Mg²⁺ | 0.76 / 0.73 | 1.01 | 0.72 | 199 of 210 never end |
| fitted, Mg²⁺ activation site | 0.68 (at 54 ms) / 0.66 | 1.03 | 0.39 | 74 of 91 large never end |
| fitted, Mg²⁺ both sites | 0.12 / 0.11 | 1.05 | 0.000 | 2 ms, one channel; none large |

The other two K_Mg,A readings bracket it and change nothing qualitative
(60 couplons). At 54 µM the activation site alone leaves the C channels
nearly shut (Po 0.012 at 0 mV). At 521 µM they reach 0.60 with no peak (200 couplons at 0 mV: peak/plateau 1.02)
and 0.03–0.05 is still open after repolarisation.

**What it settles.** Real V channels confirm Round 6.5 and sharpen it.
With the C channels fitted to the measured bell, **no Mg²⁺ arrangement gives
the release waveform.** Either the C array amplifies and does not stop when
the membrane repolarises (no Mg²⁺, or Mg²⁺ at the activation site only), or
it stops but does not amplify (Mg²⁺ at both sites, where 80 % of the C
channels are inactivated at rest). The flux never peaks: its peak/plateau
ratio is 1.0–1.1, against 2.7 with Stern's constants and the transient peak
measured in frog fibres. The peak in Stern's couplon comes from Ca²⁺
inactivation that switches on during the pulse at the tens of µM the cleft
holds (Ki 10 µM). The fitted inactivation (Ki 249 µM) is too weak to do
that, and Mg²⁺ at the inactivation site only adds a fixed block present
before the pulse. A time-dependent terminator is what is missing:
inactivation stronger than the steady-state bell shows, or SR depletion
under the couplon (Rios et al. correct their release for 50–60 % depletion
by a single conditioning pulse, their Fig. 2).

**Limits.** V channels have no Ca²⁺ or Mg²⁺ action (Stern's
simplification). Stern chose the V and C unitary currents (0.1, 0.3 pA) to
fit the peak/plateau ratio, so they are not independent evidence. There is
no Iγ (Ca²⁺ feedback onto the sensors) and no global cytosolic Ca²⁺.
Luminal depletion is Round 6.7, below. Rios's model is from frog at ~10 °C; Murayama's bell is
rabbit RyR1 at 25 °C.

## SR depletion as the terminator (Round 6.7; `physics/lumen.py`)

**The question.** Round 6.6 left the scheme fitted to Murayama's bell
without a terminator that switches on during the pulse. Emptying of the SR
was the leading candidate: Ríos et al. 1993 correct their release for
50–60 % depletion by one 100-ms pulse to +20 mV.

**The model is Stern's own** (their Fig. 20, "effects of global calcium
dynamics"). The SR is one well-mixed pool; Stern cites Shirokova & Ríos
1996 for mixing on the junctional time scale. Every unitary current, V and
C, is proportional to the content left. The content falls with the
ensemble's mean release, at a couplon density of 4.8 µm⁻³ (1 pA for 1 ms
removes 24.9 µM of myoplasmic-equivalent Ca²⁺), from Stern's 2 mM. Uptake
is first order in the deficit, with τ = 0.43 s, digitised from Stern's
Fig. 20 C recovery after the pulse (421 and 450 ms at 0 and −30 mV). The
pool is shared, so its path is made consistent with the ensemble it drives
by damped fixed-point iteration (3–10 iterations to 1 %). Given the path,
the couplon is still exact: a change of content is a scheduled event.

Left out: Stern's global cytosolic Ca²⁺ (a µM boundary), any luminal Ca²⁺
sensor on RyR1 (no source), and gradients inside the SR. Stern's pump is
driven by global Ca²⁺ and so runs faster during the pulse. First-order
uptake therefore depletes a little more than his, which favours depletion
as a terminator.

**Calibration: Stern's Fig. 20** (his constants and his 28-channel couplon,
80 couplons; Stern's values digitised on the tick positions, in brackets):

| step | content at 100 ms | corrected plateau | C open after the pulse |
|---|---|---|---|
| 0 mV | 0.54 mM [0.655] | 1.11 pA [1.2] | 0.000 |
| −30 mV | 1.19 mM [1.34] | 0.36 pA [0.4] | 0.002 |

Without depletion the same couplons give 2.88 / 1.17 pA peak / plateau at
0 mV, against his corrected 2.9 / 1.2. The model is Stern's, and it
over-depletes by about 0.15 mM, as expected from the uptake.

**Finding 1: depletion restores control, but by emptying the store**
(`python -m ip3r ec --depletion`, 28-channel couplon, 100 couplons, 2 mM):

| C channels | released in 100 ms at +20 / 0 / −30 mV | C Po after repolarisation (0 mV) | flux peak/plateau, corrected (0 mV) |
|---|---|---|---|
| Stern 1997 | 0.80 / 0.73 / 0.40 | 0.000 | 2.63 |
| fitted, no Mg²⁺ | 0.89 / 0.87 / 0.75 | 0.001 (0.72 without depletion) | 2.90 |
| fitted, Mg²⁺ activation site | 0.82 / 0.78 / 0.47 | 0.000 (0.39) | 3.47 |
| fitted, Mg²⁺ both sites | 0.69 / 0.59 / 0.11 | 0.000 | 1.33 |

Taken alone, this looks like a rescue. The fitted C array stops at
repolarisation, and its release still peaks after Schneider's correction:
the collapse of CICR as the currents shrink is a real fall in
permeability, not just a smaller driving force (Stern predicted this for
low loads). But the fitted scheme gets there by releasing 87–89 % of the
store in one pulse.

**Finding 2: at the pool the measurement allows, the peak is gone.** Even
Stern's constants release more than Ríos measured at 2 mM (80 % at +20 mV).
So the pool was scaled until they release the measured fraction
(`ec --depletion --pool-scan`, +20 mV):

| pool | Stern: released, C Po peak/plateau | fitted, no Mg²⁺ | fitted, Mg²⁺ activation site |
|---|---|---|---|
| ×1 (2 mM) | 0.80, 5.7 | 0.89, 4.0 | 0.82, 8.8 |
| ×1.5 | 0.67, 4.7 | 0.86, 2.2 | 0.77, 4.0 |
| **×2 (4 mM)** | **0.57, 4.6** | **0.83, 1.46** | **0.73, 2.7** |
| ×3 | 0.44, 4.6 | 0.74, 1.11 (0.016 open after) | 0.64, 1.65 |

Stern's peak does not depend on the pool (C Po ratio 4.6 from ×1.5 up): it
is inactivation. The fitted peaks exist only while the store is emptied
beyond the measurement, and they shrink as the release approaches it. At
the pool where Stern's constants release Ríos's 57 %, the fitted scheme
still releases 73–83 %, with a C Po ratio of 1.5–2.7 against Stern's 4.6.
Where the fitted scheme would release the measured fraction (beyond ×3),
there is no peak and control is slipping. Local depletion is smaller
still: a skeletal spark lowers free SR Ca²⁺ by at most 7.4 % (Launikonis
et al. 2006), against 54 % for cardiac blinks.

**What it settles.** SR depletion is not the missing terminator. It can
make the fitted scheme peak and stop only by emptying more store than a
fibre loses. This agrees with the experimental conclusion (Launikonis et
al. 2006; Launikonis, Zhou, Santiago, Brum & Ríos 2006, J Gen Physiol 128:45) that
depletion does not terminate release in frog skeletal muscle. What remains
is Round 6.4's alternative: a Ca²⁺ inactivation that is stronger in the
couplon than the steady-state bell shows, a two-site inactivation (Hill
1.5), or inactivation kinetics that bilayer bells do not measure.

**Limits.** The 60-channel couplon of Round 6.6 was not run through the
depletion panel. At a fixed couplon density it would double the channel
density and the release, making over-depletion worse. Ríos's 50–60 % was
measured at 10 °C in 15 mM EGTA, with a content of 1.2–1.5 mM, against
Stern's 2 mM. The pool scan absorbs the difference, but only as one factor.


## Two-site Ca²⁺ inactivation (Round 6.8; `physics/ryr_two_site.py`)

**The question.** With depletion ruled out, the remaining candidate was
inactivation shaped differently from Stern's one-Ca²⁺ gate. Murayama's bell
falls with nI 1.5, where a one-site gate gives slope 1. Note that the 1.5
was *fixed* in their fit, not measured.

**The gate.** Two Ca²⁺ bind in sequence (Stern's on rate for each step),
and the channel is inactivated only with both bound. The uninactivated
fraction is `(1 + x/K1) / (1 + x/K1 + x²/(K1 K2))`. Its Hill slope at the
half point is `2 − c_h/(K1 + c_h)`, which runs from 1 (K2 ≫ K1) to 2
(K2 ≪ K1). The activation gate is Stern's, so the scheme has six states.
Mg²⁺ joins Ca²⁺ at both binding steps, as in the one-site gate.

**Fitted** to three numbers of Murayama's 25 °C bell, measured the same way
on both: the two half-peak points and the log-log slope at half inhibition
(−0.84; the one-site fit gives −0.56). The result is Ka 4.57 µM, K1 512 µM
and K2 115 µM. The gate's own Hill slope is 1.63, and the bell's peak moves
to 34 µM (one-site 23, Murayama 44). All three targets are met to 10⁻⁶.

**What a steeper gate must do.** With the half point fixed near 320 µM, a
steeper flank means *less* inactivation below it. At 30 and 50 µM, the Ca²⁺
a C channel sees in the cleft, the two-site gate inactivates 1.4 % and
3.7 % of channels, against 11 % and 17 % for the one-site fit (tested).

**Finding** (`python -m ip3r ec --two-site`, 200 couplons; one-site fit of
Round 6.6 in brackets):

| C channels | C Po peak / plateau at 0 mV | C Po after repolarisation | flux peak/plateau |
|---|---|---|---|
| two-site, no Mg²⁺ | 0.94 / 0.88 | 0.885 [0.72] | 1.02 |
| two-site, Mg²⁺ activation site | 0.85 / 0.82 | 0.726 [0.39] | 1.02 |
| two-site, Mg²⁺ both sites | 0.091 / 0.082 [0.12] | 0.000 | 1.05 |

At −30 mV, with both sites, the C channels reach 0.016 [0.014]. At −50 mV the
events are single-channel blips: a median of 2 ms, one channel open, none
large. With no Mg²⁺, 199 of 203 events never end.
Triggered cleft sparks (`spark-mg --two-site`) tell the same story. With no
Mg²⁺, none of 20 shuts. With both sites, 26 of 30 channels are inactivated
at rest, and the four that open close in 3.5 ms.

**What it settles.** Matching the bell's slope makes the couplon worse on
every count. Stern's scheme works because its Ki is 10 µM, inside the
cleft's range. Any gate fitted to Murayama's half inhibition near 320 µM
leaves channels in the cleft almost uninactivated, whatever its slope, and
the steeper the slope, the less it inactivates. **No steady-state Ca²⁺
inactivation consistent with the bell is the couplon's terminator.** What
remains are things the bell does not see: the bell's conditions against
the fibre's ([³H]ryanodine binding as a Po index; ATP, luminal Ca²⁺,
temperature, where the 37 °C Ki of 358 µM is further away still), or
inactivation that an equilibrium bell cannot show, such as a gate whose
cycle breaks detailed balance because the flux itself drives it.

## Use-dependent inactivation (Round 6.9; `physics/ryr_use.py`)

The thing the bell does not see, made concrete. Laver & Lamb 1998 stepped
RyR1 in bilayers in voltage and in cytosolic Ca²⁺ and found an inactivation
that takes hold with τ ≈ 1–3 s, affects half to two-thirds of channels, and
proceeds at a rate set by how much the channel is *open*: it "depended on
P(OL) and not on the particular activator (Ca²⁺ (microM), ATP, caffeine,
and ryanodine), inhibitor (mM Ca²⁺ and Mg²⁺), or gating mode".

**The gate.** Stern's two gates, unchanged, plus a third that is entered
**only from the conducting state**, at a Ca²⁺-independent rate. State
`s = a + 2i + 4u`, so eight states, each gate's exit flipping its own bit.
Around the cycle C → O → OU → CU → C the forward rates are
`a k_use a₋ k_use₋` and the reverse are zero, because there is no way into
the gate except through conduction: the scheme has a net stationary cycle
flux and no detailed balance (tested, and calibrated against a variant that
may also be entered when shut, which falls into detailed balance). The
stationary open probability is therefore *not* a product of the three gates'
equilibria and is solved from the null space.

**The two inactivations partly cancel** — the result everything else
follows from. The use gate bites only while the channel conducts, so
shutting it by the Ca²⁺ route *protects* it. The share of channels the use
gate holds down rises with Ca²⁺ to the bell's peak and then falls away
(tested). So the gate does not hide inside a bell: it crushes the peak and
pushes the inhibitory flank out, widening Murayama's 1.86 decades to 2.90 at
ρ = k_use₋/k_use = 0.2 (Round 6.9's τ 2 s against recovery 0.1 s⁻¹).

**So the bell must not be counted twice.** If the channel has both gates,
Murayama's bell is already their composite, and fitting the Ca²⁺ gate alone
and then bolting a use gate on is wrong. `fit_with_use` fits Ka and Ki
*beside* the use gate, and the question becomes whether any Ca²⁺ gate
reproduces the measured bell at a given use-gate speed
(`spark-termination --scan use --bell`):

| use τ | ceiling on P_open | half-peak points | width | Ca²⁺ gate that fits |
|---|---|---|---|---|
| ≤ 0.1 s | ≤ 0.010 | no falling flank < 10 mM | — | **none** |
| 0.32 s | 0.031 | 0.87 – 8197 µM | 3.98 dec | **none** |
| 1 s | 0.091 | 1.47 – 2811 µM | 3.28 dec | Ka 15.5, Ki 21.4 µM |
| 3.2 s | 0.240 | 2.35 – 1108 µM | 2.67 dec | Ka 9.9, Ki 58.7 µM |
| 10 s | 0.500 | 3.29 – 569 µM | 2.24 dec | Ka 6.9, Ki 123.8 µM |

*Corrected in Round 6.10.* The two "none" rows were a stalled solver, not
the bell: with several starting points a Ca²⁺ gate fits beside every gate
down to a recovery ratio of 0.002. The claim that stood here — "the fastest
use gate the bell permits is τ ≈ 1 s, the fast end of Laver & Lamb's
range, an independent agreement" — is withdrawn. At the fixed recovery of
0.1 s⁻¹ used, τ and the recovery ratio were the same scan, and the
agreement was an artefact of that unsourced rate.

**Where the 30× discrepancy came from.** Crediting part of the descending
limb to the use gate lowers the Ca²⁺ gate's Ki from 249 µM to **21–59 µM**
across the τ scanned — within 2–6× of Stern's 10 µM, not 30×. The
discrepancy Rounds 6.3–6.8 chased is what you get by attributing all of a
measured bell's inhibition to Ca²⁺ binding.

**Cleft sparks** (`spark-termination --scan use`, 4 seeds × 10 s):

| C channels | Ka | Ki | sparks | median | unended | open |
|---|---|---|---|---|---|---|
| Stern 1997 | 7.1 | 10.0 | 0.95/s | 19 ms | 0 | 0.007 |
| fitted, Ca²⁺ gate only | 4.9 | 249.1 | 0.10/s | 9796 ms | 4 of 4 | 0.701 |
| + use τ 1 s | 15.5 | 21.4 | 0.05/s | 9 ms | 0 | 0.000 |
| + use τ 3.2 s | 9.9 | 58.7 | 0.90/s | 23 ms | 0 | 0.007 |
| + use τ 10 s | 6.9 | 123.8 | 0.80/s | 658 ms | 2 | 0.256 |
| no inactivation | 4.9 | ∞ | 0.10/s | 9709 ms | 4 of 4 | 0.942 |

**The couplon** (`ec --use`, 200 couplons):

| C channels | Po peak / plateau, 0 mV | after repolarisation | flux peak/plateau | −50 mV events |
|---|---|---|---|---|
| Stern 1997 | 0.644 / 0.117 | 0.0040 | 2.71 | 20 ms, 0 unended |
| fitted, Ca²⁺ gate only | 0.763 / 0.727 | 0.7265 | 1.01 | 185 ms, 199 of 199 unended |
| + use τ 1 s | 0.445 / 0.154 | 0.0002 | 1.80 | 8 ms, 0 unended |
| + use τ 3.2 s | 0.615 / 0.369 | 0.0094 | 1.34 | 22 ms, 0 unended |
| + use τ 10 s | 0.715 / 0.566 | 0.4731 | 1.11 | 184 ms, 193 of 201 unended |

For the first time in Rounds 6.3–6.9, a C-channel scheme that reproduces
Murayama's measured bell also terminates cleft sparks, stops releasing at
repolarisation, and shows a release peak. It is still not the fibre: the
measured release is 6.3 ms against 9–23 ms, and the flux peak/plateau is
1.3–1.8 against Stern's 2.71 and the measured ~4.6.

**What it rests on, and this is the finding.** The use gate's *recovery*
rate could not be sourced — Laver & Lamb report recovery only
qualitatively, and the full text is a paywalled page scan — and the result
depends on it strongly, because recovery sets how much of the bell's
descending limb the gate accounts for, and so where the refitted Ki lands
(`spark-termination --scan recovery`, τ_use 3.2 s):

| recovery τ | ceiling | Ka | Ki | cleft sparks |
|---|---|---|---|---|
| 33 s | 0.087 | 15.8 | 20.3 µM | **none start** |
| 10 s | 0.240 | 9.9 | 58.7 µM | 0.82/s, 23 ms, 0 unended |
| 3.3 s | 0.487 | 7.0 | 120.5 µM | 1.10/s, 266 ms, 1 unended |
| 1 s | 0.760 | 5.6 | 189.0 µM | 0.22/s, 3629 ms, 4 unended |
| 0.33 s | 0.905 | 5.2 | 225.3 µM | 0.10/s, 9506 ms, 4 unended |

Faster recovery leaves the use gate less of the bell to explain, Ki climbs
back towards 249 µM, and termination fails as it did in Round 6.8. Slower
recovery holds the ceiling so low that the cluster never fires. The
mechanism works in a band about a decade wide around a recovery τ of 10 s.

So Round 6.9 does not settle what terminates the couplon. It produces a
**mechanism with a falsifiable requirement**: use-dependent inactivation
is the couplon's terminator only if RyR1 recovers from it with τ of roughly
10 s (within a factor of ~2), and the same scheme must then have
τ_inactivation of 1–3 s, which is measured. A recovery rate read from
Laver & Lamb's full text, or measured afresh, decides it.

### What Laver & Lamb's full text says (Round 6.10)

The full paper settles *what* the model needed from it, but it does not
give the number.

**No recovery rate at a fixed potential is reported.** Inactivated channels
recovered only when the bilayer voltage was reversed ("reactivated in
milliseconds by voltage steps to negative potentials"); at +40 mV a channel
"remained closed until it was reactivated by a brief voltage pulse to
−40 mV" (Fig. 10), and further Ca²⁺ steps never reopened it.

**The 1–3 s is at +40 mV, not at the SR's ~0 mV.** Fig. 4 plots log₁₀ of
the inactivation rate against voltage (cardiac RyRs, Po > 0.8); the
positive limb is 20 V⁻¹ with intercept 0.056 s⁻¹ at 0 mV. The reading is
checked against the paper twice (tested): the line gives τ 2.8 s at +40 mV
(the abstract's 1–3 s) and z δ 1.18 (the text's 1.14 ± 0.25); a
natural-log reading fails both. So at 0 mV the entry rate is ~0.056 s⁻¹
(τ ≈ 18 s), an extrapolation below data that start near +20 mV, and
`ryr.k_use_on` is now that. Skeletal RyRs were not measured against voltage.

**The result rests on a ratio, not a rate.** Both use-gate rates are
seconds; a spark is milliseconds and the Ca²⁺ gate sub-millisecond. The
stationary bell, the refitted Ca²⁺ gate and spark termination therefore
depend on ρ = k_use₋/k_use alone (tested: a 7× slower gate at the same ρ
fits identically; sparks at ρ 0.32 last 19–20 ms at entry rates from 0.3
to 0.045 s⁻¹). `ryr.k_use_off` is replaced by `ryr.use_recovery_ratio`
(unverified, 0.2), and the scans run over ρ
(`spark-termination --scan ratio`, entry rate 0.056 s⁻¹, 4 seeds × 10 s):

| ρ | ceiling | Ka | Ki | cleft sparks |
|---|---|---|---|---|
| 0.02 – 0.063 | 0.02 – 0.06 | 27 – 18 | 4 – 14 µM | **none start** |
| 0.11 | 0.10 | 14.7 | 23.9 µM | 0.03/s, 9 ms |
| 0.20 | 0.17 | 11.7 | 40.3 µM | 0.35/s, 12.5 ms, 0 unended |
| 0.36 | 0.26 | 9.5 | 64.2 µM | 1.10/s, 26.5 ms, 0 unended |
| 0.63 | 0.39 | 7.8 | 95.6 µM | 1.38/s, 88 ms, 0 unended |
| 1.1 | 0.53 | 6.7 | 131.1 µM | 0.68/s, 452 ms, 1 unended |
| 2 | 0.67 | 6.0 | 165.6 µM | 0.25/s, 1160 ms, 4 unended |

Below ρ ≈ 0.1 the refitted Ka is too high for the cluster to fire; above
≈ 0.6 Ki climbs back towards 249 µM and sparks stop ending. Sparks of the
fibre's order (≤ 30 ms) need **ρ ≈ 0.1 – 0.4**. The couplon at ρ 0.2
(`ec --use`, every speed from τ 1 ms to 10 s) keeps control: C after
repolarisation 0.0004–0.0029 against the Ca²⁺-gate-only fit's 0.7265, −50 mV
events 3–14 ms, none unended, flux peak/plateau 1.1–2.4 (1.5 for τ ≥ 0.3 s,
where the speed no longer matters).

**What the paper does bound.** Fig. 8 gives each channel's ensemble
activity 5 s after a step to +40 mV relative to its peak. A residual R
bounds the ratio: a channel held at Po keeps ρ/(ρ + Po) available, so
ρ ≤ R/(1 − R). The skeletal RyRs that inactivated had R 0.03–0.34 (three of
six did not inactivate: the heterogeneity), so **at +40 mV ρ ≤ 0.03–0.52
across channels.** The band the model needs, 0.1–0.4, lies inside that
span, above its most inactivating channels. But the bound is at +40 mV,
where inactivation is voltage-driven, and recovery is steeply voltage
dependent the other way (milliseconds at −40 mV): ρ at 0 mV can lie above
the +40 mV bound. So Round 6.9's requirement is restated, not met: **the
use gate terminates the couplon if RyR1 at 0 mV keeps ρ ≈ 0.1–0.4 — a
steady residual activity of roughly 10–30 % in channels held open**, a
number a bilayer at 0 mV, or a longer record at +40 mV than Fig. 8's 5 s,
can measure.

**Heterogeneity, modelled in Round 6.11 (`physics/ryr_mixed.py`;
`spark-termination --scan fraction`).** Only some channels have the use
gate at all, stably: "only those that showed inactivation after voltage
steps inactivated after [Ca²⁺] steps, and vice versa". Laver & Lamb count
80 % of skeletal RyRs (12 of 15, Po > 0.2), 56 % of cardiac; Laver & Curtis
1996 70 % of 25 after Ca²⁺ steps; Ma 1995 50–70 %. `ryr.use_inactivating_fraction`
is 0.8. In the mixed cluster a random `round(f·n)` channels carry the gate
and all share one Ca²⁺ gate. Because Murayama's bell is [³H]ryanodine
binding by a whole preparation, **the Ca²⁺ gate is fitted to the
population's mean bell**, f·Po_use + (1 − f)·Po_plain.

At ρ 0.2: fraction 1 gives Ka 11.7, Ki 40 µM and 14 ms sparks; **0.8 gives
Ki 104 µM and 228 ms sparks (2 of 53 unended)**; 0.6 or less, sparks never
end (Ki 167–249 µM). At ρ 0.32 the same shape (22 ms → 326 ms). Sparks of
the measured order need f ≥ 0.9 (27 ms at 0.9, 19 ms at 0.95). The couplon
at 0.8 loses control: C Po after repolarisation 0.33 (0.0018 at f = 1),
and at −50 mV 148 of 216 large events unended.

**The damage comes through the bell, not through the channels.** Taken
apart at f = 0.8: the all-use Ca²⁺ gate with only 80 % carrying the use
gate gives 14 ms sparks, as at f = 1: the non-inactivating fifth is,
plausibly, shut by its neighbours: their inactivation removes its
Ca²⁺. The 0.8-fit gate with *every* channel carrying the use gate gives
151 ms. What lengthens the spark is that a population bell shared with
channels that never use-inactivate leaves the use gate less of the
descending limb, so the Ca²⁺ gate must be weaker. Round 6.9 predicted the
direction ("can only make termination harder"), not the mechanism.

**Standing after 6.11.** The use gate reconciles the bell with
termination only if nearly every channel carries it (f ≥ 0.9), against a
measured 0.5–0.8. Two assumptions carry this and both are open. (i) The
two populations share one Ca²⁺ gate; if the non-inactivating channels had
a stronger Ca²⁺ inactivation of their own, the bell could be split
differently. (ii) The bell's preparation has the bilayer's fraction.
Neither is sourced. Laver & Curtis 1996 is at +40 mV throughout. **ρ at
0 mV is still unmeasured.**

**Sitsapesan et al. 1995, in full (read after 6.11).** Sheep cardiac RyRs,
at −40 and +40 mV only, so no ρ at 0 mV. Three things it does fix:
- *+40 mV:* channels that inactivated stayed shut through repeated agonist
  steps until a brief −40 mV pulse, a residual near zero. This matches
  Fig. 8's most inactivated channels (ρ ≤ 0.03).
- *−40 mV:* held at Po 0.84–0.93 for 5 s, no decline (Po 0.841 → 0.874,
  n = 5). A use gate can lower activity by at most 1 − exp(−k·Po·t)
  whatever its ρ (`decline_bound`). The registered 0 mV rate, carried to
  −40 mV by Laver & Lamb's own slope, predicts 3.7 %, within one SD.
  **Round 6.9's 0.5 s⁻¹ predicts 29 %, beyond two SDs**, so the record
  independently rejects the reading Round 6.10 corrected. It does so
  narrowly: the largest rate at 0 mV it allows is about 0.44 s⁻¹.
- *Fraction:* 17 % (4 of 24) inactivated with Ca²⁺ alone, 56 % (5 of 9)
  with ATP or EMD 41000. This is the Po dependence that Laver & Lamb's
  Po > 0.2 selection removes, not a different population.

**The induction-decay control.** With *neither* inactivation gate, cleft
sparks never end (open fraction 0.942, 4 of 4 unended at 9.7 s). This model
therefore cannot exhibit Laver et al. 2013's "induction decay", in which
the cleft's geometry terminates release on its own, and the reason is
structural: `sparks_cleft` holds the cleft field at its steady state
between gating events, while induction decay lives in the time course of
the nanoscopic gradients — "the closed RyRs do not respond in the timescale
over which the very local [Ca²⁺] is maintained". Testing it needs a
time-dependent cleft, which is recorded as emergent.

## Copello 1997's low-activity channels (Round 6.12; `ryr_mixed.low_activity`)

Round 6.11 left two unsourced assumptions, the first being that both
populations share one Ca²⁺ gate. Copello et al. 1997 is the paper on RyR
heterogeneity of Ca²⁺ gating, and it reports a second population rather
than a second gate: **about 35 % of skeletal RyRs (7 of 20 channels, 4 of
14 with Mg²⁺/ATP) gate in a "low-activity" mode**, Po ≤ 0.1 at every Ca²⁺,
activated at 70–150 µM and inhibited at 100–300 µM — a low, narrow bump
sitting under the population bell's descending half-point. (Their
high-activity channels vary widely too: EC₅₀ 0.7–8.9 µM, IC₅₀ 0.16–1.1 mM,
Hill slopes 0.8–5.2, which is why this project's single fitted gate is a
population average, not a channel.)

`low_activity(c, reading)` is that bump as a Hill pair, with the ceiling
and both half-point ranges registered and `reading` placing the half
points at the low ends, the geometric middles or the high ends of
Copello's ranges (`LA_READINGS`; the high reading sits nearest Murayama's
half-inhibition and so does the most). `fit_mixed(...,
low_activity_reading=)` puts `ryr.la_fraction` of the population there and
fits the high-activity channels' shared Ca²⁺ gate to the whole. They are
**not** put in the cleft: uninactivating open channels there could only
lengthen a spark, so every row below is a best case.

**Measured** (`spark-termination --scan low-activity`, f = 0.8, ρ 0.2):

| Population bell | Ki | Spark |
|---|---|---|
| high-activity only (6.11) | 104 µM | 228 ms |
| + low-activity, low reading | 106 µM | 157 ms |
| + low-activity, mid reading | 95 µM | 147 ms |
| + low-activity, high reading | 65 µM | 48 ms |

The direction is right and the size is not. Channels that are open on the
descending limb without inactivating force the high-activity gate to
inactivate harder, so Ki falls and sparks shorten; but even the most
favourable placement of Copello's numbers gives 48 ms against a measured
6.3 ms, and that is with the low-activity channels kept out of the cleft.
**Copello's heterogeneity does not rescue the use gate.** The first of
Round 6.11's two assumptions is therefore addressed and survives: the
population's spread is real and documented, and it is not large enough.
What is still unsourced is a *different Ca²⁺ inactivation gate* in the
non-use-inactivating channels, which Copello does not report — its
heterogeneity is in activation and in overall activity, and its IC₅₀
spread (0.16–1.1 mM skeletal) straddles the bell's Ki rather than
reaching Stern's 10 µM.
