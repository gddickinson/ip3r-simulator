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
