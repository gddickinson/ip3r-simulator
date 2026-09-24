# SCIENCE_RYR — the ryanodine receptor models

The RyR1 half of `SCIENCE.md`: structures, conductance and charge mutants,
gating, and sparks (mean-field and in the junctional cleft).

## Ryanodine receptor 1

**Why.** RyR1 is IP3R's closest relative (both share the trefoil, MIR, RIH,
RIH-associated and PF00520 domains), with many more deposits, a larger
conductance and measured pore mutants. It is the natural control for
anything measured on the ITPR3 pore.

**Resource** (`scripts/curate_ryr.py` → `resources/ryr1.json`, SHA-256 of
every response). Reference rabbit P11716 (5,037 residues), because the
deposits are rabbit; the GGGIGD filter is at 4894–4899 and the gate
isoleucine at 4937, as in the literature. Domains are InterPro's Pfam
matches; PF00520 (4789–4947) is named `channel`, so the pore-domain span is
found the same way as for the ITPRs. RyR1 has no conservation, sites or
variants here, so those colourings are grey on it.

**Panel, by rule**, from the 158 PDB entries mapped to P11716 (on 2026-09-24):
single-particle EM of four full-length chains (not domain crystals or local
refinements); wild type; only Ca²⁺/Zn²⁺/Mg²⁺, ATP/ACP/ADP, caffeine or
lipid as ligands; ≤ 4.0 Å; the state the title names, best resolution per
state. That gives 9OL6 closed (3.11 Å), 8RRX primed (3.10), 9HEO open
(3.40), 7TDG inactivated (3.80) and 7TDI closed-inactivated (3.30). The morph
needs one preparation, so the best primed deposit of the open state's own
paper, 9R8O (3.30, same DOI as 9HEO), is added. 140 entries fail a rule,
each recorded in the resource with the first rule it fails; 18 pass, and
the best per state is taken from those.

**Measured.** Every deposit is 100 % in P11716 numbering. The four shut
states gate at I4937 (2.6–3.3 Å r_min) and 7TDI at Q4933 (2.21 Å). Only
9HEO opens (5.05 Å, lined by Q4933). Primed 9R8O → open 9HEO morphs over
4,235 residues × 4 subunits (RMSD 2.72 Å after the pore fit). The start's
elastic network explains little of it: the cumulative overlap of 20 modes is
0.174 against a null of 0.031, where ITPR3's is 0.64. Most of the raw
displacement (71.5 %) is rigid-body, removed before the overlap.

**Conductance.** Symmetric 250 mM KCl, the bath of Xu et al. 2006's
recombinant RyR1 (801 ± 7 pS, n = 17, planar bilayer). 9HEO: 134 pS series,
136 pS neutral, 180 pS charged, 161 pS paired; the sweep reaches at most
308 (neutral) and 449 pS (charged). The neutral model is 5.9× short. So
ITPR3's 2.4× shortfall is not particular to ITPR3. The same continuum
under-predicts both receptors, and RyR1's acidic wall (−32 e) *raises* its
conductance where ITPR3's alternating rings lowered it. Primed 9R8O
(r_free 1.62 Å) conducts 4.5 pS neutral. Its charged solve does not
converge (Cl⁻ is excluded, and K⁺ alone cannot neutralise a basic zone), and
it is reported as n.c., never as 0.

**Charge mutants** (`physics/ryr_mutants.py`, `python -m ip3r mutants`). Each
mutant drops one residue's charge on all four subunits of 9HEO. The ratio
mutant / wild type cancels the transport constants, so it tests the wall
charge alone:

| mutant | measured | model, charged | model, paired | in 9HEO |
|---|---|---|---|---|
| D4899Q | 0.20 | 0.90 | 1.00 | lining, bridged to R4892′ |
| E4900N | 0.63 | 0.76 | 0.41 | lining |
| D4938N | 0.65 | 0.78 | 0.80 | lining |
| D4945N | 0.92 | 0.90 | 0.91 | lining |
| E4955Q | 1.01 | 1.00 | 1.00 | not lining |

The direction is right for all four lining residues, the size right for
D4945N, and the null right for E4955Q. The filter aspartate is 4× too weak
in the charged model and absent in the paired one. A continuum Donnan
partition spread over 3 Å does not carry what a ring of four carboxylates
does at a 5 Å filter.

## RyR1 gating and sparks

**Sources, and why these.** A search for RyR1 models whose constants can be
read from the source's own text (many papers are behind a script
challenge) found one measured bell and one kinetic scheme.

- *Murayama et al. 2015* (PLoS One, S1 Table): recombinant rabbit RyR1,
  Ca²⁺-dependent [³H]ryanodine binding, A = Amax fA (1 − fI) with
  KA 5.5 µM, nA 1.2, KI 0.27 mM, nI 1.5 (wild type, 25 °C). Binding is an
  activity index, not P_open, so it is compared by shape (flanks), never by
  height.
- *Stern, Pizarro & Ríos 1997* (J Gen Physiol, Table I): the skeletal "C
  channel", two gates in series. Activation opens on two Ca²⁺
  (k_o = 10 µM⁻² s⁻¹, k_o− = 500 s⁻¹) and inactivation closes on one
  (k_i, k_i− = 20 s⁻¹). k_i is printed "2 × 10⁻⁶ M⁻¹ s⁻¹", a sign typo: the
  text gives the inactivation Kd as 10 µM, so k_i = 2 µM⁻¹ s⁻¹. The printed
  value would remove inhibition entirely (tested). The authors say the
  model "has not been objectively 'fitted' to data".

The best-constrained readable kinetic scheme, Zahradníková et al. 1999, is
cardiac RyR2 and has no Ca²⁺ inhibition, so it was not used.

**Bells on one ruler** (`python -m ip3r ryr-gating`). Half-activation: 3.9 µM
(scheme) vs 4.4 µM (measured). Half-inhibition: 48 vs 320 µM. The scheme
activates where RyR1 does, but inactivates 6.7× too readily. Stern et al.
said as much of their inactivation site ("one or two orders of magnitude
lower than in bilayers").

**Sparks** (`physics/sparks.py`, `python -m ip3r sparks [--scan]`). The
cluster is the 30 Ca²⁺-gated C channels of Stern's 60-channel couplon, with
the mean-field coupling the IP3R puffs use. The coupling is derived:
free-diffusion Ca²⁺ one channel spacing (30 nm) from one open channel
(0.3 pA, D 5 × 10⁻⁶ cm² s⁻¹), which is 8.25 µM, an unbuffered upper bound.
Each step is exact (expm of the generator, tabulated per number open). The
step was measured: 1e-5 and 2.5e-5 s agree within noise, 1e-4 s lengthens
sparks by 15–20 %, and a 5 ms step is caught.

Read with the puff ruler (10 s):

- Uncoupled: only single-channel blips.
- Coupled: 1.5 sparks per second, most reaching 25–30 of 30 channels,
  with a size gap between blips and sparks; Fano 4.
- Coupling scan (× derived value): sparks switch on between 0.05× and 0.1×
  (0.44–0.80 µM per open channel) and saturate by ~0.2×. So buffering
  could cut the coupling five-fold without changing the answer.

**Where it fails.** Sparks last ~120 ms (median), against a measured
release of 6.3 ms (Ríos et al. 1999, frog). The trace shows why. After the
first near-whole-cluster peak, the mean-field cluster settles at the point
where 30 · P_open(0.1 + 8.25 n µM) = n (between 5 and 6 open). It stays
there until the number open falls to zero by chance. Stern et al.'s
termination relied on local geometry and voltage-sensor coupling, which a
single cluster Ca²⁺ cannot represent. Measured spark durations need a
spatial Ca²⁺ field, as they had, not a better constant.

(Measured later, Round 6.3: at the end of a mean-field spark 27 of 30
channels are inactivated, so the self-sustaining point sits on a mostly
inactivated cluster, fed by the few channels that recover.)

## Sparks in the junctional cleft (Round 6.3)

**The geometry is Stern et al.'s** (Table I, Fig. 7 B, Appendix;
`physics/cleft.py`). The cleft is 60 nm wide, 15 nm high and as long as the
couplon, with two rows of sites 30 nm apart. V and C channels alternate
like a chessboard, and only the 30 C channels are simulated. An open
channel releases 0.3 pA uniformly over a 30 nm disc (the foot). Ca²⁺ leaves
across every edge at `D_inf c`, with `D_inf = 2πD / (h ln(R_max/h))`
(their Eq. 13, R_max = 1 µm). The steady field is solved by finite volumes
on a 1 nm grid. It conserves the release exactly, and 0.5 nm changes no
coupling by 1 %. Stern et al. show the cleft settles in microseconds
(Fig. 10), so each channel's Ca²⁺ is the superposition of the open
channels' steady fields: `c_i = c_rest + Σ_j G_ij open_j`.

**Against their Fig. 9** (read from the figure, per pA, ±5 µM per pixel):

| Position | Stern | This solve |
|---|---|---|
| Source centre | ~185 µM | 171 µM |
| Across the row | ~75 µM | 65 µM |
| 30 nm along | ~73 µM | 56 µM |
| 60 nm along | ~25 µM | 15 µM |

The shapes agree, but this solve falls off faster. At 0.3 pA a channel sees
51 µM from its own release and 11.4 µM from each diagonal C neighbour
(42 nm). With every other channel open it sees 36 µM, against 239 µM in the
mean-field cluster.

**Exact simulation** (`physics/sparks_cleft.py`,
`python -m ip3r sparks --cleft`). The field is constant between gating
events, so the array is a Markov process. It is simulated by Gillespie's
method, as Stern et al. did, with no step to converge. An uncoupled array
at a clamped 5 µM reproduces the stationary P_open within 5 %. The
"coupling" knob is the nearest-neighbour value: the off-diagonal of `G` is
scaled and a channel's own release is not.

**Measured** (8 seeds × 20 s):

- Sparks last a median 21 ms (IQR 17–26), against ~130 ms in the
  mean-field cluster. They start with 0 inactivated channels and end with
  16 of 30 inactivated: the regenerative front burns out by local
  inactivation, as Stern et al. describe ("terminated by local
  inactivation").
- There are 0.9 sparks per second reaching half the cluster (median peak
  18), with a broad size distribution rather than the mean-field
  cluster's all-or-none. Uncoupled, only blips.
- The coupling is a threshold. At 0.5× the solved coupling no event
  reaches half the cluster; at 0.75–2× sparks last 18–21 ms.
- Self-coupling, the number Stern et al. call least certain: a source
  diameter of 15–45 nm (own release 74–40 µM) gives 16–22 ms. A 4× finer
  observation bin gives 20 ms.

**What it settles, and what it does not.** The ~120 ms mean-field spark was
an artefact of one cluster Ca²⁺. With a spatial field the duration falls
six-fold and no geometric uncertainty moves it far. It is still 3× the
6.3 ms frog release. Since no geometric number closes that gap, the next
suspect is the gating scheme itself: Stern's inactivation, which the
Murayama bell says is 6.7× too sensitive at steady state. Mobile buffers
(their fura-2 runs) need time-dependent diffusion and are not modelled.
Fixed fast buffers leave a steady field unchanged.

## What ends a cleft spark (Round 6.4)

**The suspicion.** Round 6.3 left cleft sparks at ~20 ms, 3× the measured
release, and named Stern's inactivation as the suspect. At steady state
it is 6.7× more sensitive than Murayama's measured bell.

**The refit** (`ryr_gating.fit_to_bell`). The steady state fixes only two
ratios, Ka = √(k_o−/k_o) and Ki = k_i−/k_i. The fit moves the off rates,
keeps Stern's on rates, and solves for both half-peak points exactly. A
one-Ca²⁺ gate has a Hill slope of 1 where Murayama fixed nI at 1.5, so the
slopes are not matched.

| Target | Ka (µM) | Ki (µM) |
|---|---|---|
| Stern 1997 as published | 7.1 | 10 |
| Murayama, 25 °C (KA 5.5 µM, KI 0.27 mM) | 4.9 | 249 |
| Murayama, 37 °C (KA 20.5 µM, KI 0.41 mM; S1 Table) | 16.9 | 358 |

**The cleft array with each** (`physics/spark_termination.py`,
`python -m ip3r spark-termination`; 4 seeds × 10 s). A spark still running
when the trace ends is counted as unended, and its duration is only a lower
bound.

- Stern: 38 sparks, median 19 ms, all ended, 16 of 30 channels
  inactivated at the end.
- Fitted to 25 °C: after the first spark the array never shuts. All 4
  sparks are unended (one per run), and 70 % of channel-time is open.
- Fitted to 37 °C: no sparks. With Ka 17 µM, a neighbour's ~11 µM cannot
  recruit.

**Ki scan** (Stern's Ka and on rate): median duration 19/18/17 ms at Ki
3/5.5/10 µM, then 25, 33 and 89 ms at 19, 35 and 64 µM, then 1.4 s at 118 µM
with most sparks unended. From 217 µM up, every spark is unended.
Termination by inactivation needs a Ki comparable to the Ca²⁺ a channel
sees in the cleft (tens of µM). The measured bell puts it 5–10× higher.

**Rate scan** (Ki fixed; steady state leaves the rate free). At Stern's
Ki, slower inactivation lengthens sparks (162 ms at 1/30×). Faster
inactivation shortens them only to ~13–15 ms, and by then almost none
start (0.1/s at 3×, none at 30×): inactivation wins before recruitment
does. At the fitted Ka and Ki, no rate from 1/30× to 30× ends a spark.

**What this settles.** The 6.3 ms release cannot come from this scheme's
Ca²⁺ inactivation. Values that agree with the measured bell do not end
sparks at all. Values that end them never reach 6 ms while sparks still
fire. Stern's too-sensitive inactivation was not an error that made sparks
long. It was what let them end. Murayama's bell was measured without Mg²⁺
(1 mM AMP, 0.17 M NaCl). In the fibre, ~1 mM free Mg²⁺ both competes at the
activation site and binds the low-affinity inhibitory site, so the bell a
spark sees is not the one measured here. Mg²⁺, luminal Ca²⁺ depletion and
the V channels are the mechanisms left. Each is recorded as emergent in
the roadmap, not assumed.
