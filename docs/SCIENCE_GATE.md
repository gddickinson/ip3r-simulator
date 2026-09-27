# The gate's own geometry (Round 7.21)

`python -m ip3r gate [PDB ...] [--delta D] [--spacing H]`;
`ip3r/physics/gate_geometry.py`; `tests/test_gate_geometry.py`. Continues
`SCIENCE_SEL3D.md` (Round 7.19) and the shortfall in `SCIENCE_PERM.md`
(Round 7.6).

## The question

Round 7.19 left two readings pointing at the gate. First, 9HEO's uncharged
gate held 53 % of Ca²⁺'s resistance under `pb + csc`, and P_Ca:P_K stayed
at 1.08 against Xu's 7.0. Second, the deposits' gates are narrower than a
pore that conducts the measured currents (Round 7.6). The candidate was
therefore: the deposited gate is too narrow. Widen it until it no longer
limits, and see whether conductance and selectivity follow.

## The move

`widen_gate` moves every atom within `gate.widen_half_width` (8 Å) of the
gate constriction radially outward by Δ·cos²(π(z − z_gate)/2w). The shift
is full at the gate and zero at ±w. Every atom at one height moves the same
distance, so the move keeps C4 and the axis, and the gate's charges ride
with their residues. 9HEO's filter lies 10 Å luminal of its gate, outside
the window. The scan runs Δ = 0, 2, 4 Å (`gate.widen_max`,
`gate.widen_step`). Each widened copy is measured afresh. The shares are
read at the *deposited* constrictions, so "the gate" stays the same place
while it opens.

**Calibrated first**:
- on synthetic atoms, each atom moves outward by exactly the taper (1e-4
  Å) with no turn and no axial move, atoms outside the window are
  untouched bit for bit, and a C4 copy stays C4;
- a ring of radius 3 Å at the gate ends at 3 + Δ;
- on 9HEO the measured gate radius grows by Δ (±0.3 Å), the filter moves
  < 0.05 Å and the axis < 0.1°.

## Measured

Grid 0.5 Å; P_Ca:P_K by Round 7.19's linear response; in brackets, the
share of Ca²⁺'s resistance at the deposited gate / filter.

| deposit | +Δ Å | gate r | K⁺ g (pS) | neutral | pb | local + csc | pb + csc |
|---|---|---|---|---|---|---|---|
| 9HEO (RyR1, 7.0) | 0 | 5.05 | 249 | 0.63 | 1.01 (51/5 %) | 0.85 | 1.08 (53/3 %) |
| | 2 | 6.90 | 299 | 0.60 | 0.99 (6/8 %) | 0.98 | 1.18 (6/5 %) |
| | 4 | 8.74 | 308 | 0.60 | 1.01 (7/9 %) | 1.10 | 1.27 (6/6 %) |
| 8TKF (ITPR3, 15.2) | 0 | 5.85 | 98 | 0.54 | 1.03 | 0.88 | 1.39 (31/6 %) |
| | 4 | 9.51 | 109 | 0.57 | 0.97 | 1.01 | 1.50 (16/9 %) |
| 7T3T (ITPR3, 15.2) | 0 | 5.34 | 82 | 0.84 | 1.41 | 1.05 | 1.79 (0/9 %) |
| | 4 | 9.10 | 98 | 0.64 | 1.28 | 0.94 | 1.53 (1/9 %) |

## What it shows

- **No gate radius conducts the measured current.** A 4 Å wider gate lifts
  K⁺ conductance by 24 % (9HEO), 11 % (8TKF) and 20 % (7T3T), and the
  gain has saturated by 2 Å. Once the gate is wider than the filter, the
  resistance is spread over the rest of the ~50 Å pore (Round 7.6), which
  the gate never held. Round 7.6's favourable corner (bulk diffusivity, a
  small K⁺) is still what brings RyR1 to 801 pS. Shape at the gate is not.
- **The gate is not what caps P_Ca:P_K.** In 9HEO the gate's share of
  Ca²⁺'s resistance falls from 53 % to 6 %, yet P_Ca:P_K goes only from 1.08
  to 1.27 against 7.0. In 8TKF it goes from 1.39 to 1.50 against 15.2, and
  in 7T3T it *falls*. Round 7.19's reading, "the gate holds the resistance,
  so the gate caps the ratio", was a correlation. Remove the gate and the
  resistance moves to the uncharged stretches either side. There neither
  ion is enriched, so the ratio returns toward the neutral pore's D_Ca/D_K
  (0.54–0.64).
- **What linear response can give.** P_Ca:P_K here is D_Ca/D_K times Ca²⁺'s
  enrichment over K⁺, averaged harmonically along each ion's resistance. A
  ratio of 7–15 needs Ca²⁺ enriched 12–25× over K⁺ along most of that
  resistance, not in one ring. The deposits place their acidic charge in
  the filter and the luminal vestibule, and that stretch holds ≤ 10 % of
  Ca²⁺'s resistance.

## What remains

The continuum is exhausted at the gate as well: its field (7.19) and its
shape (this round). Left are the two candidates outside the deposit's
wall:
- charges the deposit does not place (unresolved lining residues,
  groups nearer the cytosolic vestibule, protonation);
- the reversal protocol itself. Xu's 7.0 and Vais's 15.2 are bi-ionic
  reversal potentials, where binding in the filter sets E_rev without
  setting the conductance ratio. The 1-D model already reads 0.87 in
  linear response against 0.64 at reversal (Round 7.19), a difference in
  the opposite direction, so a 3-D reversal solve is the test.
  *Tested in Round 7.23 (`SCIENCE_REVERSAL.md`), negative:* in 3-D the
  reversal reads within 15 % of linear response (9HEO 0.90, 8TKF 1.15).
