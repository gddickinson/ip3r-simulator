"""RyR1 gating and sparks: the scheme against itself, the source and the ruler.

* The stationary occupancy from the generator's null space equals the
  two-gates-in-series product. A wrong generator entry breaks it.
* The inactivation constant is the text's 10 µM. The printed k_i
  (2 x 10^-6 M^-1 s^-1) would put Ki at 10^7 µM, and the bell would never
  fall.
* The derived coupling is recomputed by hand from Stern's Table I.
* A clamped, uncoupled cluster reproduces the stationary P_open, so the
  sampling is right. The step converges, and a 5 ms step visibly does not.
* Coupling makes sparks and removing it leaves only blips. That is the
  mechanism, with the ruler IP3R puffs are read with.
"""

import math

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import ryr_gating as rg
from ip3r.physics import sparks as sk
from ip3r.physics.puff_compare import recruitment, spark_scan


@pytest.mark.parametrize("c", [0.01, 0.3, 3.0, 12.0, 100.0, 2000.0])
def test_stationary_is_the_product_of_two_independent_gates(c):
    pi = rg.stationary(c)
    sp = rg.SternParams()
    f_a = c ** 2 / (c ** 2 + sp.k_a ** 2)
    f_i = c / (c + sp.k_i)
    assert pi == pytest.approx([(1 - f_a) * (1 - f_i), f_a * (1 - f_i),
                                (1 - f_a) * f_i, f_a * f_i], rel=1e-9, abs=1e-15)
    assert rg.open_probability(c) == pytest.approx(pi[rg.OPEN], rel=1e-9)


def test_constants_are_the_texts_and_the_printed_typo_would_remove_inhibition():
    sp = rg.SternParams()
    assert sp.k_i == pytest.approx(10.0)            # "the K d for inactivation, also 10 uM"
    assert sp.k_a == pytest.approx(math.sqrt(50.0))  # 7.1 uM (the text rounds to 10)
    typo = rg.SternParams(k_inact_on=2e-6 * 1e-6)    # as printed, in uM^-1 s^-1
    assert np.isnan(rg.bell_at(typo).c_half_inh)
    assert np.isfinite(rg.bell_at().c_half_inh)


def test_the_scheme_activates_like_rabbit_ryr1_but_inactivates_too_readily():
    s, m = rg.compare_bells().values()
    assert s.c_half_act == pytest.approx(m.c_half_act, rel=0.2)
    assert m.c_half_inh / s.c_half_inh > 5


def test_murayama_bell_is_their_equation():
    mp = rg.MurayamaParams()
    c = mp.k_a
    f_i = c ** mp.n_i / (c ** mp.n_i + mp.k_i ** mp.n_i)
    assert rg.murayama_activity(c) == pytest.approx(mp.a_max * 0.5 * (1 - f_i))


def test_derived_coupling_by_hand():
    # 0.3 pA / (2 F) over 4 pi x 5e-10 m^2/s x 30 nm, in uM.
    expected = 0.3e-12 / (2 * 96485.33212) / (4 * math.pi * 5e-10 * 30e-9) * 1e3
    assert sk.diffusion_coupling() == pytest.approx(expected, rel=1e-9)
    assert 8.0 < expected < 8.5


def test_clamped_uncoupled_cluster_reproduces_the_stationary_open_probability():
    pp = sk.SparkParams()
    pp.ca_per_open, pp.ca_rest = 0.0, 5.0
    tr = sk.simulate_sparks(0.0, 20.0, 1, pp)
    frac = tr.n_open.mean() / pp.n_channels
    assert frac == pytest.approx(float(rg.open_probability(5.0)), rel=0.05)


def _rate(dt, seeds=range(3), duration=10.0):
    out = []
    for s in seeds:
        pp = sk.SparkParams()
        pp.dt = dt
        out.append(recruitment(sk.simulate_sparks(0.0, duration, s, pp))["large_per_s"])
    return float(np.mean(out))


def test_the_step_converges_and_a_coarse_step_is_caught():
    base = _rate(_P.value("spark.dt"))
    assert _rate(_P.value("spark.dt") / 2) == pytest.approx(base, rel=0.25)
    assert _rate(5e-3) < 0.5 * base


def test_coupling_makes_sparks_and_its_absence_leaves_blips():
    rows = spark_scan(duration=10.0, seed=0)
    off, on = rows[0], rows[-1]
    assert off["coupling"] == 0.0 and off["max_open"] <= 2 and off["large"] == 0
    assert on["large"] > 5 and on["large_share"] > 0.5 and on["fano"] > 2
    assert rows[1]["large"] == 0          # 3 % of the derived coupling: no sparks
