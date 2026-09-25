"""The thresholds' uncertainty (Paper 5 §8): the exact interval for a median
by hand, its coverage, why not a percentile bootstrap, and which VUS a
threshold inside the intervals could move."""

import math

import numpy as np
import pytest

from conftest import needs_genes
from ip3r.analysis.stats import median_interval
from ip3r.analysis.vus_strata import Bands


def test_order_statistic_ranks_by_hand():
    # Conover's table: n = 10 → x(2), x(9); n = 50 → x(18), x(33) at 95 %
    assert median_interval(np.arange(10)) == (1.0, 8.0)
    assert median_interval(np.arange(50)) == (17.0, 32.0)
    assert median_interval(np.arange(6)) == (0.0, 5.0)     # 2/64 = 0.031 ≤ 0.05


@pytest.mark.parametrize("n", [1, 2, 5])
def test_too_few_is_unbounded(n):
    assert median_interval(np.arange(n)) == (-math.inf, math.inf)


def test_coverage_on_a_skewed_distribution():
    rng = np.random.default_rng(0)
    hits = sum(lo <= math.log(2) <= hi for lo, hi in
               (median_interval(rng.exponential(size=17)) for _ in range(3000)))
    assert hits / 3000 >= 0.94


def test_a_percentile_bootstrap_would_call_one_point_certain():
    """The reason this is not a bootstrap: resampling ITPR2's single P/LP
    position gives a zero-width interval, the most certain-looking threshold
    for the least-known median."""
    x = np.array([0.7])
    boot = np.median(np.random.default_rng(0).choice(x, (1000, 1)), axis=1)
    assert np.ptp(np.percentile(boot, [2.5, 97.5])) == 0.0
    assert median_interval(x) == (-math.inf, math.inf)


def test_unsure_by_hand():
    b = Bands((0.7, 0.8), (0.3, 0.4), 0.95)
    assert b.unsure(0.85) == () and b.unsure(0.2) == () and b.unsure(0.5) == ()
    assert b.unsure(0.75) == ("P/LP",) and b.unsure(0.35) == ("B/LB",)
    assert b.unsure(0.8) == () and b.unsure(0.7) == ("P/LP",)   # edges: ≥ counts
    assert b.unsure(float("nan")) == ()
    u = Bands((-math.inf, math.inf), (0.3, 0.4), 0.95)
    assert u.unsure(0.99) == ("P/LP",) and u.unsure(0.1) == ("P/LP",)


@needs_genes
def test_real_thresholds():
    from ip3r.analysis.vus_strata import bands
    from ip3r.render.variant_spheres import resource_stratification
    s1, s2, s3 = (resource_stratification(g, "deep") for g in ("ITPR1", "ITPR2", "ITPR3"))
    b1, b2, b3 = bands(s1, 0.95), bands(s2, 0.95), bands(s3, 0.95)
    assert b1.bounded and not b2.bounded and not b3.bounded   # 39 vs 1 vs 5 positions
    assert len(s2.pathogenic) == 1 and len(s3.pathogenic) == 5
    assert b1.pathogenic[0] <= s1.median_pathogenic <= b1.pathogenic[1]
    assert len(b2.near(s2, "P/LP")) == s2.n_vus               # no VUS is firm about it
    assert len(b1.near(s1, "P/LP")) < s1.n_vus
