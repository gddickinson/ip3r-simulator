"""Round 7.25: the series bound B = b a² ≤ 1, its escapes, and the well
ceiling on a tube."""

from types import SimpleNamespace

import numpy as np
import pytest

from ip3r.physics.ohmic3d import cylinder_volume, geometric_conductance
from ip3r.physics.selectivity_bound import (band_mask, bound_product,
                                            linear_ratios, parallel_ratios,
                                            profile_ratios, well_ceiling)
from ip3r.physics.selectivity3d import bath_species


def _product(r):
    a, b = r[2] / r[1], r[-1] / r[1]
    return b * a * a


def test_uniform_potential_by_hand_and_equality():
    for c in (-3.0, -0.5, 0.0, 1.2):
        r = profile_ratios(np.full(50, c), np.full(50, 7.0))
        for z in (1, 2, -1):
            assert r[z] == pytest.approx(np.exp(-z * c), rel=1e-12)
        assert _product(r) == pytest.approx(1.0, rel=1e-12)


def test_no_series_profile_exceeds_the_bound():
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(2000):
        n = int(rng.integers(2, 60))
        psi = rng.normal(0, rng.uniform(0.1, 4.0), n) + rng.uniform(-5, 5)
        area = rng.uniform(1.0, 200.0, n)
        worst = max(worst, _product(profile_ratios(psi, area)))
    assert worst <= 1.0 + 1e-12


def test_bound_product_is_the_ratio_form():
    # 0.27 / 0.29 and 15.2 / 0.26: Vais's pair over the uncharged pore
    assert bound_product(0.27, 15.2, 0.29, 0.26) == pytest.approx(
        (0.27 / 0.29) * (15.2 / 0.26) ** 2)


def test_parallel_paths_of_opposite_sign_escape_it():
    area = np.full(20, 10.0)
    r = parallel_ratios([(np.full(20, -2.0), area), (np.full(20, 2.0), area)])
    assert _product(r) > 10.0            # cosh(4)²/cosh(2)³ × … ≫ 1
    # one path alone, or two alike, stays on the bound
    same = parallel_ratios([(np.full(20, -2.0), area), (np.full(20, -2.0), area)])
    assert _product(same) == pytest.approx(1.0, rel=1e-12)


def _tube_pore(spacing=1.0):
    vol = cylinder_volume(4.0, 20.0, spacing, 10.0, 6.0)
    species = bath_species(False)
    cons = {"filter": SimpleNamespace(z=-5.0), "gate": SimpleNamespace(z=5.0)}
    summary = SimpleNamespace(constrictions=cons, span=(-10.0, 10.0))
    return SimpleNamespace(st=None, summary=summary, species=species,
                           vols={s.name: vol for s in species}, elec=vol,
                           ryr=False, spacing=spacing)


def test_uncharged_tube_reads_the_diffusivity_ratio():
    pore = _tube_pore()
    r = linear_ratios(pore)
    d = {s.name: s.diffusivity for s in pore.species}
    assert r.pca_pk == pytest.approx(d["Ca2+"] / d["K+"], rel=1e-9)
    assert r.pcl_pk == pytest.approx(d["Cl-"] / d["K+"], rel=1e-9)


def test_well_ceiling_is_one_over_the_share_left_outside():
    """In a tube each plane is an equipotential, so the ceiling is
    1 / (1 - the band's share of the neutral drop). The well's edge faces
    conduct as the well (Scharfetter-Gummel weights), so the band counts
    from half a voxel to a voxel beyond its last planes."""
    pore = _tube_pore()
    vol = pore.elec
    lap = geometric_conductance(vol)
    grid = np.full(vol.mask.shape, np.nan)
    grid[vol.mask] = lap.potential
    means = np.nanmean(grid.reshape(-1, grid.shape[2]), axis=0)
    total = means[-1] - means[0]

    def expected(edge):
        inside = np.interp(edge, vol.zs, means) - np.interp(-edge, vol.zs, means)
        return 1.0 / (1.0 - abs(inside / total))

    ceiling = well_ceiling(pore, -4.0, 4.0)
    assert expected(4.5) < ceiling < expected(5.0) * 1.01
    assert well_ceiling(pore, -4.0, 4.0, 20.0) == pytest.approx(ceiling, rel=0.01)
    assert band_mask(pore, -4.0, 4.0).sum() < vol.mask.sum()
    # a well leaves the anion alone
    deep = linear_ratios(pore, extra={"Ca2+": -30.0 * band_mask(pore, -4, 4)})
    assert deep.pcl_pk == pytest.approx(linear_ratios(pore).pcl_pk, rel=1e-9)
