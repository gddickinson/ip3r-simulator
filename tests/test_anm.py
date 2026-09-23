"""Elastic-network modes on an exact C4 tetramer."""

import numpy as np

from ip3r.physics.anm import ANM, build_hessian
from conftest import c4_tetramer


def _modes(n=18):
    blocks = c4_tetramer(n_per=40, seed=3)
    anm = ANM(np.vstack(blocks), cutoff=12.0, gamma=1.0, d0=7.5)
    return anm, anm.label_symmetry(anm.calc_modes(n))


def test_hessian_has_six_zero_modes_per_component():
    """Six rigid-body modes per connected piece — the count calc_modes drops."""
    coords = np.vstack(c4_tetramer(n_per=30, seed=2))
    anm = ANM(coords, cutoff=12.0, gamma=1.0, d0=7.5)
    h = build_hessian(coords, 12.0, 1.0, 7.5).toarray()
    w = np.linalg.eigvalsh(h)
    assert np.sum(np.abs(w) < 1e-8) == 6 * anm.n_components()
    assert np.allclose(h, h.T)
    far = np.vstack([coords, coords + 500.0])          # two separate copies
    w2 = np.linalg.eigvalsh(build_hessian(far, 12.0, 1.0, 7.5).toarray())
    assert np.sum(np.abs(w2) < 1e-8) == 6 * ANM(far, cutoff=12.0).n_components()


def test_every_mode_gets_a_c4_irrep():
    _, ms = _modes()
    assert set(ms.symmetry) <= {"A", "B", "E"}
    assert np.all(ms.eigenvalues > 0)


def test_e_modes_come_in_degenerate_pairs():
    _, ms = _modes(24)
    e = np.flatnonzero(ms.symmetry == "E")
    ev = ms.eigenvalues
    # every E eigenvalue has a partner within 1e-6 relative
    for i in e:
        partners = [j for j in e if j != i and abs(ev[j] - ev[i]) < 1e-6 * ev[i]]
        assert partners, f"E mode {i} has no degenerate partner"


def test_label_can_fail_on_an_asymmetric_network():
    """Scramble one subunit: modes must stop being clean irreps."""
    blocks = c4_tetramer(n_per=40, seed=3)
    rng = np.random.default_rng(9)
    blocks[1] = blocks[1] + rng.normal(size=blocks[1].shape) * 2.0
    anm = ANM(np.vstack(blocks), cutoff=12.0, gamma=1.0, d0=7.5)
    ms = anm.label_symmetry(anm.calc_modes(12), tolerance=0.05)
    assert (ms.symmetry == "mixed").any()
