"""Round 7.19: the charge–space fluid on the voxelised lumen, calibrated
before it is read.

The wall's groups carry the fixed map's charge; with its terms off the
fluid is Donnan (local) and Poisson–Boltzmann (pb) exactly; in a long
charged tube the pb closure reaches the local one (the Donnan limit), and
beyond the charge it carries a field the local closure cannot. The ruler
(P = D x the weighted Laplace conductance) is the 1-D series by hand on a
tube.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import csc3d
from ip3r.physics.charge3d import (donnan_field, local_density,
                                   poisson_boltzmann, slice_density)
from ip3r.physics.csc_readings import _mixed
from ip3r.physics.ohmic3d import bernoulli_weight, cylinder_volume
from ip3r.physics.pore_charge import ChargedGroup
from ip3r.physics.selectivity3d import Pore, _read

from test_charged3d import _tube

SP = _mixed(0.25, 0.01)


def _acid(z, res, x=0.0):
    return ChargedGroup(res, "ASP", "A", z, x, 0.0, -1.0)


def _band(vol, lo=-20.0, hi=20.0, n=12):
    groups = [_acid(z, i) for i, z in enumerate(np.linspace(lo, hi, n))]
    return groups, slice_density(vol, groups)


def test_wall_fluid_carries_the_fixed_charge():
    vol = cylinder_volume(4.0, 30.0, 0.5, 10.0, 6.0)
    groups = [_acid(-5.0, 1), ChargedGroup(2, "LYS", "A", 4.0, 0.0, 0.0, 1.0),
              _acid(9.0, 3)]
    x = slice_density(vol, groups)
    w = csc3d.wall_fluid(vol, groups, placement="slice")
    assert np.allclose(w.charge(), x[vol.mask], rtol=1e-12, atol=1e-9)
    assert w.names == ("O", "N")
    pos = np.array([[4.5, 0.0, -5.0], [0.0, 4.5, 4.0], [-4.5, 0, 9.0]])
    xl, _ = local_density(vol, pos, np.array([-1.0, 1.0, -1.0]))
    wl = csc3d.wall_fluid(vol, groups, pos, placement="local")
    assert np.allclose(wl.charge(), xl[vol.mask], rtol=1e-12, atol=1e-9)
    # two oxygens per acid, each carrying half its charge
    o = wl.density[0] > 0
    assert np.allclose(wl.valence[0][o], -0.5)
    with pytest.raises(ValueError):
        csc3d.wall_fluid(vol, groups, placement="pb")


def test_terms_off_is_donnan_and_pb_exactly():
    vol = _tube(n=12, length=80.0, h=0.5)
    groups, x = _band(vol, -10, 10, 6)
    w = csc3d.wall_fluid(vol, groups, placement="slice")
    off = csc3d.local_csc(vol, x, SP, w, hs=False, use_msa=False)
    d = donnan_field(vol, x, SP)
    assert np.abs(off.potential - d).max() < 1e-10
    assert np.abs(off.energies["Ca2+"] - 2 * d).max() < 1e-10
    u, _, _ = poisson_boltzmann(vol, x, SP, permittivity=40.0)
    pb = csc3d.pb_csc(vol, x, SP, w, permittivity=40.0, hs=False, use_msa=False)
    assert pb.converged
    assert np.abs(pb.potential - u).max() < 1e-9


def test_pb_reaches_the_local_fluid_in_a_long_charged_tube():
    vol = _tube(n=12, length=120.0, h=0.5)
    groups, x = _band(vol)
    w = csc3d.wall_fluid(vol, groups, placement="slice")
    loc = csc3d.local_csc(vol, x, SP, w)
    pb = csc3d.pb_csc(vol, x, SP, w, permittivity=40.0, reference=loc)
    assert loc.converged and pb.converged
    mid = int(np.argmin(np.abs(vol.zs)))
    for name in ("K+", "Ca2+", "Cl-"):
        assert pb.energies[name][5, 5, mid] == pytest.approx(
            loc.energies[name][5, 5, mid], rel=1e-5)
    # uniform in-plane: the tube is 1-D
    assert np.ptp(loc.energies["Ca2+"][:, :, mid]) < 1e-9


def test_the_fluid_binds_calcium_beyond_donnan():
    """Gillespie's point, in the voxels: at the same wall the csc fluid
    favours Ca2+ over K+ by more than the mean potential alone."""
    vol = _tube(n=12, length=120.0, h=0.5)
    groups, x = _band(vol)
    w = csc3d.wall_fluid(vol, groups, placement="slice")
    loc = csc3d.local_csc(vol, x, SP, w)
    d = donnan_field(vol, x, SP)
    mid = int(np.argmin(np.abs(vol.zs)))
    csc_adv = 2 * loc.energies["K+"][5, 5, mid] - loc.energies["Ca2+"][5, 5, mid]
    don_adv = 2 * d[5, 5, mid] - 2 * d[5, 5, mid]       # zero: z u for both
    assert csc_adv > don_adv + 1.0


def test_pb_carries_the_field_past_the_charge():
    """Uncharged tube beyond a charged band: local neutrality leaves it at
    the bath; Poisson's field reaches into it (the gate's question)."""
    vol = _tube(n=12, length=120.0, h=0.5)
    groups, x = _band(vol, -40.0, -20.0, 6)
    w = csc3d.wall_fluid(vol, groups, placement="slice")
    loc = csc3d.local_csc(vol, x, SP, w)
    pb = csc3d.pb_csc(vol, x, SP, w, permittivity=40.0, reference=loc)
    # the first plane past the band whose charge is 1e-4 of the band's
    line = np.abs(x[5, 5, :])
    k = int(np.flatnonzero((vol.zs > -20.0) & (line < 1e-4 * line.max()))[0])
    assert pb.energies["Ca2+"][5, 5, k] < -0.1
    assert abs(pb.energies["Ca2+"][5, 5, k]) > 20 * abs(loc.energies["Ca2+"][5, 5, k])


def _pore(vol):
    summary = SimpleNamespace(constrictions={"filter": SimpleNamespace(z=-5.0),
                                             "gate": SimpleNamespace(z=5.0)})
    return Pore(None, summary, SP, {s.name: vol for s in SP}, vol, True, vol.spacing)


def test_the_ruler_is_the_1d_series_on_a_tube():
    vol = _tube(n=4, length=30.0, h=0.5)
    pore = _pore(vol)
    neutral = _read(pore, "neutral", None)
    d = {s.name: s.diffusivity for s in SP}
    assert neutral.ratio == pytest.approx(d["Ca2+"] / d["K+"], rel=1e-8)
    rng = np.random.default_rng(3)
    planes = {s.name: rng.normal(0, 1.0, len(vol.zs)) for s in SP}
    energies = {n: np.broadcast_to(e, vol.mask.shape).copy() for n, e in planes.items()}
    got = _read(pore, "x", energies)

    def series(name):
        w = bernoulli_weight(planes[name][:-1], planes[name][1:])
        return d[name] / np.sum(1.0 / w)
    assert got.ratio == pytest.approx(series("Ca2+") / series("K+"), rel=1e-7)
    # the shares of a uniform tube are its length fraction
    half = _P.value("lumen.constriction_half_width")
    assert neutral.shares["gate"]["K+"] == pytest.approx(
        2 * half / (vol.zs[-1] - vol.zs[0]), rel=0.05)


# ------------------------------------------------------------- real deposits
def _deposit(pdb):
    from ip3r.io import loader
    try:
        return loader.load(pdb)
    except Exception:
        pytest.skip(f"{pdb} not fetched")


def test_real_9heo_the_gate_still_caps_p_ca():
    """RyR1: the filter binds Ca2+ (csc), the 3-D field reaches the gate
    and lifts P_Ca:P_K past the 1-D reading, and still the gate window holds
    over half of Ca2+'s resistance and the ratio stays far below Xu's 7.0.
    (1 Å grid: within 2 % of 0.5 Å.)"""
    from ip3r.physics.selectivity3d import selectivity_3d
    r = selectivity_3d(_deposit("9HEO"), readings=("neutral", "pb", "pb + csc"),
                       spacing=1.0)
    assert r.converged and r.ryr and r.unreached == []
    neutral, pb, full = (r.readings[k] for k in ("neutral", "pb", "pb + csc"))
    assert 0.5 < neutral.ratio < 0.7
    assert full.peak["Ca2+"] > 5 * full.peak["K+"]       # the filter binds Ca2+
    assert pb.peak["Ca2+"] < 1.5 * pb.peak["K+"]         # point ions do not
    assert r.one_d["1-D csc"] < full.ratio < 1.5
    assert full.shares["gate"]["Ca2+"] > 0.5
    assert full.ratio < r.measured / 5


def test_a_wall_that_fills_its_voxels_is_refused():
    """Four acids' oxygens in a 2 Å tube over a narrow Gaussian leave no
    room: refused by name, before any solve."""
    from ip3r.physics.csc_readings import moved
    vol = _tube(n=4, length=20.0, h=0.5)
    groups = [_acid(0.0, i) for i in range(4)]
    with moved("pore_charge.smoothing", 0.5):
        x = slice_density(vol, groups)
        w = csc3d.wall_fluid(vol, groups, placement="slice")
    with pytest.raises(ValueError, match="no fluid fits.*smoothing"):
        csc3d.local_csc(vol, x, SP, w)
