"""Round 7.30: the lumen box's luminal CaCl₂.

The box's default is Vais's own experiment, and every other choice is one
of Round 7.29's sweep points with the same baths the sweep used, so a site
drawn at 100 mM reverses where the sweep's point did. A choice that has no
Ca²⁺ experiment to replace (Cl⁻, RyR1's fixed protocol) is refused by name.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from ip3r.parameters import PARAMETERS as _P
from ip3r.physics import mole_fraction as mf
from ip3r.physics.lumen_reversal import ReversalLumen, _with_calcium
from ip3r.physics.reversal3d import experiments

from conftest import needs_structure


def _vais_ca():
    return next(e for e in experiments(False) if e.name == "Ca2+")


def test_vais_experiment_is_the_sweeps_at_10_mm():
    vais = _vais_ca()
    e = mf.experiment(_P.value("selectivity.cacl2_lumen"))
    assert (e.name, e.lumen, e.cytosol, e.impermeant) == \
        (vais.name, vais.lumen, vais.cytosol, vais.impermeant)


def test_calcium_replaced_only_in_ip3rs_ca_experiment():
    ip3r = SimpleNamespace(ryr=False)
    e = _with_calcium(ip3r, _vais_ca(), 0.1)
    assert e.lumen["Ca2+"] == 0.1 and e.cytosol == _vais_ca().cytosol
    cl = next(x for x in experiments(False) if x.name == "Cl-")
    with pytest.raises(ValueError, match="Ca2\\+ experiment only"):
        _with_calcium(ip3r, cl, 0.1)
    with pytest.raises(ValueError, match="Ca2\\+ experiment only"):
        _with_calcium(SimpleNamespace(ryr=True), experiments(True)[0], 0.1)
    for bad in (0.0, -1e-3, float("nan")):
        with pytest.raises(ValueError, match="positive"):
            _with_calcium(ip3r, _vais_ca(), bad)


def test_a_reading_knows_its_luminal_calcium():
    def rev(baths):
        return ReversalLumen(None, "Ca2+ site", "Ca2+", 0.0, {}, {}, {}, baths)
    assert rev({"K+": (0.14, 0.14), "Ca2+": (0.003, 0.0)}).ca == 0.003
    assert rev({"K+": (0.14, 0.14)}).ca == 0.0


def test_the_boxs_levels_are_the_sweeps_with_vais_marked():
    from ip3r.ui.lumen_view import _same_ca, _vais
    c = mf.concentrations()
    assert sum(_same_ca(x, _vais()) for x in c) == 1
    assert _same_ca(0.01, 0.01 * (1 + 1e-12)) and not _same_ca(0.01, 0.0101)
    assert not _same_ca(True, 1.0) and not _same_ca("0.01", 0.01)


# ------------------------------------------------------------- 8TKF, pinned
@needs_structure("8TKF")
def test_8tkf_site_at_100_mm_is_the_sweeps_point():
    """The site drawn at 100 mM reverses at the sweep's +34.38 mV
    (data/molefrac/mf_8tkf_site.out) and holds more Ca²⁺ than at 10 mM
    (2.01, test_wall_candidates). Six minutes."""
    from ip3r.io import loader
    from ip3r.physics.lumen_reversal import reversal_lumen
    from ip3r.structure.channel import measure_channel
    st = loader.load("8TKF")
    r = reversal_lumen(st, "Ca2+ site", "Ca2+", measure_channel(st), ca=0.1)
    assert r.converged and r.ca == 0.1
    assert r.v * 1e3 == pytest.approx(34.38, abs=0.05)
    assert r.held > 2.01 + 0.5
    assert np.nanmax(r.occupancy) <= 1.0
    assert "luminal CaCl2 100 mM" in r.summary()
