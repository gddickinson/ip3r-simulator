"""The De Young-Keizer / Li-Rinzel gating curves have the shapes the
literature measured, and the model constants are the registered ones."""

import numpy as np

from ip3r.physics.gating import (GatingParams, bell_peak, h_inf,
                                 hill_fit_left_flank, open_probability)


def test_bell_shape():
    c = np.logspace(-3, 3, 500)
    po = open_probability(c, 1.0)
    i = int(np.argmax(po))
    assert 0 < i < len(c) - 1                       # interior maximum
    assert po[0] < 1e-3 * po[i] and po[-1] < 1e-3 * po[i]


def test_ip3_raises_and_shifts_the_peak():
    peaks = [bell_peak(p) for p in (0.1, 1.0, 10.0)]
    assert peaks[0][1] < peaks[1][1] < peaks[2][1]
    assert peaks[0][0] < peaks[1][0] < peaks[2][0]


def test_ip3_tunes_inhibition_more_than_activation():
    """Right flank moves more than the left as IP3 rises 100-fold.

    Measured: 2.4x against 1.8x for 0.1 -> 10 uM. The first version of this
    test demanded 3x the left shift, encoding the experimental statement that
    IP3 tunes inhibition *only* (Mak et al. 1998) — which the DYK model does
    not reproduce, and the prose was corrected to say so.
    """
    left = [hill_fit_left_flank(p) for p in (0.1, 10.0)]
    right = []
    for p in (0.1, 10.0):
        c = np.logspace(-1, 3, 4000)
        po = open_probability(c, p)
        cp, pk = bell_peak(p)
        after = c > cp
        right.append(float(np.interp(-pk / 2, -po[after], c[after])))
    assert right[1] / right[0] > 1.2 * (left[1] / left[0]) > 1.2


def test_h_inf_falls_with_calcium():
    assert h_inf(10.0, 0.5) < h_inf(0.1, 0.5)


def test_constants_are_registered():
    g = GatingParams()
    assert (g.d1, g.d2, g.d3, g.d5, g.a2) == (0.13, 1.049, 0.9434, 0.08234, 0.2)
