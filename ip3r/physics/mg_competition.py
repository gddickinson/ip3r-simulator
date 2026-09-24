"""Mg2+ and Na+ competing with Ca2+ at the RyR1 activation site.

Meissner, Rios, Tripathy & Pasek 1997 (JBC 272:1628) measured the Mg2+
constant of the activation site in the [3H]ryanodine assay Murayama's bell
uses, and fitted it as a competitive inhibitor (their Eq. 4):

    K_eff^na = Ka^na (1 + ([I]/Ki)^ni)

so half-activation moves by the na-th root of that factor. The same table
gives the monovalent cations, which compete at the same site (their Scheme
2). With both present, and each excluding the other as competitors of one
site do, the factor becomes ``1 + N + M`` with ``N = (Na/K_Na)^n_Na`` and
``M = (Mg/K_Mg)^n_Mg``. The shift Mg2+ adds in a medium that already holds
Na+ is therefore

    shift = ((1 + N + M) / (1 + N)) ^ (1/na)

Murayama's bell was measured in 0.17 M NaCl, where N is about 28: Na+
already holds most of the site, and 1 mM Mg2+ moves half-activation 2.3x,
not the 11x it does in Meissner's choline medium.

The scheme here (:class:`ryr_gating.SternParams`) shifts half-activation by
``1 + Mg/K_Mg,A`` instead. :func:`equivalent_k_mg_a` is the K_Mg,A that
gives Meissner's shift at a chosen Mg2+ (default the fibre's
``ryr.mg_free``). The two forms differ in shape, so the equivalent moves
with the Mg2+ it is matched at; :func:`equivalent_k_mg_a` reports that
directly.

The competition form is checked against the paper's own Table II
(``tests/test_mg_competition.py``): each monovalent row of Table IV predicts
the half-activation measured in 0.25 M of that salt to within 1.5x.
"""

from __future__ import annotations

from ..parameters import PARAMETERS as _P

__all__ = ["hill_ka", "ka_shift", "equivalent_k_mg_a"]


def hill_ka(ka: float, n_ca: float, conc: float, ki: float, ni: float) -> float:
    """Meissner's Eq. 4 for one competitor: the effective half-activation."""
    return ka * (1.0 + (conc / ki) ** ni) ** (1.0 / n_ca)


def _term(conc: float, ki: float, ni: float) -> float:
    return (conc / ki) ** ni if conc > 0 else 0.0


def ka_shift(mg: float, na: float | None = None) -> float:
    """Factor by which ``mg`` µM moves half-activation in ``na`` µM Na+.

    ``na`` defaults to Murayama's assay (``ryr.murayama_sodium``); pass 0
    for Meissner's own choline medium.
    """
    na = _P.value("ryr.murayama_sodium") if na is None else float(na)
    n = _term(na, _P.value("ryr.meissner_ki_na"), _P.value("ryr.meissner_ni_na"))
    m = _term(mg, _P.value("ryr.meissner_ki_mg"), _P.value("ryr.meissner_ni_mg"))
    return ((1.0 + n + m) / (1.0 + n)) ** (1.0 / _P.value("ryr.meissner_na_ca"))


def equivalent_k_mg_a(mg: float | None = None, na: float | None = None) -> float:
    """K_Mg,A for which ``1 + mg/K`` equals Meissner's shift at ``mg``."""
    mg = _P.value("ryr.mg_free") if mg is None else float(mg)
    if mg <= 0:
        raise ValueError("the equivalent is matched at a positive Mg2+")
    return mg / (ka_shift(mg, na) - 1.0)
