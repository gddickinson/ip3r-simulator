"""The GUI smoke test's steps in named groups, so a change can be checked
by the groups it touches (``screenshot_app.py --steps lumen,extras``).

Each group begins by starting the work its first step checks (the
``enter_*`` functions), so it runs the same alone or in the full sequence.
A group that needs another's state names it in ``requires``, and it is
added. Every run that leaves out ``core`` first loads the default deposit
(core's first two steps: the load and the fit checks).

A step is ``(label, fn)`` with ``fn(win, app, out, ctx) -> bool``: True
means "a worker is still running, call again".
"""

from __future__ import annotations

from dataclasses import dataclass

import screenshot_core as core
import screenshot_extras as se
import screenshot_ip3r as si
import screenshot_sparks as sk

__all__ = ["Group", "GROUPS", "select", "schedule", "describe"]


@dataclass(frozen=True)
class Group:
    name: str
    about: str
    steps: tuple
    requires: tuple = ()


def _core(*indices):
    return tuple((f"core {s}", lambda w, a, o, c, s=s: core.core_step(s, w, a, o, c))
                 for s in indices)


def _enter(fn):
    return ((fn.__name__, lambda w, a, o, c: bool(fn(w, a, o))),)


def _plain(fns):
    return tuple((f.__name__.lstrip("_"), lambda w, a, o, c, f=f: bool(f(w, a, o)))
                 for f in fns)


def _sparks():
    return tuple((f"spark {k}", lambda w, a, o, c, k=k: sk.spark_step(k, w, a, o))
                 for k in range(sk.SPARK_STEPS))


GROUPS = (
    Group("core", "load, fit, conservation, modes, findings (verdicts only with "
          "--checks), modules and shells shown", _core(*range(7))),
    Group("transition", "the 8TKG -> 8TKF morph, its gate and A-subspace headline",
          _enter(core.enter_transition) + _core(7)),
    Group("publication", "Tree, Genomes and Range tabs opened from their checks",
          _enter(core.enter_publication) + _core(8, 9, 10)),
    Group("dynamics", "Mak's bell and park/drive puffs",
          _enter(core.enter_dynamics) + _core(11)),
    Group("unitary", "the ITPR3 state panel's conductances",
          _enter(core.enter_unitary) + _core(12)),
    Group("params", "the parameter editor and the modified banner", _core(13)),
    Group("session", "a session saved and restored, variants, AlphaFold fills, "
          "9YKK refused", _core(14, 15, 16, 17, 18, 19), requires=("transition",)),
    Group("ryr", "RyR1 9HEO: channel, morph preset, mutants, gating",
          _enter(core.enter_ryr) + _core(20, 21)),
    Group("sparks", "RyR1 sparks: mean field, cleft, fitted, Mg2+",
          _enter(core.enter_sparks) + _sparks()),
    Group("models", "the gating comparison and the microdomain puffs",
          _plain((si._gating, si._domain))),
    Group("views", "tree beside --bnni, lesion layer, S23 genomes, VUS bands",
          _plain((si._tree_pair_start, si._tree_pair, si._lesion_start, si._lesion,
                  si._range_genomes, si._vus_bands)), requires=("publication",)),
    Group("ratfill", "7LHF filled through the alignment",
          _plain((si._rat_fill_start, si._rat_fill))),
    Group("lumen", "8TKF's lumen: neutral, paired PB, dielectric, + image, K+ energy "
          "(the slow group: minutes)",
          _plain((si._lumen_start, si._lumen, si._lumen_charged, si._lumen_dielectric,
                  si._lumen_image))),
    Group("extras", "HUD, selection, sequence, context menu, distance, layout, "
          "Analyses menu, full screen", _plain(se._STEPS)),
)

_BY_NAME = {g.name: g for g in GROUPS}
# Every step of the three step modules belongs to exactly one group.
assert sum(len(g.steps) for g in GROUPS if g.name in (
    "models", "views", "ratfill", "lumen")) == si.IP3R_STEPS
assert len(_BY_NAME["extras"].steps) == se.EXTRAS_STEPS


def select(spec: str | None) -> list[Group]:
    """The groups ``spec`` names (comma-separated; None or "all" = every
    group), with what they require, in the full run's order."""
    if not spec or spec == "all":
        return list(GROUPS)
    names = {n.strip() for n in spec.split(",") if n.strip()}
    unknown = names - set(_BY_NAME)
    if unknown:
        raise SystemExit(f"unknown step group(s) {sorted(unknown)}; "
                         f"known: {', '.join(_BY_NAME)}")
    todo = list(names)
    while todo:
        for r in _BY_NAME[todo.pop()].requires:
            if r not in names:
                names.add(r)
                todo.append(r)
    return [g for g in GROUPS if g.name in names]


def schedule(groups: list[Group]) -> list[tuple[str, str, object]]:
    """``(group, label, fn)`` in order; the default deposit is loaded first
    when ``core`` is not among the groups."""
    out = []
    if not any(g.name == "core" for g in groups):
        out += [("load", label, fn) for label, fn in _core(0, 1)]
    for g in groups:
        out += [(g.name, label, fn) for label, fn in g.steps]
    return out


def describe() -> str:
    return "\n".join(f"  {g.name:12s} {g.about}"
                     + (f" (adds {', '.join(g.requires)})" if g.requires else "")
                     for g in GROUPS)
