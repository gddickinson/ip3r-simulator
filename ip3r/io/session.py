"""Saving and restoring what was being looked at.

A session records the *view*: which deposit, how it was drawn, which sites
were marked, where the camera was, which tab was open, and which transition
was built and at what frame. Reopening it re-derives everything from the
same inputs. It never stores coordinates or results. A file carrying its own
copy of the numbers would drift silently out of step with the code that made
them (the rule is the PIEZO1 simulator's, where this was first written).

Panels beyond the viewport keep their *controls* here too (``dynamics``,
``modes``, ``variants``): the Dynamics settings, which normal mode was
animating and at what amplitude, and the variants view. Controls, not
results: a restored Dynamics tab shows the settings a run used, and nothing
is re-run until asked. A mode animation is re-computed (the ANM on the
deposit) and restarted, like a transition.

One input *is* stored: the parameter overrides in force when it was saved.
Every measurement the view shows depends on them, so a session reopened
under a different set would show different numbers under the same name.
They are recorded here and compared on restore (:func:`parameter_differences`);
whether to apply them is the user's decision, never silent.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime, timezone
from pathlib import Path

from .. import __version__

__all__ = ["Session", "SESSION_FORMAT", "save_session", "load_session",
           "parameter_differences"]

#: Bumped whenever the stored fields change incompatibly.
SESSION_FORMAT = 1

#: The transition fields a session carries (the spec, never the path).
TRANSITION_KEYS = ("end", "fit", "method", "frame", "paint")

#: Panel views whose controls a session carries: flat ``name → scalar``
#: dicts written and read by the panels (``view_state`` / ``restore``).
#: Added without a format bump: an older file lacks them and opens with
#: the panels as they are; an older build drops them as unknown keys.
PANEL_VIEWS = ("dynamics", "modes", "variants")


@dataclass
class Session:
    """The state needed to put the application back where it was."""

    structure: str = ""
    #: Atom count of the deposit when saved: a re-deposited file with a
    #: different model is reported, since the camera addresses coordinates.
    n_atoms: int = 0
    style: str = "cartoon"
    color_by: str = "element"
    layer: str = ""
    show_ligands: bool = True
    #: Chains drawn; empty = all of them.
    visible_chains: list[str] = field(default_factory=list)
    sites: list[str] = field(default_factory=list)
    show_pore: bool = False
    #: The Completeness choice (``structure.graft.FILL_MODES`` key): which
    #: AlphaFold fill is drawn. The choice only; the fill is rebuilt on restore.
    completeness: str = "none"
    tab: str = ""

    #: Camera: unit quaternion (w, x, y, z), pivot, distance, pan (Å).
    camera_rotation: list[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 0.0])
    camera_pivot: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    camera_distance: float = 300.0
    camera_pan: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    #: Front clip depth from the pivot (Å) of a site view; negative = none.
    camera_slab: float = -1.0
    orthographic: bool = False

    #: ``{end, fit, method, frame, paint}`` when a transition was built.
    transition: dict = field(default_factory=dict)
    #: Dynamics tab controls (gating model, oscillation, puffs, microdomain).
    dynamics: dict = field(default_factory=dict)
    #: ``{index, amplitude}`` of the animating normal mode, if one was.
    modes: dict = field(default_factory=dict)
    #: Variants view: ``{paralog, class, layer, draw}``.
    variants: dict = field(default_factory=dict)
    #: Parameter overrides in force when saved (key → value).
    parameters: dict = field(default_factory=dict)
    notes: str = ""

    format_version: int = SESSION_FORMAT
    software_version: str = __version__
    saved_at: str = ""

    def as_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "Session":
        """Build from a parsed file, refusing what cannot be a session.

        Unknown keys are dropped, so a file from a newer minor version opens
        with what it can; a newer *format* is refused, and a known key of the
        wrong type is refused by name rather than failing later mid-restore.
        """
        if not isinstance(data, dict):
            raise ValueError("a session file holds a JSON object")
        version = data.get("format_version", 0)
        if not isinstance(version, int) or version > SESSION_FORMAT:
            raise ValueError(f"session format {version!r} is newer than this "
                             f"build understands ({SESSION_FORMAT})")
        out = cls()
        for f in fields(cls):
            if f.name in data:
                setattr(out, f.name, _checked(f.name, data[f.name],
                                              getattr(out, f.name)))
        out._validate()
        return out

    def _validate(self) -> None:
        if len(self.camera_rotation) != 4 or len(self.camera_pivot) != 3 \
                or len(self.camera_pan) != 3:
            raise ValueError("camera_rotation needs 4 numbers, pivot and pan 3")
        norm = math.sqrt(sum(v * v for v in self.camera_rotation))
        if not norm > 1e-6:
            raise ValueError("camera_rotation is not a rotation (zero quaternion)")
        self.camera_rotation = [v / norm for v in self.camera_rotation]
        if not self.camera_distance > 0:
            raise ValueError("camera_distance must be positive")
        unknown = set(self.transition) - set(TRANSITION_KEYS)
        if unknown:
            raise ValueError(f"transition has unknown keys {sorted(unknown)}")
        if self.transition and not isinstance(self.transition.get("end"), str):
            raise ValueError("transition needs an 'end' deposit id")
        if not all(isinstance(v, (int, float)) and not isinstance(v, bool)
                   for v in self.parameters.values()):
            raise ValueError("parameters must map a key to a number")
        for name in PANEL_VIEWS:
            bad = [k for k, v in getattr(self, name).items()
                   if not isinstance(v, (str, int, float, bool, type(None)))
                   or (isinstance(v, float) and not math.isfinite(v))]
            if bad:
                raise ValueError(f"{name} holds values that are not plain "
                                 f"settings: {sorted(bad)}")
        if self.modes:
            i, a = self.modes.get("index"), self.modes.get("amplitude")
            if not (isinstance(i, int) and not isinstance(i, bool) and i >= 0
                    and isinstance(a, (int, float)) and not isinstance(a, bool)
                    and a > 0):
                raise ValueError("modes needs an 'index' ≥ 0 and a positive "
                                 "'amplitude'")

    def describe(self) -> str:
        bits = [self.structure or "no structure", f"{self.style}/{self.color_by}"]
        if self.sites:
            bits.append("sites: " + ", ".join(self.sites))
        if self.transition:
            t = self.transition
            bits.append(f"transition → {t['end']} frame {t.get('frame', 0)}")
        if self.modes:
            bits.append(f"mode #{self.modes['index'] + 1} animating")
        bits.append(f"{len(self.parameters)} parameter(s) modified"
                    if self.parameters else "default parameters")
        return " · ".join(bits)


def _checked(name: str, value, default):
    """``value`` if it has the type of the field's default, else a refusal."""
    if isinstance(default, bool):
        ok = isinstance(value, bool)
    elif isinstance(default, int):
        ok = isinstance(value, int) and not isinstance(value, bool)
    elif isinstance(default, float):
        ok = isinstance(value, (int, float)) and not isinstance(value, bool) \
            and math.isfinite(value)
        value = float(value) if ok else value
    elif isinstance(default, list) and not default:        # names (chains, sites)
        ok = isinstance(value, list) and all(isinstance(v, str) for v in value)
    elif isinstance(default, list):                        # camera vectors
        ok = isinstance(value, list) and all(
            isinstance(v, (int, float)) and not isinstance(v, bool)
            and math.isfinite(v) for v in value)
        value = [float(v) for v in value] if ok else value
    else:
        ok = isinstance(value, type(default))
    if not ok:
        raise ValueError(f"session field {name!r} has the wrong type: {value!r}")
    return value


def save_session(session: Session, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    session.saved_at = datetime.now(timezone.utc).isoformat(timespec="seconds")
    session.software_version = __version__
    session.format_version = SESSION_FORMAT
    path.write_text(json.dumps(session.as_dict(), indent=1) + "\n")
    return path


def load_session(path: str | Path) -> Session:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"no session at {path}")
    return Session.from_dict(json.loads(path.read_text()))


def parameter_differences(saved: dict, current: dict) -> dict[str, tuple]:
    """Keys whose effective value differs between two override sets.

    Returns ``key → (saved, current)``, ``None`` standing for "at its
    default". Two sets that differ only in how a default is spelt (absent vs
    listed at the default) cannot occur: the registry drops an override equal
    to the default.
    """
    keys = set(saved) | set(current)
    return {k: (saved.get(k), current.get(k)) for k in sorted(keys)
            if saved.get(k) != current.get(k)}
