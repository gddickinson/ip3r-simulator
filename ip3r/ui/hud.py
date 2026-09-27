"""Heads-up display over the viewport (Round 7.18, ported from PIEZO1).

What a screenshot needs and a bare 3-D view does not give it:

* a **scale bar** in round Ångström. In perspective it is exact only in the
  plane through the pivot and says "at the pivot"; in orthographic (O) it
  holds at every depth and the caveat is dropped;
* an **orientation gnomon**: the four-fold axis with its cytosolic end
  labelled, beside the deposit's x, y, z. "Cytosol up" is the one
  orientation this receptor is read in;
* **readouts** (the selection, a measured distance, the morph frame), so
  the number and the thing it describes are in the same frame;
* an amber **provenance banner** whenever part of what is drawn is an
  AlphaFold prediction. Not switchable: it is the difference between a
  measurement and a model.

QPainter on a transparent child widget, as the labels are: mixing 2-D
drawing into the GL pass fights the impostors' depth writes.
"""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import QWidget

from ..render.camera import quat_to_matrix

__all__ = ["HudOverlay", "HudSettings", "nice_scale_length", "NICE_LENGTHS"]

#: Lengths the scale bar may take (Å): always a round number.
NICE_LENGTHS = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
_AMBER, _TEXT, _DIM = "#ffb454", "#f0f3f8", "#c8ccd4"


def nice_scale_length(target: float) -> float:
    """Largest tidy length not exceeding ``target`` (Å); the smallest if none."""
    usable = [v for v in NICE_LENGTHS if v <= target]
    return float(usable[-1] if usable else NICE_LENGTHS[0])


@dataclass
class HudSettings:
    """What the overlay shows; each element switchable (View → HUD)."""

    title: bool = True
    scale_bar: bool = True
    gnomon: bool = True
    readouts: bool = True

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, data: dict | None) -> "HudSettings":
        names = {f.name for f in fields(cls)}
        return cls(**{k: bool(v) for k, v in (data or {}).items() if k in names})


class HudOverlay(QWidget):
    """Transparent overlay: title, scale bar, gnomon, readouts, provenance."""

    MARGIN = 14

    def __init__(self, viewport) -> None:
        super().__init__(viewport)
        self.viewport = viewport
        self.settings = HudSettings()
        self.title = ""
        #: key -> text, drawn top-left under the title in insertion order.
        self.readouts: dict[str, str] = {}
        #: Set while part of the drawing is prediction; not a preference.
        self.provenance = ""
        #: The four-fold axis (unit, deposit coordinates, toward the cytosol).
        self.axis: np.ndarray | None = None
        for attr in (Qt.WidgetAttribute.WA_TransparentForMouseEvents,
                     Qt.WidgetAttribute.WA_NoSystemBackground,
                     Qt.WidgetAttribute.WA_TranslucentBackground):
            self.setAttribute(attr, True)
        self.setGeometry(viewport.rect())

    # ------------------------------------------------------------------ data

    def set_readout(self, key: str, text: str) -> None:
        """Record a readout; empty text removes it."""
        if text:
            self.readouts[key] = text
        else:
            self.readouts.pop(key, None)
        self.update()

    def set_title(self, text: str, axis=None) -> None:
        self.title = text
        self.axis = None if axis is None else np.asarray(axis, float)
        self.readouts.clear()
        self.update()

    def set_provenance(self, text: str) -> None:
        self.provenance = text
        self.update()

    # ---------------------------------------------------------------- scale

    def world_per_pixel(self) -> float:
        """Å per logical pixel in the plane through the camera's pivot."""
        scene = getattr(self.viewport, "scene", None)
        if scene is None or self.height() <= 0:
            return 0.0
        cam = scene.camera
        return float(2.0 * cam.distance * np.tan(np.radians(cam.fov) / 2.0) / self.height())

    def scale_bar(self) -> tuple[float, float] | None:
        """(length Å, length px) of the bar, or None when none fits."""
        scale = self.world_per_pixel()
        if scale <= 0:
            return None
        length = nice_scale_length(0.2 * self.width() * scale)
        pixels = length / scale
        if pixels < 24 or pixels > 0.8 * self.width():
            return None
        return length, pixels

    # ---------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:                # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        font = QFont()
        font.setPointSize(11)
        p.setFont(font)
        s = self.settings
        top = self.MARGIN
        if s.title and self.title:
            top = self._title(p)
        if s.readouts and self.readouts:
            for i, text in enumerate(self.readouts.values()):
                self._text(p, QRectF(self.MARGIN, top + 19 * i, self.width(), 19),
                           Qt.AlignmentFlag.AlignLeft, text, _DIM)
        if s.scale_bar:
            self._scale(p)
        if s.gnomon:
            self._gnomon(p)
        if self.provenance:
            font.setBold(True)
            p.setFont(font)
            self._text(p, QRectF(self.MARGIN, self.height() - 34,
                                 self.width() - 2 * self.MARGIN, 20),
                       Qt.AlignmentFlag.AlignRight, self.provenance, _AMBER)
        p.end()

    def _text(self, p, rect, align, text, colour=_TEXT) -> None:
        """Text over a dark offset copy, legible on any background."""
        flags = align | Qt.AlignmentFlag.AlignTop
        p.setPen(QPen(QColor(0, 0, 0, 190)))
        p.drawText(rect.translated(1.0, 1.0), flags, text)
        p.setPen(QPen(QColor(colour)))
        p.drawText(rect, flags, text)

    def _title(self, p) -> float:
        font = p.font()
        font.setPointSize(15)
        font.setBold(True)
        p.setFont(font)
        self._text(p, QRectF(self.MARGIN, self.MARGIN, self.width(), 26),
                   Qt.AlignmentFlag.AlignLeft, self.title)
        font.setPointSize(11)
        font.setBold(False)
        p.setFont(font)
        return self.MARGIN + 28

    def _scale(self, p) -> None:
        bar = self.scale_bar()
        if bar is None:
            return
        length, pixels = bar
        x, y = self.MARGIN, self.height() - self.MARGIN - 24
        for pen in (QPen(QColor(0, 0, 0, 190), 5.0), QPen(QColor(_TEXT), 3.0)):
            pen.setCapStyle(Qt.PenCapStyle.FlatCap)
            p.setPen(pen)
            p.drawLine(QPointF(x, y), QPointF(x + pixels, y))
            for end in (x, x + pixels):
                p.drawLine(QPointF(end, y - 5), QPointF(end, y + 5))
        self._text(p, QRectF(x, y + 5, max(pixels, 140), 18), Qt.AlignmentFlag.AlignLeft,
                   self.scale_label(length))

    def scale_label(self, length: float) -> str:
        """In perspective the bar is exact only in the pivot's plane, and says
        so; an orthographic view has one scale at every depth."""
        scene = getattr(self.viewport, "scene", None)
        ortho = scene is not None and getattr(scene.camera, "orthographic", False)
        return f"{length:.0f} Å" + ("" if ortho else " at the pivot")

    def _gnomon(self, p) -> None:
        """The deposit's x, y, z and the four-fold axis, turned with the camera."""
        scene = getattr(self.viewport, "scene", None)
        if scene is None:
            return
        rot = quat_to_matrix(scene.camera.rotation)
        size = 28.0
        o = QPointF(self.width() - self.MARGIN - size - 18,
                    self.height() - self.MARGIN - size - 40)
        arrows = [(np.eye(3)[k], c, n) for k, c, n in
                  ((0, "#f26d6d", "x"), (1, "#7ed67e", "y"), (2, "#6fb1ff", "z"))]
        if self.axis is not None:
            arrows.append((self.axis * 1.35, _AMBER, "cyt"))
        for vec, colour, name in arrows:
            local = rot @ vec
            tip = QPointF(o.x() + float(local[0]) * size, o.y() - float(local[1]) * size)
            p.setPen(QPen(QColor(colour), 2.5 if name == "cyt" else 1.8))
            p.drawLine(o, tip)
            self._text(p, QRectF(tip.x() - 12, tip.y() - 16, 30, 16),
                       Qt.AlignmentFlag.AlignLeft, name, colour)
