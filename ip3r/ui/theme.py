"""The dark application palette, matched to the viewport background."""

from __future__ import annotations

from PyQt6.QtGui import QColor, QPalette

__all__ = ["apply_dark_theme", "STATUS_COLORS"]

#: Colour per check outcome, used by the findings panel.
STATUS_COLORS = {"confirmed": "#72cc80", "discrepancy": "#f28c4d",
                 "not_run": "#8a8f99", "error": "#e05561"}


def apply_dark_theme(app) -> None:
    app.setStyle("Fusion")
    p = QPalette()
    base, text = QColor(24, 27, 34), QColor(215, 219, 227)
    p.setColor(QPalette.ColorRole.Window, QColor(30, 33, 41))
    p.setColor(QPalette.ColorRole.WindowText, text)
    p.setColor(QPalette.ColorRole.Base, base)
    p.setColor(QPalette.ColorRole.AlternateBase, QColor(34, 38, 47))
    p.setColor(QPalette.ColorRole.Text, text)
    p.setColor(QPalette.ColorRole.Button, QColor(40, 44, 54))
    p.setColor(QPalette.ColorRole.ButtonText, text)
    p.setColor(QPalette.ColorRole.Highlight, QColor(70, 120, 200))
    p.setColor(QPalette.ColorRole.HighlightedText, QColor(255, 255, 255))
    p.setColor(QPalette.ColorRole.ToolTipBase, base)
    p.setColor(QPalette.ColorRole.ToolTipText, text)
    app.setPalette(p)
