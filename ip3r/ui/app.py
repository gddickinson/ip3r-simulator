"""Application launcher: argument parsing and the Qt event loop."""

from __future__ import annotations

import argparse
import sys

__all__ = ["main"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m ip3r gui",
                                     description="IP3R structural simulator")
    parser.add_argument("--structure", metavar="PDB", default=None,
                        help="structure to load at startup (default 6DQN)")
    parser.add_argument("--geometry", metavar="WxH")
    parser.add_argument("--no-load", action="store_true",
                        help="start without loading a structure")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    from PyQt6.QtCore import QTimer
    from PyQt6.QtWidgets import QApplication

    from ..config import SETTINGS
    from .gl_widget import configure_surface_format
    from .main_window import MainWindow
    from .theme import apply_dark_theme

    configure_surface_format(SETTINGS.render)
    app = QApplication(sys.argv[:1])
    apply_dark_theme(app)
    win = MainWindow()
    if args.geometry:
        w, h = (int(v) for v in args.geometry.lower().split("x"))
        win.resize(w, h)
    win.show()
    if not args.no_load:
        pdb = (args.structure or SETTINGS.default_structure).upper()
        QTimer.singleShot(300, lambda: win.structure_panel.select(pdb))
    return app.exec()
