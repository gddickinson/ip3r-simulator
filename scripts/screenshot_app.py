#!/usr/bin/env python
"""Scripted GUI smoke test: launch, load a structure, drive panels, screenshot.

Mechanical refactors of Qt code break the GUI silently (the PIEZO1 project
learned this twice), so this runs the real application with a real OpenGL
context and exits non-zero if any step fails. It also produces the README
screenshots.

    python scripts/screenshot_app.py [--structure 6DQN] [--out docs/img]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", default="6DQN")
    ap.add_argument("--out", default="docs/img")
    ap.add_argument("--checks", action="store_true", help="also run the checks")
    args = ap.parse_args()

    from PyQt6.QtCore import QTimer
    from PyQt6.QtWidgets import QApplication

    from ip3r.config import SETTINGS
    from ip3r.ui.gl_widget import configure_surface_format
    from ip3r.ui.main_window import MainWindow
    from ip3r.ui.theme import apply_dark_theme

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    configure_surface_format(SETTINGS.render)
    app = QApplication(sys.argv[:1])
    apply_dark_theme(app)
    win = MainWindow()
    win.resize(1600, 980)
    win.show()
    state = {"errors": [], "step": 0}

    def fail(msg):
        state["errors"].append(msg)
        print("FAIL:", msg)

    def step():
        s = state["step"]
        state["step"] += 1
        try:
            if s == 0:
                win.structure_panel.select(args.structure)
            elif s == 1:
                if win.scene.structure is None:
                    state["step"] -= 1           # still loading
                    return QTimer.singleShot(500, step)
                win.structure_panel.site_boxes["ip3_contact"].setChecked(True)
                win.grab().save(str(out / "gui_element.png"))
                win.viewport.grabFramebuffer().save(str(out / "viewport_element.png"))
            elif s == 2:
                sp = win.structure_panel
                sp.color.setCurrentIndex(sp.color.findText("Conservation (JSD)"))
                win.channel.show_pore.setChecked(True)
                win.grab().save(str(out / "gui_conservation.png"))
            elif s == 3:
                win.compute_modes()
            elif s == 4:
                if win.modes.modes is None:
                    state["step"] -= 1
                    return QTimer.singleShot(500, step)
                win.modes.table.selectRow(win.modes.modes.first("A") or 0)
            elif s == 5:
                win.grab().save(str(out / "gui_modes.png"))
                win.scene.stop_animation()
                if args.checks:
                    win.findings.run_all()
            elif s == 6:
                if args.checks and not win.findings.results:
                    state["step"] -= 1
                    return QTimer.singleShot(1000, step)
                win.findings.tree.setCurrentItem(win.findings._items["S0.pore_profile"])
                win.grab().save(str(out / "gui_findings.png"))
                win.findings.tree.setCurrentItem(win.findings._items["P6.module_contrast"])
                if not win.findings.show_btn.isEnabled():
                    raise RuntimeError("P6.module_contrast cannot be shown on the structure")
                win.findings.show_btn.click()
                view = win.scene.view
                cols = {tuple(c) for c in view.highlight_rgb[view.highlight]}
                if len(cols) < 2:
                    raise RuntimeError(f"module highlight drew {len(cols)} colour(s), not 2")
                win.grab().save(str(out / "gui_modules.png"))
                tabs = win.transition.parentWidget().parentWidget()
                tabs.setCurrentWidget(win.transition)
                win.transition.preset.click()           # loads 8TKG, builds -> 8TKF
            elif s == 7:
                if win.transition.result is None:
                    if win.transition.status.text().startswith("not built"):
                        raise RuntimeError(win.transition.status.text())
                    state["step"] -= 1
                    return QTimer.singleShot(1000, step)
                win.transition.paint.setChecked(True)
                win.transition.slider.setValue(win.transition.slider.maximum())
                drawn = win.scene.view.structure.xyz
                if abs(drawn - win.scene.structure.xyz).max() < 1.0:
                    raise RuntimeError("the morph end frame did not move the model")
                win.grab().save(str(out / "gui_transition.png"))
            else:
                print("screenshots written to", out)
                return app.quit()
        except Exception as exc:                 # report and stop
            fail(f"step {s}: {type(exc).__name__}: {exc}")
            return app.quit()
        QTimer.singleShot(1500, step)

    QTimer.singleShot(800, step)
    QTimer.singleShot(240_000, lambda: (fail("timed out"), app.quit()))
    app.exec()
    return 1 if state["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
