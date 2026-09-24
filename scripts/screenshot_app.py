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

    import tempfile

    import numpy as np

    from ip3r.render.representations import Style

    tmp = tempfile.TemporaryDirectory()
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
                win.findings.tree.setCurrentItem(win.findings._items["P6.shell_trend"])
                win.findings.show_btn.click()
                view = win.scene.view
                if view.color_by.value != "ligand_shell":
                    raise RuntimeError(f"shell check painted by {view.color_by}")
                n = len({tuple(c) for c in view.atom_colors()})
                if n < 5:
                    raise RuntimeError(f"ligand shells drew {n} colours, not 4 + grey")
                app.processEvents()                     # let the legend re-lay out
                win.grab().save(str(out / "gui_shells.png"))
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
                win.findings.tree.setCurrentItem(
                    win.findings._items["P2.cyclostome_lineages"])
                win.findings.show_btn.click()           # opens the Tree tab
                if win.tabs.currentWidget() is not win.tree:
                    raise RuntimeError("P2.cyclostome_lineages did not open the Tree tab")
            elif s == 8:
                if win.tree.root is None:
                    if not win.tree.status.text().startswith("Reading"):
                        raise RuntimeError(win.tree.status.text())
                    state["step"] -= 1
                    return QTimer.singleShot(500, step)
                info = win.tree.info
                if sorted(info["clades"].values()) != [13, 19, 19] or \
                        info["cyclostome_clades"] != [4, 2]:
                    raise RuntimeError(f"tree drew {info}")
                app.processEvents()
                win.grab().save(str(out / "gui_tree.png"))
                win.tree.zoom_vertebrates()
                app.processEvents()
                win.tree.canvas.grab().save(str(out / "tree_vertebrates.png"))
                win.findings.tree.setCurrentItem(
                    win.findings._items["P3.miss_by_contiguity"])
                win.findings.show_btn.click()           # opens the Genomes tab
                if win.tabs.currentWidget() is not win.genomes:
                    raise RuntimeError("P3.miss_by_contiguity did not open the Genomes tab")
            elif s == 9:
                g = win.genomes
                if g.grid is None or g.layer.currentData() != "miss":
                    if g.status.text() != "Not loaded." and not g.status.text().startswith("Reading"):
                        raise RuntimeError(g.status.text())
                    state["step"] -= 1
                    return QTimer.singleShot(500, step)
                info = g.info
                if info["genomes"] != 309 or info["counts"].get("missed") != 182 \
                        or info["bar_row"] != 189:
                    raise RuntimeError(f"genome grid drew {info}")
                app.processEvents()
                win.grab().save(str(out / "gui_genomes.png"))
                g.layer.setCurrentIndex(g.layer.findData("recovery"))
                g.order.setCurrentText("class, then N50")
                from ip3r.analysis.genome_grid import REACHABLE
                if g.info["counts"].get(REACHABLE) != 292:
                    raise RuntimeError(f"recovery layer drew {g.info['counts']}")
                app.processEvents()
                g.canvas.grab().save(str(out / "genomes_recovery.png"))
                win.findings.tree.setCurrentItem(win.findings._items["P1.absences"])
                win.findings.show_btn.click()           # opens the Range tab
                if win.tabs.currentWidget() is not win.range:
                    raise RuntimeError("P1.absences did not open the Range tab")
            elif s == 10:
                r = win.range
                if r.data is None:
                    if r.status.text() != "Not loaded." and not r.status.text().startswith("Reading"):
                        raise RuntimeError(r.status.text())
                    state["step"] -= 1
                    return QTimer.singleShot(500, step)
                info = r.info
                held = sum(a.holds for a in r.data.absences)
                if (info["clades"], info["with_call"], info["present"], info["proteomes"],
                        held) != (70, 45, 662, 6928, 35):
                    raise RuntimeError(f"range drew {dict(info, rows=len(info['rows']))}, "
                                       f"{held} absences held")
                app.processEvents()
                win.grab().save(str(out / "gui_range.png"))
                win.tabs.setCurrentWidget(win.dynamics)
                gt = win.dynamics.gating
                gt.model.setCurrentIndex(1)             # Mak et al. 1998
                if "Mak 1998" not in gt.text.text() or "K_inh" not in gt.note.text():
                    raise RuntimeError(f"Mak gating drew: {gt.text.text()[:120]}")
                app.processEvents()
                win.grab().save(str(out / "gui_gating_mak.png"))
                pz = win.dynamics.puffs
                pz.parent().parent().setCurrentWidget(pz)
                pz.model.setCurrentIndex(pz.model.findData("park-drive"))
                if abs(pz.coupling.value() - 0.1) > 1e-9:
                    raise RuntimeError(f"park/drive coupling {pz.coupling.value()}")
                pz.duration.setValue(5.0)
                pz.run()
            elif s == 11:
                pz = win.dynamics.puffs
                if pz.result is None:
                    if not pz.text.text().startswith("Simulating"):
                        raise RuntimeError(pz.text.text())
                    state["step"] -= 1
                    return QTimer.singleShot(500, step)
                if pz.result["coupled"]["fano"] <= pz.result["uncoupled"]["fano"]:
                    raise RuntimeError(f"park/drive puffs drew {pz.result}")
                app.processEvents()
                win.grab().save(str(out / "gui_puffs_pd.png"))
                win.tabs.setCurrentWidget(win.channel)
                win.channel.unitary_btn.click()
            elif s == 12:
                rows = win.channel.unitary_rows
                if rows is None:
                    if win.channel.unitary_btn.isEnabled():
                        raise RuntimeError("the unitary conductance worker failed")
                    state["step"] -= 1
                    return QTimer.singleShot(1000, step)
                open_ = [u.name.upper() for u in rows if u.neutral.is_conducting]
                if open_ != ["8TKF"]:
                    raise RuntimeError(f"conducting states drawn: {open_}")
                app.processEvents()
                win.grab().save(str(out / "gui_unitary.png"))
            elif s == 13:
                from ip3r.analysis.checks import all_checks, run_check
                from ip3r.parameters import PARAMETERS
                from ip3r.ui.params_dialog import ParametersDialog
                if win.params_strip.isVisible():
                    raise RuntimeError("banner shown at the documented defaults")
                d = win.params_dialog = ParametersDialog(win)
                d.show()
                p = PARAMETERS.get("gating.d1")
                if d.edit("gating.d1", "not a number") or PARAMETERS.modified:
                    raise RuntimeError("the editor accepted a non-number")
                d.edit("gating.d1", f"{p.default * 2:g}")
                app.processEvents()
                if not win.params_strip.isVisible() or "gating.d1" not in \
                        win.params_banner.text.text():
                    raise RuntimeError("no banner after an edit")
                check = next(c for c in all_checks() if c.id == "S0.c4_symmetry")
                if run_check(check).outcome.status != "not_run":
                    raise RuntimeError("a check ran against a modified registry")
                d.filter_edit.setText("gating")
                d.resize(1180, 520)
                app.processEvents()
                d.grab().save(str(out / "gui_parameters.png"))
                win.grab().save(str(out / "gui_params_banner.png"))
                d.reset_all()
                app.processEvents()
                if PARAMETERS.modified or win.params_strip.isVisible():
                    raise RuntimeError("Reset all left the registry or banner modified")
                d.close()
            elif s == 14:                            # a view worth restoring
                from ip3r.parameters import PARAMETERS
                sp = win.structure_panel
                sp.style.setCurrentIndex(sp.style.findData(Style.BACKBONE))
                sp.site_boxes["gate_lining"].setChecked(True)
                list(sp.chain_boxes.values())[-1].setChecked(False)
                win.transition.slider.setValue(5)
                cam = win.viewport.scene.camera
                cam.orbit(0.3, 0.1)
                cam.zoom(0.8)
                win.tabs.setCurrentWidget(win.channel)
                PARAMETERS.set_value("display.displacement_max", 20.0)
                state["session"] = win.sessions.save_to(Path(tmp.name) / "s.json")
                PARAMETERS.reset()
                win.structure_panel.select("6DQN")      # somewhere else entirely
            elif s == 15:
                if win.scene.structure is None or win.scene.structure.name != "6DQN":
                    state["step"] -= 1
                    return QTimer.singleShot(500, step)
                from ip3r.io.session import load_session
                if not win.sessions.apply(load_session(Path(tmp.name) / "s.json"),
                                          parameters="apply"):
                    raise RuntimeError("the session was not started")
            elif s == 16:
                ss = win.sessions
                if ss.pending is not None or ss._frame is not None:
                    if win.transition.status.text().startswith("not built"):
                        raise RuntimeError(win.transition.status.text())
                    state["step"] -= 1
                    return QTimer.singleShot(1000, step)
                skip = {"saved_at", "software_version", "format_version", "notes"}
                want = {k: v for k, v in state["session"].as_dict().items() if k not in skip}
                got = {k: v for k, v in ss.capture().as_dict().items() if k not in skip}
                def same(a, b):                     # camera floats to 1e-9 Å
                    if isinstance(a, list) and a and isinstance(a[0], float):
                        return np.allclose(a, b, rtol=0, atol=1e-9)
                    return a == b
                bad = [k for k in want if not same(want[k], got[k])]
                if win.structure_panel.current_id() != want["structure"]:
                    bad.append("deposition list selection")
                if bad:
                    raise RuntimeError("session not restored: " + "; ".join(
                        f"{k} {want[k]!r} → {got[k]!r}" if k in want else k
                        for k in bad))
                app.processEvents()
                win.grab().save(str(out / "gui_session.png"))
                from ip3r.parameters import PARAMETERS
                PARAMETERS.reset()
            else:
                print("screenshots written to", out)
                return app.quit()
        except Exception as exc:                 # report and stop
            fail(f"step {s}: {type(exc).__name__}: {exc}")
            return app.quit()
        QTimer.singleShot(1500, step)

    QTimer.singleShot(800, step)
    QTimer.singleShot(420_000, lambda: (fail("timed out"), app.quit()))
    app.exec()
    return 1 if state["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
