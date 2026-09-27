"""The GUI smoke test's first 22 steps (screenshot_app.py): load,
fit, modes, findings, the morph, the publication tabs, puffs, unitary
conductances, the parameter editor, a session round trip, variants,
AlphaFold fills, and RyR1 9HEO.

core_step(s, win, app, out, ctx) runs step s and returns True
while a worker is still running (the caller retries the same step). The
enter_* functions start the work a group's first step checks, so a
group can run on its own (screenshot_groups.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import screenshot_session as ss_
import screenshot_view as sv
from screenshot_ip3r import check_transition_headline
from ip3r.render.representations import Style

__all__ = ["core_step", "CORE_STEPS"] + [
    "enter_transition", "enter_publication", "enter_dynamics", "enter_unitary",
    "enter_ryr", "enter_sparks"]

CORE_STEPS = 22


def core_step(s: int, win, app, out, ctx) -> bool:
    """ctx: args (structure, checks), state (the saved session),
    tmp (a temporary directory)."""
    args, state, tmp = ctx.args, ctx.state, ctx.tmp
    if s == 0:
        win.structure_panel.select(args.structure)
    elif s == 1:
        if win.scene.structure is None:
            return True
        print(sv.check_fit(win, app), file=sys.stderr)
        sv.check_list(win)
        sv.check_subunit_refit(win, app)
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
            return True
        win.modes.table.selectRow(win.modes.modes.first("A") or 0)
    elif s == 5:
        win.grab().save(str(out / "gui_modes.png"))
        win.scene.stop_animation()
        if args.checks:
            win.findings.run_all()
    elif s == 6:
        if args.checks and not win.findings.results:
            return True
        win.findings.tree.setCurrentItem(win.findings._items["S0.pore_profile"])
        if args.checks:                          # only with verdicts drawn
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
        sv.check_site_view(win, app)
        win.grab().save(str(out / "gui_shells.png"))
    elif s == 7:
        if win.transition.result is None:
            if win.transition.status.text().startswith("not built"):
                raise RuntimeError(win.transition.status.text())
            return True
        win.transition.paint.setChecked(True)
        win.transition.slider.setValue(win.transition.slider.maximum())
        drawn = win.scene.view.structure.xyz
        if abs(drawn - win.scene.structure.xyz).max() < 1.0:
            raise RuntimeError("the morph end frame did not move the model")
        g = win.transition.result.gate
        if not (g.gate[0] < 3.0 < 5.5 < g.gate[-1]) or "gate" not in \
                win.transition.frame_label.text():
            raise RuntimeError(f"gate along the morph: {g.gate[[0, -1]]}")
        check_transition_headline(win)
        win.grab().save(str(out / "gui_transition.png"))
    elif s == 8:
        if win.tree.root is None:
            if not win.tree.status.text().startswith("Reading"):
                raise RuntimeError(win.tree.status.text())
            return True
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
            return True
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
            return True
        info = r.info
        held = sum(a.holds for a in r.data.absences)
        if (info["clades"], info["with_call"], info["present"], info["proteomes"],
                held) != (70, 45, 662, 6928, 35):
            raise RuntimeError(f"range drew {dict(info, rows=len(info['rows']))}, "
                               f"{held} absences held")
        app.processEvents()
        win.grab().save(str(out / "gui_range.png"))
    elif s == 11:
        pz = win.dynamics.puffs
        if pz.result is None:
            if not pz.text.text().startswith("Simulating"):
                raise RuntimeError(pz.text.text())
            return True
        if pz.result["coupled"]["fano"] <= pz.result["uncoupled"]["fano"]:
            raise RuntimeError(f"park/drive puffs drew {pz.result}")
        app.processEvents()
        win.grab().save(str(out / "gui_puffs_pd.png"))
    elif s == 12:
        rows = win.channel.unitary_rows
        if rows is None:
            if win.channel.unitary_btn.isEnabled():
                raise RuntimeError("the unitary conductance worker failed")
            return True
        open_ = [u.name.upper() for u in rows if u.neutral.is_conducting]
        if open_ != ["8TKF"]:
            raise RuntimeError(f"conducting states drawn: {open_}")
        u = next(u for u in rows if u.neutral.is_conducting)
        if not u.paired_charge.bridged or not (
                0 < u.paired.conductance_pS < u.charged.conductance_pS):
            raise RuntimeError(f"salt-bridged reading not drawn: {u.row()}")
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
        ss_.set_panel_view(win)
        state["session"] = win.sessions.save_to(Path(tmp.name) / "s.json")
        ss_.disturb_panel_view(win)
        PARAMETERS.reset()
        win.structure_panel.select("6DQN")      # somewhere else entirely
    elif s == 15:
        if win.scene.structure is None or win.scene.structure.name != "6DQN":
            return True
        from ip3r.io.session import load_session
        if not win.sessions.apply(load_session(Path(tmp.name) / "s.json"),
                                  parameters="apply"):
            raise RuntimeError("the session was not started")
    elif s == 16:
        ss = win.sessions
        if ss.pending is not None or ss._frame is not None \
                or ss._mode is not None:
            if win.transition.status.text().startswith("not built"):
                raise RuntimeError(win.transition.status.text())
            return True
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
    elif s == 17:                            # variants, on the restored view
        if ss_.mode_round_trip(win, state):
            return True
        ss_.check_follow(win)
        v = win.variants
        win.tabs.setCurrentWidget(v)
        if v.paralog.currentText() != "ITPR3":
            raise RuntimeError("the Variants tab did not follow the ITPR3 deposit")
        v.bucket.setCurrentText("all")
        v.layer.setCurrentIndex(v.layer.findData("family"))
        v.draw.setChecked(True)
        batch = win.scene.scene.get("variants")
        st, view = win.scene.structure, win.scene.view
        idx = win.scene._variant_atoms
        if batch is None or batch.count == 0:
            raise RuntimeError("no variant spheres drawn")
        shown = set(st.chain[idx])
        if shown != set(view.visible_chains or st.chains):
            raise RuntimeError(f"spheres on chains {sorted(shown)}, "
                               f"not the visible {sorted(view.visible_chains)}")
        if abs(view.structure.xyz[idx] - st.xyz[idx]).max() < 0.1:
            raise RuntimeError("spheres at the deposit, not the morph frame")
        strata = {row for row in range(v.table.rowCount())
                  if v.table.item(row, 6).text()}
        if not strata:
            raise RuntimeError("no VUS stratum in the table with a layer chosen")
        app.processEvents()
        win.grab().save(str(out / "gui_variants.png"))
        v.paralog.setCurrentText("ITPR1")       # wrong numbering: must refuse
        if win.scene.scene.get("variants") is not None or \
                "not in human ITPR1" not in v.status.text():
            raise RuntimeError("variants drawn on a deposit in another numbering")
        v.draw.setChecked(False)
        win.structure_panel.set_completeness("gaps")   # AlphaFold fills
    elif s == 18:
        fc, sc = win.fills, win.scene
        if fc.model is None:
            if fc.message:
                raise RuntimeError(fc.message)
            return True
        m = fc.model
        if m.prediction != "AF-Q14573-F1" or len(m.fills) != 56 or m.skipped:
            raise RuntimeError(f"8TKG fill: {m.summary()}")
        fv, seams = sc.fill.view, sc.scene.get("seams")
        visible = sc.view.visible_chains or set(sc.structure.chains)
        n_seams = sum(len(f.seam_fill) for f in m.fills if f.stretch.chain in visible)
        if fv is None or seams is None or seams.count != n_seams:
            raise RuntimeError(f"fill or seams not drawn ({n_seams} seams expected)")
        if abs(fv.structure.xyz - m.atoms.xyz).max() < 0.5:
            raise RuntimeError("the fill sits at the deposit, not the morph frame")
        if "AlphaFold fill" not in win.structure_panel.legend.text():
            raise RuntimeError("no pLDDT legend with a fill drawn")
        held = win.sessions.capture().completeness
        if held != "gaps":
            raise RuntimeError(f"session holds completeness {held!r}")
        sp = win.structure_panel
        sp.style.setCurrentIndex(sp.style.findData(Style.CARTOON))
        sp.color.setCurrentIndex(sp.color.findText("Uniform"))   # fill stands out
        win.tabs.setCurrentWidget(win.channel)
        app.processEvents()
        win.grab().save(str(out / "gui_alphafold.png"))
        win.viewport.grabFramebuffer().save(str(out / "viewport_alphafold.png"))
        win.transition.slider.setValue(0)
        if abs(fv.structure.xyz - m.atoms.xyz).max() > 0.05:
            raise RuntimeError("at frame 0 the fill is not where it was built")
        win.structure_panel.select("9YKK")      # ITPR2: no model on either route
    elif s == 19:
        fc = win.fills
        if win.scene.structure is None or win.scene.structure.name != "9YKK" \
                or (fc.model is None and "9YKK" not in fc.message):
            return True
        if not fc.message.startswith("not filled") or fc.model is not None \
                or win.scene.scene.get("fill:ribbon") is not None:
            raise RuntimeError(f"9YKK was filled: {fc.message}")
    elif s == 20:
        if win.scene.structure is None or win.scene.structure.name != "9HEO":
            return True
        ch, tp = win.channel, win.transition
        if ch.panel_paralog != "RYR1" or not ch.mutants_btn.isVisibleTo(ch):
            raise RuntimeError(f"channel panel not on RyR1: {ch.panel_paralog}")
        if "9R8O" not in tp.preset.text() or tp.end.findData("9R8O") < 0:
            raise RuntimeError(f"RyR1 morph preset missing: {tp.preset.text()}")
        app.processEvents()
        win.viewport.grabFramebuffer().save(str(out / "viewport_ryr1.png"))
        win.tabs.setCurrentWidget(ch)
        ch.mutants_btn.click()
    elif s == 21:
        rows = win.channel.mutant_rows
        if rows is None:
            if win.channel.mutants_btn.isEnabled():
                raise RuntimeError("the RyR1 mutant worker failed")
            return True
        d = next(r for r in rows if r.name == "D4899Q")
        if not (d.bridged and d.measured_ratio < 0.3 < d.paired_ratio):
            raise RuntimeError(f"RyR1 mutants drawn wrong: {d.row()}")
        app.processEvents()
        win.grab().save(str(out / "gui_ryr_mutants.png"))
        win.tabs.setCurrentWidget(win.dynamics)
        gt = win.dynamics.gating
        gt.parent().parent().setCurrentWidget(gt)
        gt.select("ryr1")                            # RyR1 bells
        if "never end" not in gt.text.text():
            raise RuntimeError(f"RyR1 gating drew: {gt.text.text()[:120]}")
        app.processEvents()
        win.grab().save(str(out / "gui_ryr_gating.png"))
    return False


def enter_transition(win, app, out) -> None:
    """Build the 8TKG -> 8TKF morph (the Transition tab's preset)."""
    tabs = win.transition.parentWidget().parentWidget()
    tabs.setCurrentWidget(win.transition)
    win.transition.preset.click()           # loads 8TKG, builds -> 8TKF


def enter_publication(win, app, out) -> None:
    """Open the Tree tab from P2.cyclostome_lineages' Show."""
    win.findings.tree.setCurrentItem(
        win.findings._items["P2.cyclostome_lineages"])
    win.findings.show_btn.click()           # opens the Tree tab
    if win.tabs.currentWidget() is not win.tree:
        raise RuntimeError("P2.cyclostome_lineages did not open the Tree tab")


def enter_dynamics(win, app, out) -> None:
    """Mak 1998's bell, then start a park/drive puff run."""
    win.tabs.setCurrentWidget(win.dynamics)
    gt = win.dynamics.gating
    gt.select("mak")                        # Mak et al. 1998
    if "Mak 1998" not in gt.text.text() or "K_inh" not in gt.note.text():
        raise RuntimeError(f"Mak gating drew: {gt.text.text()[:120]}")
    app.processEvents()
    win.grab().save(str(out / "gui_gating_mak.png"))
    pz = win.dynamics.puffs
    pz.parent().parent().setCurrentWidget(pz)
    pz.model.setCurrentIndex(pz.model.findData("park-drive"))
    sv.check_puff_rows(pz)
    if abs(pz.coupling.value() - 0.1) > 1e-9:
        raise RuntimeError(f"park/drive coupling {pz.coupling.value()}")
    pz.duration.setValue(5.0)
    pz.run()


def enter_unitary(win, app, out) -> None:
    """Start the ITPR3 state panel's unitary conductances."""
    win.tabs.setCurrentWidget(win.channel)
    win.channel.unitary_btn.click()


def enter_ryr(win, app, out) -> None:
    """Load RyR1 9HEO with no fill."""
    win.structure_panel.set_completeness("none")
    win.structure_panel.select("9HEO")     # a ryanodine receptor


def enter_sparks(win, app, out) -> None:
    """Start a RyR1 mean-field spark run (the spark steps check it)."""
    pz = win.dynamics.puffs
    pz.parent().parent().setCurrentWidget(pz)
    pz.model.setCurrentIndex(pz.model.findData("ryr1"))
    if pz.n.value() != 30 or pz.p.isEnabled():
        raise RuntimeError("the spark receptor did not take its cluster")
    pz.result = None
    pz.duration.setValue(10.0)
    pz.run()
