"""The Round 7.3 steps of the GUI smoke test (``screenshot_app.py``): the
three IP3R gating models side by side, the A-subspace headline of the
Transition tab, and the park/drive cluster in its microdomain; Round
7.7's reported tree beside the ``--bnni`` re-search; and Round 7.8's
publication views (the Genomes lesion layer opened from its check, a Range
clade's S23 genomes, the VUS thresholds' intervals); Round 7.9's rat fill;
Round 7.10's lumen and where the voltage falls; Round 7.12's charged lumen;
Round 7.14's dielectric closure in it, and Round 7.16's image cost on that.

``ip3r_step`` follows the spark steps' contract: True means a worker is
still running (call again), and a wrong result raises.
"""

from __future__ import annotations

import numpy as np

from ip3r.render.colormaps import MISSING, ramp
from ip3r.render.lumen_mesh import image_colors, wall_colors

__all__ = ["IP3R_STEPS", "ip3r_step", "check_transition_headline"]


def check_transition_headline(win) -> None:
    """The first plot is the collective A modes together, titled with the
    value the report prints."""
    tp = win.transition
    sub = tp.result.overlap.subspace("A")
    title = tp.canvas.axes[0, 0].get_title()
    if not (len(sub.index) and f"{sub.total:.2f}" in title and "A subspace" in title):
        raise RuntimeError(f"transition headline: '{title}'")
    if f"{sub.total:.3f}" not in tp.report.text():
        raise RuntimeError("the headline and the report disagree")
    for i in reversed(range(tp.lower.count())):      # ends on the gate
        tp.lower.setCurrentIndex(i)
        lower = tp.canvas.axes[1, 0].get_title()
        if tp.canvas.axes[0, 0].get_title() != title or not lower:
            raise RuntimeError(f"lower plot {tp.lower.currentData()}: '{lower}'")
    if tp._marker is None:
        raise RuntimeError("the gate plot lost its frame marker")


def _gating(win, app, out) -> None:
    gt = win.dynamics.gating
    win.tabs.setCurrentWidget(win.dynamics)
    gt.parent().parent().setCurrentWidget(gt)
    gt.select("pd")
    if "Park/drive" not in gt.text.text() or "parked" not in \
            gt.canvas.axes[0, 1].get_ylabel():
        raise RuntimeError(f"park/drive gating drew: {gt.text.text()[:120]}")
    gt.select("compare")
    titles = [ax.get_title() for ax in gt.canvas.axes[0]]
    if len(gt.canvas.axes[0, 0].get_lines()) < 3 or not all("IP3" in t for t in titles):
        raise RuntimeError(f"gating comparison drew {titles}")
    app.processEvents()
    win.grab().save(str(out / "gui_gating_compare.png"))
    pz = win.dynamics.puffs
    pz.parent().parent().setCurrentWidget(pz)
    pz.model.setCurrentIndex(pz.model.findData("ryr1"))
    if pz.domain.isVisibleTo(pz):
        raise RuntimeError("microdomain controls shown for RyR1")
    pz.model.setCurrentIndex(pz.model.findData("park-drive"))
    if not pz.domain.isVisibleTo(pz):
        raise RuntimeError("microdomain controls hidden for park/drive")
    pz.result = None
    pz.domain.duration.setValue(30.0)
    pz.run_domain()


def _domain(win, app, out) -> bool:
    pz = win.dynamics.puffs
    if pz.result is None:
        if pz.text.text().startswith("Simulating"):
            return True
        raise RuntimeError(pz.text.text())
    s = pz.result["stats"]
    if not (s.n_puffs >= 5 and 0.4 < s.blip_mean < 1.0 and len(s.ipis) >= 5
            and s.lam > 0 and pz.result["open_fraction"] < 0.2):
        raise RuntimeError(f"microdomain drew {s}")
    if "Thurley" not in pz.text.text():
        raise RuntimeError(pz.text.text()[:160])
    app.processEvents()
    win.grab().save(str(out / "gui_puffs_domain.png"))
    return False


def _tree_pair_start(win, app, out) -> None:
    win.tabs.setCurrentWidget(win.tree)
    win.tree.beside.setChecked(True)


def _tree_pair(win, app, out) -> bool:
    tr = win.tree
    if tr.pair is None or tr.root is None:
        if tr.status.text().startswith("Reading"):
            return True
        raise RuntimeError(tr.status.text())
    titles = [ax.get_title() for ax in tr.canvas.axes[0]]
    if titles != ["reported", "--bnni"] or "9 of 10 clade claims held" not in tr.status.text():
        raise RuntimeError(f"tree pair drew {titles}: {tr.status.text()[:120]}")
    app.processEvents()
    win.grab().save(str(out / "gui_tree_bnni.png"))
    tr.beside.setChecked(False)
    if tr.canvas.axes.shape != (1, 1):
        raise RuntimeError("unchecking --bnni left two trees drawn")
    return False


def _lesion_start(win, app, out) -> None:
    win.findings.tree.setCurrentItem(win.findings._items["P3.lesion_strata"])
    win.findings.show_btn.click()                   # opens the Genomes tab on Aves
    if win.tabs.currentWidget() is not win.genomes:
        raise RuntimeError("P3.lesion_strata did not open the Genomes tab")


def _lesion(win, app, out) -> bool:
    g = win.genomes
    if g.grid is None or g.layer.currentData() != "lesion":
        if g.grid is None and g.status.text().startswith("Reading"):
            return True
        raise RuntimeError(f"lesion layer not shown: {g.status.text()[:120]}")
    itpr3 = g.shown.counts("lesion", ["ITPR3"])
    if g.vclass.currentText() != "Aves" or (itpr3.get("excess"), itpr3.get("deficit")) != (25, 2):
        raise RuntimeError(f"lesion layer on {g.vclass.currentText()}: ITPR3 {itpr3}")
    if g.info["bar_row"] is None:
        raise RuntimeError("the bird rows are in N50 order but no bar was drawn")
    app.processEvents()
    win.grab().save(str(out / "gui_genomes_lesion.png"))
    g.vclass.setCurrentIndex(0)
    return False


def _range_genomes(win, app, out) -> None:
    from ip3r.analysis import range_genomes as RG
    r = win.range
    win.tabs.setCurrentWidget(r)
    clade = max((c for c, _ in RG.clades_with_genomes(r.genomes)),     # one with genes
                key=lambda c: sum(g.copies for g in RG.in_clade(r.genomes, c)))
    if not r.show_genomes(clade):
        raise RuntimeError(f"{clade} has no genomes to show")
    want = RG.in_clade(r.genomes, clade)
    if r.info.get("genomes") != len(want) or r.min_n.isEnabled() or not r.clade.isEnabled():
        raise RuntimeError(f"range genomes drew {r.info} for {clade}")
    app.processEvents()
    win.grab().save(str(out / "gui_range_genomes.png"))
    r.view.setCurrentIndex(0)
    if r.shown_genomes or not r.rows:
        raise RuntimeError("back to clades left the genome rows")


def _vus_bands(win, app, out) -> None:
    v = win.variants
    win.tabs.setCurrentWidget(v)
    v.bucket.setCurrentText("VUS")
    v.layer.setCurrentIndex(v.layer.findData("deep"))
    v.paralog.setCurrentText("ITPR1")
    if "Unbounded" in v.status.text() or not v.bands.bounded:
        raise RuntimeError(f"ITPR1's thresholds: {v.status.text()[-160:]}")
    v.paralog.setCurrentText("ITPR2")
    cells = {v.table.item(i, 6).text() for i in range(v.table.rowCount())}
    if "Unbounded: P/LP (n=1)" not in v.status.text() or \
            not any("near P/LP" in c for c in cells):
        raise RuntimeError(f"ITPR2's thresholds: {v.status.text()[-200:]}")
    app.processEvents()
    win.grab().save(str(out / "gui_variants_bands.png"))


def _rat_fill_start(win, app, out) -> None:
    win.structure_panel.set_completeness("gaps")
    win.structure_panel.select("7LHF")


def _rat_fill(win, app, out) -> bool:
    """Round 7.9: rat 7LHF filled from rat isoform 8 through an alignment."""
    fc, sc = win.fills, win.scene
    if sc.structure is None or sc.structure.name != "7LHF" or fc.model is None:
        if fc.message.startswith("not filled"):
            raise RuntimeError(fc.message)
        return True
    m = fc.model
    if m.deposit != "7LHF":
        return True
    if m.prediction != "AF-P29994-8-F1" or m.numbering.route != "alignment" or \
            m.numbering.unmapped != ((318, 332), (1693, 1732)):
        raise RuntimeError(f"7LHF fill: {m.summary()}")
    if sc.fill.view is None or "the model lacks" not in win.structure_panel.fill_info.text():
        raise RuntimeError("7LHF's fill not drawn, or its splice segments not named")
    app.processEvents()
    win.grab().save(str(out / "gui_alphafold_rat.png"))
    return False


def _lumen_start(win, app, out) -> None:
    win.structure_panel.set_completeness("none")
    win.tabs.setCurrentWidget(win.channel)
    win.channel.show_lumen.setChecked(True)
    win.structure_panel.select("8TKF")


def _lumen(win, app, out) -> bool:
    """Round 7.10: 8TKF's lumen drawn, coloured by the potential, and the
    Channel panel's plot of where the voltage falls; hidden off the deposit."""
    from ip3r.physics.lumen_field import lumen_field
    lc, sc = win.lumen, win.scene
    if sc.structure is None or sc.structure.name != "8TKF" or lc.field is None:
        if lc.message.startswith("lumen not built"):
            raise RuntimeError(lc.message)
        return True
    f, batch = lc.field, sc.scene.get("lumen")
    if f.name != "8TKF" or not f.conducts or batch is None or not batch.count \
            or not batch.visible:
        raise RuntimeError(f"8TKF lumen not drawn: {lc.message}")
    if not np.ptp(lc.mesh.phi) > 0.5:
        raise RuntimeError("the lumen's colour does not span the drop")
    labels = [ln.get_label() for ln in win.channel.canvas.axes[1, 0].get_lines()]
    if not any("3-D" in t for t in labels) or not any("1-D" in t for t in labels):
        raise RuntimeError(f"lumen plot drew {labels}")
    head = lumen_field(sc.structure, sc.summary)      # the panel = the headless field
    if f"{head.half_z('3d'):+.1f}" not in win.channel.lumen_info.text() \
            or "filter" not in win.channel.lumen_info.text():
        raise RuntimeError(f"lumen text: {win.channel.lumen_info.text()[:160]}")
    app.processEvents()
    win.grab().save(str(out / "gui_lumen.png"))
    xyz = sc.structure.xyz
    sc.move_overlays(xyz + 1.0)                          # a morph or mode frame
    if batch.visible:
        raise RuntimeError("the deposit's lumen stayed up on a moved frame")
    sc.move_overlays(xyz)
    if not batch.visible:
        raise RuntimeError("the lumen did not come back with the deposit")
    box = win.channel.lumen_box                          # on to Round 7.12
    box.colour.setCurrentIndex(box.colour.findData("wall"))
    if not np.allclose(lc.mesh.colors, wall_colors(np.zeros(1))):
        raise RuntimeError("the neutral pore's wall potential is not zero")
    box.pairs.setChecked(True)
    box.charge.setCurrentIndex(box.charge.findData("pb"))
    return False


def _lumen_charged(win, app, out) -> bool:
    """Round 7.12: 8TKF's paired wall under PB, the surface coloured by its
    potential, the panel's reading equal to the headless one."""
    from ip3r.physics.lumen_charge import charged_lumen
    lc, sc = win.lumen, win.scene
    if lc.busy:
        return True
    c = lc.charged
    if c is None or c.closure != "pb" or not c.pair_bridges:
        raise RuntimeError(f"8TKF's charged lumen not built: {lc.message}")
    if not c.ratio < 1.0:
        raise RuntimeError(f"paired PB should lower K+ g: x{c.ratio:.2f}")
    if not np.allclose(lc.mesh.colors, wall_colors(lc.mesh.sample(c.u))):
        raise RuntimeError("the surface is not coloured by the wall potential")
    if win.channel.canvas.axes.shape[0] != 3 or \
            f"×{c.ratio:.2f}" not in win.channel.lumen_info.text():
        raise RuntimeError(f"charged lumen text: {win.channel.lumen_info.text()[:200]}")
    head = charged_lumen(sc.structure, lc.field, "pb", sc.summary, True)
    if abs(head.g - c.g) > 1e-9 * abs(c.g):
        raise RuntimeError("the panel's charged reading is not the headless one")
    win.channel.show_pore.setChecked(False)              # the spheres hide the lumen
    app.processEvents()
    win.grab().save(str(out / "gui_lumen_charged.png"))
    box = win.channel.lumen_box
    box.colour.setCurrentIndex(box.colour.findData("drop"))
    if not np.allclose(lc.mesh.colors, ramp(lc.mesh.sample(c.mu))):
        raise RuntimeError("the charged drop colouring is not the K+ drop")
    box.pairs.setChecked(False)                          # on to Round 7.14
    box.charge.setCurrentIndex(box.charge.findData("dielectric"))
    return False


def _lumen_dielectric(win, app, out) -> bool:
    """Round 7.14: 8TKF's wall under the dielectric closure, the salt bridge
    a dipole, equal to the headless reading; then "+ image" ticked."""
    from ip3r.physics.lumen_charge import charged_lumen
    lc, sc = win.lumen, win.scene
    if lc.busy:
        return True
    c = lc.charged
    if c is None or c.closure != "dielectric" or c.pair_bridges:
        raise RuntimeError(f"8TKF's dielectric lumen not built: {lc.message}")
    if not c.ratio > 1.0 or "in the box" not in win.channel.lumen_info.text():
        raise RuntimeError(f"dielectric text: {win.channel.lumen_info.text()[:200]}")
    head = charged_lumen(sc.structure, lc.field, "dielectric", sc.summary)
    if abs(head.g - c.g) > 1e-9 * abs(c.g):
        raise RuntimeError("the panel's dielectric reading is not the headless one")
    app.processEvents()
    win.grab().save(str(out / "gui_lumen_dielectric.png"))
    box = win.channel.lumen_box                          # on to Round 7.16
    if not box.image_box.isEnabled():
        raise RuntimeError("'+ image' is not offered under the dielectric closure")
    box.colour.setCurrentIndex(box.colour.findData("image"))
    if not np.allclose(lc.mesh.colors, MISSING):
        raise RuntimeError("the image colouring is not grey before W is solved")
    box.image_box.setChecked(True)
    return False


def _lumen_image(win, app, out) -> bool:
    """Round 7.16: 8TKF's dielectric wall with the image cost (Round 7.15's
    cached W, the lumen re-cut at 1 Å), the surface coloured by W, equal to
    the headless reading; then everything switched off."""
    from ip3r.parameters import PARAMETERS as _P
    from ip3r.physics.lumen_charge import charged_lumen
    lc, sc = win.lumen, win.scene
    if lc.busy:
        return True
    c = lc.charged
    if c is None or not c.image:
        raise RuntimeError(f"8TKF's image reading not built: {lc.message}")
    if lc.field.volume.spacing != _P.value("born.lumen_spacing"):
        raise RuntimeError("the image reading is not on its grid")
    if not np.allclose(lc.mesh.colors, image_colors(lc.mesh.sample(c.w))):
        raise RuntimeError("the surface is not coloured by the image cost")
    labels = [ln.get_label() for ln in win.channel.canvas.axes[2, 0].get_lines()]
    if not any("image cost" in t for t in labels) or \
            "+ image" not in win.channel.lumen_info.text():
        raise RuntimeError(f"image plot / text: {labels}")
    head = charged_lumen(sc.structure, lc.field, "dielectric", sc.summary,
                         image=True)
    if abs(head.g - c.g) > 1e-9 * abs(c.g):
        raise RuntimeError("the panel's image reading is not the headless one")
    app.processEvents()
    win.grab().save(str(out / "gui_lumen_image.png"))
    box = win.channel.lumen_box
    win.channel.show_lumen.setChecked(False)     # first: no re-solve behind us
    if sc.scene.get("lumen") is not None:
        raise RuntimeError("unticking left the lumen drawn")
    box.image_box.setChecked(False)
    for w, v in ((box.charge, "none"), (box.colour, "drop")):
        w.setCurrentIndex(w.findData(v))
    if lc.busy:
        raise RuntimeError("resetting the unticked lumen box started a solve")
    return False


_STEPS = (_gating, _domain, _tree_pair_start, _tree_pair, _lesion_start, _lesion,
          _range_genomes, _vus_bands, _rat_fill_start, _rat_fill, _lumen_start, _lumen,
          _lumen_charged, _lumen_dielectric, _lumen_image)
IP3R_STEPS = len(_STEPS)


def ip3r_step(k: int, win, app, out) -> bool:
    """Run step ``k``; True means "not ready yet, call again"."""
    return bool(_STEPS[k](win, app, out))
