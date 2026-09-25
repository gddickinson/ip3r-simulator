"""The Round 7.3 steps of the GUI smoke test (``screenshot_app.py``): the
three IP3R gating models side by side, the A-subspace headline of the
Transition tab, and the park/drive cluster in its microdomain; Round
7.7's reported tree beside the ``--bnni`` re-search; and Round 7.8's
publication views (the Genomes lesion layer opened from its check, a Range
clade's S23 genomes, the VUS thresholds' intervals).

``ip3r_step`` follows the spark steps' contract: True means a worker is
still running (call again), and a wrong result raises.
"""

from __future__ import annotations

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


_STEPS = (_gating, _domain, _tree_pair_start, _tree_pair, _lesion_start, _lesion,
          _range_genomes, _vus_bands)
IP3R_STEPS = len(_STEPS)


def ip3r_step(k: int, win, app, out) -> bool:
    """Run step ``k``; True means "not ready yet, call again"."""
    return bool(_STEPS[k](win, app, out))
