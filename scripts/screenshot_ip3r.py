"""The Round 7.3 steps of the GUI smoke test (``screenshot_app.py``): the
three IP3R gating models side by side, the A-subspace headline of the
Transition tab, and the park/drive cluster in its microdomain.

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


_STEPS = (_gating, _domain)
IP3R_STEPS = len(_STEPS)


def ip3r_step(k: int, win, app, out) -> bool:
    """Run step ``k``; True means "not ready yet, call again"."""
    return bool(_STEPS[k](win, app, out))
