"""The smoke test's Round 7.5 steps: panel views in sessions, controls that
follow the parameters.

``set_panel_view`` puts non-default Dynamics and Variants controls on the
window before the session is saved, and ``disturb_panel_view`` moves them
after, so the whole-session comparison in ``screenshot_app`` proves they
come back. ``mode_round_trip`` saves a view with a normal mode animating,
stops it, and restores; it has to wait for the ANM, so it returns True
while waiting. ``check_follow`` edits the registered values that seed the
Puffs panel and the gating plot, and checks what follows and what does not.
"""

from __future__ import annotations

__all__ = ["set_panel_view", "disturb_panel_view", "mode_round_trip", "check_follow"]


def set_panel_view(win) -> None:
    dy, v = win.dynamics, win.variants
    dy.gating.select("compare")
    pz = dy.puffs
    pz.model.setCurrentIndex(pz.model.findData("park-drive"))
    pz.n.setValue(12)                              # typed: not the registry's
    pz.domain.duration.setValue(90)
    pz.domain.clamp.setCurrentIndex(pz.domain.clamp.findData("none"))
    dy.tabs.setCurrentWidget(pz)
    v.bucket.setCurrentText("VUS")
    v.layer.setCurrentIndex(v.layer.findData("deep"))
    v.draw.setChecked(True)


def disturb_panel_view(win) -> None:
    dy, v = win.dynamics, win.variants
    dy.gating.select("dyk")
    dy.puffs.model.setCurrentIndex(dy.puffs.model.findData("dyk"))
    dy.puffs.domain.duration.setValue(30)
    dy.tabs.setCurrentWidget(dy.gating)
    v.bucket.setCurrentText("all")
    v.layer.setCurrentIndex(0)
    v.draw.setChecked(False)


def mode_round_trip(win, state: dict) -> bool:
    """True while waiting (the ANM, then the restore's transition rebuild)."""
    ss = win.sessions
    if win.modes.modes is None:
        if not state.get("modes_started"):
            state["modes_started"] = True
            win.compute_modes()
        return True
    if "mode_saved" not in state:
        index = win.modes.modes.first("A") or 0
        win.modes.amp.setValue(20)
        win.modes.table.selectRow(index)
        if win.scene.animated_mode != (index, 20.0):
            raise RuntimeError(f"animating {win.scene.animated_mode}, not mode {index}")
        saved = state["mode_saved"] = ss.capture()
        if saved.modes != {"index": index, "amplitude": 20.0}:
            raise RuntimeError(f"session holds modes {saved.modes}")
        win.scene.stop_animation()
        win.modes.amp.setValue(5)
        if not ss.apply(saved, parameters="keep"):
            raise RuntimeError("the mode session was not started")
        return True
    if ss.pending is not None or ss._frame is not None or ss._mode is not None:
        return True
    want = state["mode_saved"].modes
    if win.scene.animated_mode != (want["index"], want["amplitude"]):
        raise RuntimeError(f"restored animation {win.scene.animated_mode}, "
                           f"saved {want}")
    win.transition.slider.setValue(5)          # a frame stops the animation
    if ss.capture().modes:
        raise RuntimeError("a stopped animation is still in the session")
    return False


def check_follow(win) -> None:
    """Seeded spins follow an edit, a typed one does not; the gating plot is
    redrawn; a reset puts the seeded values back."""
    from ip3r.parameters import PARAMETERS
    pz, gt = win.dynamics.puffs, win.dynamics.gating
    pz.model.setCurrentIndex(pz.model.findData("park-drive"))
    seeded = pz.coupling.value()
    pz.n.setValue(12)                                     # typed over
    typed = pz.domain.duration.value()           # 90 s, restored from the session
    domain = PARAMETERS.value("domain.gui_duration")
    if typed == domain:
        raise RuntimeError("the restored duration is the seeded one: no test")
    PARAMETERS.set_value("domain.gui_duration", domain + 30)
    if pz.domain.duration.value() != typed:
        raise RuntimeError("a restored (typed) duration followed an edit")
    PARAMETERS.reset()
    pz.domain.duration.setValue(domain)           # back on the seeded value
    PARAMETERS.set_value("puff.pd_ca_per_open", 2 * seeded)
    PARAMETERS.set_value("puff.n_channels", 25.0)
    PARAMETERS.set_value("domain.gui_duration", domain + 30)
    try:
        if abs(pz.coupling.value() - 2 * seeded) > 1e-9:
            raise RuntimeError(f"coupling {pz.coupling.value()} did not follow")
        if pz.n.value() != 12:
            raise RuntimeError("a typed cluster size was overwritten by an edit")
        if pz.domain.duration.value() != domain + 30:
            raise RuntimeError("the microdomain duration did not follow")
        win.tabs.setCurrentWidget(win.dynamics)
        win.dynamics.tabs.setCurrentWidget(gt)
        gt.select("compare")
        before = gt.text.text()
        PARAMETERS.set_value("mak.compare_ip3_high", 3.0)
        if gt.text.text() == before or "→ 3 µM" not in gt.text.text():
            raise RuntimeError(f"the gating plot did not follow: {gt.text.text()[:100]}")
    finally:
        PARAMETERS.reset()
    if abs(pz.coupling.value() - seeded) > 1e-9 or pz.n.value() != 12:
        raise RuntimeError("reset did not put the seeded coupling back")
