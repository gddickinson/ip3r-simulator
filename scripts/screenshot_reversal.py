"""The GUI smoke test's Round 7.24 steps (group ``reversal``): 8TKF's lumen
at the Ca²⁺ experiment's reversal under ``pb + csc``, coloured by Ca²⁺'s
concentration and then Cl⁻'s drop, equal to the headless reading; the
equilibrium colourings grey there; then (Round 7.26) Round 7.25's span well
drawn in its place, equal to the headless candidate, with the deposit's own
reading dashed beside it; then (Round 7.28) Round 7.27's Ca²⁺ site,
coloured by its occupancy and K⁺'s block, at the headless reading's V and
ions held (pinned in ``test_wall_candidates``);
then everything switched off.
"""

from __future__ import annotations

import numpy as np

from ip3r.render.colormaps import ramp
from ip3r.render.lumen_mesh import block_colors, conc_colors

__all__ = ["STEPS"]

_READING, _EXPERIMENT = "pb + csc", "Ca2+"
_CANDIDATE, _SITE = "span well", "Ca2+ site"


def _pick(combo, value) -> None:
    i = combo.findData(value)
    if i < 0:
        raise RuntimeError(f"no {value!r} in the lumen box")
    combo.setCurrentIndex(i)


def _reversal_start(win, app, out) -> None:
    box = win.channel.lumen_box
    win.structure_panel.set_completeness("none")
    win.tabs.setCurrentWidget(win.channel)
    _pick(box.experiment, _EXPERIMENT)
    _pick(box.colour, "conc:Ca2+")
    _pick(box.reversal_box, _READING)
    win.channel.show_lumen.setChecked(True)
    if win.scene.structure is None or win.scene.structure.name != "8TKF":
        win.structure_panel.select("8TKF")


def _reversal(win, app, out) -> bool:
    from ip3r.physics.lumen_reversal import reversal_lumen
    lc, sc = win.lumen, win.scene
    if sc.structure is None or sc.structure.name != "8TKF":
        return True                           # a refusal on 6DQN (shut) is right
    if lc.message.startswith("lumen not built"):
        raise RuntimeError(lc.message)
    if lc.busy:
        return True
    rev = lc.reversal
    if rev is None or rev.reading != _READING or rev.experiment != _EXPERIMENT:
        raise RuntimeError(f"8TKF's reversal not built: {lc.message}")
    batch = sc.scene.get("lumen")
    if batch is None or not batch.count or lc.field is not rev.lumen:
        raise RuntimeError("the reversal's own lumen is not the one drawn")
    ca = lc.mesh.sample(rev.conc["Ca2+"])
    if not np.allclose(lc.mesh.colors, conc_colors(ca)):
        raise RuntimeError("the surface is not coloured by Ca2+'s concentration")
    peak, _ = rev.peak("Ca2+")
    if not peak > rev.baths["Ca2+"][0]:
        raise RuntimeError(f"no Ca2+ well at reversal: peak {peak:.3g} M")
    rows = win.channel.canvas.axes
    labels = [ln.get_label() for ln in rows[2, 0].get_lines()]
    if not all(any(ion in t for t in labels) for ion in ("K+", "Cl-", "Ca2+")):
        raise RuntimeError(f"reversal drop plot drew {labels}")
    head = reversal_lumen(sc.structure, _READING, _EXPERIMENT, sc.summary)
    if abs(head.v - rev.v) > 1e-9 or f"{rev.v * 1e3:+.2f} mV" \
            not in win.channel.lumen_info.text():
        raise RuntimeError("the panel's reversal is not the headless one")
    app.processEvents()
    win.grab().save(str(out / "gui_lumen_reversal.png"))
    box = win.channel.lumen_box
    _pick(box.colour, "rdrop:Cl-")
    if not np.allclose(lc.mesh.colors, ramp(lc.mesh.sample(rev.drop["Cl-"]))):
        raise RuntimeError("the surface is not coloured by Cl-'s drop")
    _pick(box.colour, "wall")                 # no wall potential at reversal
    if not np.allclose(lc.mesh.colors, ramp(np.array([np.nan]))):
        raise RuntimeError("the wall colouring is not grey at reversal")
    _pick(box.colour, "conc:Ca2+")
    _pick(box.reversal_box, _CANDIDATE)
    return False


def _candidate(win, app, out) -> bool:
    from ip3r.physics.lumen_reversal import reversal_lumen
    lc, sc = win.lumen, win.scene
    if lc.message.startswith("lumen not built"):
        raise RuntimeError(lc.message)
    if lc.busy:
        return True
    rev = lc.reversal
    if rev is None or rev.reading != _CANDIDATE or not rev.candidate:
        raise RuntimeError(f"the candidate wall was not drawn: {lc.message}")
    if lc.own is None or lc.own.reading != _READING:
        raise RuntimeError("the deposit's own reading was not kept beside it")
    ca = lc.mesh.sample(rev.conc["Ca2+"])
    if not np.allclose(lc.mesh.colors, conc_colors(ca)):
        raise RuntimeError("the candidate is not coloured by Ca2+'s concentration")
    labels = [ln.get_label() for ln in win.channel.canvas.axes[1, 0].get_lines()]
    if not any("deposit (pb + csc)" in t for t in labels):
        raise RuntimeError(f"the deposit's reading is not beside it: {labels}")
    text = win.channel.lumen_info.text()
    head = reversal_lumen(sc.structure, _CANDIDATE, _EXPERIMENT, sc.summary)
    if abs(head.v - rev.v) > 1e-9 or "Beside the deposit's own wall" not in text:
        raise RuntimeError("the panel's candidate is not the headless one")
    app.processEvents()
    win.grab().save(str(out / "gui_lumen_candidate.png"))
    box = win.channel.lumen_box
    _pick(box.colour, "occupancy")            # grey: the well has no site
    if not np.allclose(lc.mesh.colors, ramp(np.array([np.nan]))):
        raise RuntimeError("the occupancy colouring is not grey without a site")
    _pick(box.reversal_box, _SITE)
    return False


def _site(win, app, out) -> bool:
    lc, sc = win.lumen, win.scene
    if lc.message.startswith("lumen not built"):
        raise RuntimeError(lc.message)
    if lc.busy:
        return True
    rev = lc.reversal
    if rev is None or rev.reading != _SITE or rev.occupancy is None:
        raise RuntimeError(f"the Ca2+ site was not drawn: {lc.message}")
    if not np.allclose(lc.mesh.colors, ramp(lc.mesh.sample(rev.occupancy))):
        raise RuntimeError("the site is not coloured by its occupancy")
    labels = [ln.get_label() for ln in win.channel.canvas.axes[2, 0].get_lines()]
    if not any("occupancy" in t for t in labels):
        raise RuntimeError(f"the site's occupancy is not plotted: {labels}")
    # the headless solve is six minutes; test_wall_candidates pins it at
    # +18.15 mV and 2.01 held (Round 7.27's crossing), so hold the panel to that
    if abs(rev.v * 1e3 - 18.15) > 0.05 or abs(rev.held - 2.01) > 0.02 \
            or "the site holds" not in win.channel.lumen_info.text():
        raise RuntimeError(f"the panel's site is not Round 7.27's: {lc.message}")
    app.processEvents()
    win.grab().save(str(out / "gui_lumen_site.png"))
    box = win.channel.lumen_box
    _pick(box.colour, "block")
    if not np.allclose(lc.mesh.colors, block_colors(lc.mesh.sample(rev.k_block))):
        raise RuntimeError("the site is not coloured by K+'s block")
    win.channel.show_lumen.setChecked(False)  # first: no re-solve behind us
    _pick(box.reversal_box, "equilibrium")
    _pick(box.colour, "drop")
    if lc.busy or sc.scene.get("lumen") is not None:
        raise RuntimeError("resetting the unticked lumen box started a solve")
    return False


STEPS = (_reversal_start, _reversal, _candidate, _site)
