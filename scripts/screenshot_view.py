"""The viewport and navigation checks of the GUI smoke test (Round 7.1).

Each function raises on a failure; ``screenshot_app.py`` calls them at the
step where their state exists. What they hold the GUI to:

- ``check_fit``: the loaded molecule reaches the edge of the view along one
  axis (NDC extent ≥ ``FILL``), and the viewport has a real share of the
  window (≥ ``VIEWPORT_SHARE``). Before Round 7.1 the docks left it a fifth.
  The Structure panel fits its dock (no sideways scroll hiding a subunit).
- ``check_subunit_refit``: hiding three subunits re-fits to the one left.
- ``check_site_view``: a shell check's "Show" centres the camera on an IP3
  site and the pocket fills the view.
- ``check_list``: the deposition list is grouped by family, RyR1 collapsed.
- ``check_puff_rows``: the RyR-only Puffs controls are hidden for IP3R.
"""

from __future__ import annotations

import numpy as np

__all__ = ["check_fit", "check_subunit_refit", "check_site_view", "check_list",
           "check_puff_rows"]

#: The molecule must reach at least this NDC extent along one axis.
FILL = 0.85
#: The IP3 pocket must reach this: the site view pulls back SITE_MARGIN.
SITE_FILL = 0.7
#: The viewport must be at least this share of the window's width.
VIEWPORT_SHARE = 0.3


def _extent(win, atoms=None) -> float:
    xyz = win.scene.view.structure.xyz
    xyz = xyz if atoms is None else xyz[atoms]
    return max(win.scene.scene.camera.screen_extent(xyz))


def check_fit(win, app) -> str:
    app.processEvents()
    share = win.viewport.width() / win.width()
    if share < VIEWPORT_SHARE:
        raise RuntimeError(f"viewport only {share:.0%} of the window")
    sp = win.structure_panel                 # no sideways scroll hides a control
    if sp.width() > sp.parentWidget().width():
        raise RuntimeError(f"Structure panel {sp.width()} px in a {sp.parentWidget().width()} px dock")
    e = _extent(win)
    if not FILL <= e <= 1.0:
        raise RuntimeError(f"molecule spans {e:.2f} of the view, not {FILL}-1")
    return f"viewport {share:.0%} of the window, molecule extent {e:.2f}"


def check_subunit_refit(win, app) -> None:
    boxes = win.structure_panel.chain_boxes
    keep = next(iter(boxes))
    for c, b in boxes.items():
        b.setChecked(c == keep)
    app.processEvents()
    st = win.scene.structure
    e = _extent(win, np.flatnonzero(st.chain == keep))
    for b in boxes.values():
        b.setChecked(True)
    app.processEvents()
    if not FILL <= e <= 1.0:
        raise RuntimeError(f"one subunit spans {e:.2f} of the view after the refit")


def check_site_view(win, app) -> None:
    app.processEvents()
    if win.scene.fit_target != "site":
        raise RuntimeError(f"shell check framed {win.scene.fit_target!r}, not an IP3 site")
    atoms, _ = win.scene._site_atoms()
    e = _extent(win, atoms)
    if not SITE_FILL <= e <= 1.0:
        raise RuntimeError(f"IP3 pocket spans {e:.2f} of the view")
    cam = win.scene.scene.camera
    if cam.slab_front is None:
        raise RuntimeError("site view does not clip what lies in front of the pocket")


def check_list(win) -> None:
    lst = win.structure_panel.list
    heads = {lst.topLevelItem(i).text(0).split(" (")[0]: lst.topLevelItem(i)
             for i in range(lst.topLevelItemCount())}
    ryr = [h for k, h in heads.items() if k.startswith("RyR1")]
    if len(heads) != 2 or len(ryr) != 1:
        raise RuntimeError(f"deposition list headings {list(heads)}")
    if ryr[0].isExpanded():
        raise RuntimeError("RyR1 group open while an IP3R deposit is shown")


def check_puff_rows(pz) -> None:
    ip3r = not pz.form.isRowVisible(pz.mg) and not pz.form.isRowVisible(pz.reading) \
        and pz.mg_scan_btn.isHidden() and pz.form.isRowVisible(pz.p)
    if not ip3r:
        raise RuntimeError("RyR-only Puffs controls shown for an IP3R receptor")
