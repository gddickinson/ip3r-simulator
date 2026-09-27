"""The right-click menu on the viewport (Round 7.18, after PIEZO1's).

Two rules, both PIEZO1's:

* **Nothing is implemented twice.** Style and colour entries set the
  Structure panel's own combo boxes, so the panel never says "cartoon"
  while the model shows spheres; views call the same scene methods the
  View menu does; selection is :class:`~ip3r.ui.selection.SelectionController`.
* **Opening the menu is not a selection.** The residue under the cursor is
  named in a disabled header; nothing changes unless an entry is chosen.

With an atom under the cursor the menu gains that residue's entries;
without one it keeps what acts on the view, so it is never empty.
"""

from __future__ import annotations

from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QMenu

__all__ = ["build_context_menu"]


def _add(menu: QMenu, text: str, slot, tip: str = "", checked: bool | None = None) -> QAction:
    a = QAction(text, menu)
    if tip:
        a.setToolTip(tip)
    if checked is not None:
        a.setCheckable(True)
        a.setChecked(checked)
    a.triggered.connect(slot)
    menu.addAction(a)
    return a


def build_context_menu(win, index: int) -> QMenu:
    """The menu for a right-click on atom ``index`` (-1: empty space)."""
    menu = QMenu(win)
    menu.setToolTipsVisible(True)
    st = win.scene.structure
    if st is not None and 0 <= index < st.n_atoms:
        _residue_entries(win, menu, index)
    _selection_entries(win, menu)
    _appearance_entries(win, menu)
    _view_entries(win, menu)
    return menu


def _residue_entries(win, menu: QMenu, index: int) -> None:
    st, sel = win.scene.structure, win.selection
    ch, r, name = str(st.chain[index]), int(st.res_seq[index]), str(st.res_name[index])
    header = QAction(win.scene.describe_atom(index) or f"{ch}:{name}{r}", menu)
    header.setEnabled(False)
    menu.addAction(header)
    menu.addSeparator()
    say = win.statusBar().showMessage
    if not st.hetero[index]:
        _add(menu, f"Select {name}{r} (chain {ch})", lambda: sel.select({(ch, r)}))
        _add(menu, f"Add / remove {name}{r} (chain {ch})",
             lambda: sel.select({(ch, r)} ^ sel.residues))
        _add(menu, f"Select {r} on every subunit",
             lambda: say(f"{r} selected on {sel.select_everywhere(r)} subunits"),
             "The same residue number on each chain that resolves it (C4 copies).")
        _add(menu, f"Select chain {ch}", lambda: sel.select_chain(ch))
        _add(menu, f"Show {r} in the sequence window",
             lambda: win.show_sequence(chain=ch, residue=r))
    _add(menu, "Centre here", lambda: sel.centre_on(index),
         "Pivot about this atom and zoom to 40 Å; the automatic fit stops "
         "(Ctrl+0 resumes it).")
    _add(menu, "Measure from this atom",
         lambda: say(sel.add_point(index)),
         "The next atom added (here, or by a click while measuring) closes a "
         "distance.")
    menu.addSeparator()


def _selection_entries(win, menu: QMenu) -> None:
    sel = win.selection
    has = bool(sel.residues)
    _add(menu, "Centre on selection", sel.centre_on_selection).setEnabled(has)
    _add(menu, "Clear selection", sel.clear).setEnabled(has)
    _add(menu, "Measure distances (click two atoms)", sel.arm,
         checked=sel.measuring)
    _add(menu, "Clear distances", sel.clear_distances).setEnabled(bool(sel.distances))
    menu.addSeparator()


def _combo_menu(menu: QMenu, title: str, combo) -> None:
    sub = menu.addMenu(title)
    for i in range(combo.count()):
        a = QAction(combo.itemText(i), sub)
        a.setCheckable(True)
        a.setChecked(i == combo.currentIndex())
        a.triggered.connect(lambda _=False, i=i: combo.setCurrentIndex(i))
        sub.addAction(a)


def _appearance_entries(win, menu: QMenu) -> None:
    sp = win.structure_panel
    _combo_menu(menu, "Style", sp.style)
    _combo_menu(menu, "Colour by", sp.color)
    menu.addSeparator()


def _view_entries(win, menu: QMenu) -> None:
    sc = win.scene
    _add(menu, "Side view (cytosol up)", sc.side_view)
    _add(menu, "Top view (down the pore)", sc.top_view)
    _add(menu, "IP3 site", win._site_view)
    _add(menu, "Fit to view", sc.fit_view)
