"""The application window: viewport in the middle, panels around it.

Everything slow — parsing a tetramer, measuring its channel, solving for
normal modes, running the checks — happens on a worker thread
(:mod:`ip3r.ui.workers`); this module only moves results onto the screen.
Drawing extras (pore, highlights, mode animation) is in
:mod:`ip3r.ui.scene_controller` so this file stays about layout and wiring.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QKeySequence
from PyQt6.QtWidgets import (QComboBox, QDockWidget, QFileDialog, QMainWindow, QMessageBox,
                             QScrollArea, QTabWidget, QToolBar)

from .. import __version__
from ..config import DEFAULT_STRUCTURE, SETTINGS, genes_results
from ..io import loader
from ..parameters import PARAMETERS
from ..io.fetch import fetch_all
from ..io.registry import get_entry
from ..structure.channel import measure_channel
from .channel_panel import ChannelPanel
from .dynamics_panel import DynamicsPanel
from .fill_controller import FillController
from .findings_panel import FindingsPanel
from .genomes_panel import GenomesPanel
from .gl_widget import ViewportWidget
from .modes_panel import ModesPanel
from .params_banner import ParametersBanner
from .params_dialog import ParametersDialog
from .range_panel import RangePanel
from ..core.modules import MODULES_KEY
from .scene_controller import SceneController
from .session_controller import SessionController
from .structure_panel import StructurePanel
from .transition_controller import TransitionController, build_transition
from .transition_panel import TransitionPanel
from .tree_panel import TreePanel
from .variants_panel import VariantsPanel
from .workers import run_async

__all__ = ["MainWindow"]

#: Which measured sites a structural check highlights when shown.
CHECK_SITES = {"S0.ip3_contacts": "ip3_contact", "S0.selectivity_filter": "filter_lining",
               "S0.gate": "gate_lining", "P6.shell_agreement": "ip3_contact",
               "P6.contacts_heavy_atom": "ip3_contact", "P6.module_map": MODULES_KEY,
               "P6.module_contrast": MODULES_KEY, "P6.loop_reverses": MODULES_KEY}
#: Checks whose "Show on structure" is a colouring rather than a site set.
#: Checks whose "Show" opens the Tree tab rather than the structure.
CHECK_TREE = frozenset(("P2.sister_pair", "P2.paralog_clades", "P2.cyclostome_lineages",
                        "P2.support_bar"))
#: Checks whose "Show" opens the Genomes tab, and the layer it shows.
CHECK_GENOMES = {"P3.no_absent_cells": "state", "P3.false_negatives": "miss",
                 "P3.miss_by_contiguity": "miss", "P3.contiguity_tests": "miss",
                 "P4.unreachable": "recovery", "P4.recovery_channels": "recovery"}
#: Checks whose "Show" opens the Range tab (Paper 1).
CHECK_RANGE = frozenset(("P1.presence_range", "P1.kingdom_absences", "P1.relaxed_controls",
                         "P1.absence_targets", "P1.absences", "P1.copy_number",
                         "P1.record_chase"))
CHECK_COLOURS = {k: "ligand_shell" for k in ("P6.shell_distances", "P6.shell_constraint",
                                             "P6.shell_trend", "P6.no_contact_step")}
#: Checks about one IP3 pocket: "Show" centres the camera on a site.
CHECK_SITE_VIEW = frozenset(CHECK_COLOURS) | frozenset(
    k for k, v in CHECK_SITES.items() if v == "ip3_contact")
#: Dock widths on first show (px); the viewport takes the rest. Without
#: them Qt sizes the docks from their contents and leaves the molecule a
#: strip a fifth of the window wide.
DOCK_WIDTHS = (380, 580)


def compact_combos(root, chars: int = 10, wide: int = 160) -> None:
    """Let every wide combo box under ``root`` shrink below its longest item.

    A combo's default minimum is its longest entry, and one row of long
    layer titles held the whole Analysis dock at 850 px.
    """
    for combo in root.findChildren(QComboBox):
        if combo.minimumSizeHint().width() <= wide:
            continue
        combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        combo.setMinimumContentsLength(chars)
        combo.setMinimumWidth(wide)       # the policy alone leaves minimumSizeHint wide
        combo.view().setMinimumWidth(combo.view().sizeHintForColumn(0) + 24)  # popup: full text


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"IP3R Structural Simulator {__version__}")
        self.resize(1560, 960)
        self.viewport = ViewportWidget(SETTINGS.render)
        self.setCentralWidget(self.viewport)
        # Across the whole window, above the docks: the centre is too narrow.
        self.params_banner = ParametersBanner()
        self.params_banner.edit_requested.connect(self.edit_parameters)
        self.params_strip = strip = QToolBar("Parameters modified")
        strip.setMovable(False)
        strip.setFloatable(False)
        strip.toggleViewAction().setVisible(False)
        strip.addWidget(self.params_banner)
        self.params_banner.visibility_changed.connect(strip.setVisible)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, strip)
        strip.setVisible(PARAMETERS.modified)
        self.scene = SceneController(self.viewport)
        self._pending_check: str | None = None
        self._pending_transition: tuple | None = None
        self.morph = TransitionController(self.scene)
        self.sessions = SessionController(self)

        self.structure_panel = StructurePanel()
        compact_combos(self.structure_panel)
        self._dock("Structure", self.structure_panel, Qt.DockWidgetArea.LeftDockWidgetArea)
        self.channel = ChannelPanel()
        self.modes = ModesPanel()
        self.transition = TransitionPanel()
        self.dynamics = DynamicsPanel()
        self.findings = FindingsPanel()
        self.variants = VariantsPanel()
        self.tree = TreePanel()
        self.genomes = GenomesPanel()
        self.range = RangePanel()
        self.tabs = tabs = QTabWidget()
        for w, name in ((self.findings, "Findings"), (self.channel, "Channel"),
                        (self.modes, "Modes"), (self.transition, "Transition"),
                        (self.dynamics, "Dynamics"),
                        (self.range, "Range"), (self.tree, "Tree"),
                        (self.genomes, "Genomes"),
                        (self.variants, "Variants")):
            tabs.addTab(w, name)
        compact_combos(tabs)
        self._dock("Analysis", tabs, Qt.DockWidgetArea.RightDockWidgetArea, scroll=False)

        sp = self.structure_panel
        sp.load_requested.connect(self.load_structure)
        sp.style_changed.connect(self.restyle)
        sp.sites_toggled.connect(lambda *_: self._apply_sites())
        self.fills = FillController(self.scene, sp, self.statusBar().showMessage)
        self.channel.pore_toggled.connect(self.scene.show_pore)
        self.channel.states_requested.connect(self._compare_states)
        self.channel.unitary_requested.connect(self._unitary_states)
        self.channel.mutants_requested.connect(self._ryr_mutants)
        self.modes.compute_requested.connect(self.compute_modes)
        self.modes.animate_requested.connect(self._animate_mode)
        self.modes.stop_requested.connect(self.scene.stop_animation)
        tp = self.transition
        tp.build_requested.connect(self.build_transition)
        tp.preset_requested.connect(self._transition_preset)
        tp.frame_requested.connect(self.morph.show_frame)
        tp.play_requested.connect(lambda: self.morph.play(on_frame=tp.follow))
        tp.stop_requested.connect(self.morph.stop)
        tp.paint_toggled.connect(self._paint_displacement)
        self.findings.show_structure.connect(self._show_check)
        self.findings.showable = (frozenset(CHECK_SITES) | CHECK_TREE
                                  | frozenset(CHECK_GENOMES) | CHECK_RANGE)
        tabs.currentChanged.connect(self._tab_shown)
        self.variants.highlight_residue.connect(self._highlight_variant)
        self.variants.draw_requested.connect(self._draw_variants)
        self.viewport.atom_picked.connect(self._picked)
        self.viewport.scene_ready.connect(lambda _: self.scene.attach())
        self.viewport.status.connect(self.statusBar().showMessage)
        self.viewport.navigated.connect(self.scene.navigated)
        self.viewport.resized.connect(self.scene.resized)
        self._docks_sized = False
        self._menus()
        self.statusBar().showMessage(f"ip3r_genes results: {genes_results()}")

    def showEvent(self, event) -> None:
        super().showEvent(event)
        if not self._docks_sized:
            self._docks_sized = True
            self.resizeDocks(self._docks, list(DOCK_WIDTHS), Qt.Orientation.Horizontal)

    def _dock(self, title, widget, area, scroll=True) -> None:
        dock = QDockWidget(title, self)
        if scroll:
            area_w = QScrollArea()
            area_w.setWidget(widget)
            area_w.setWidgetResizable(True)
            # Fit the width, scroll only vertically: a sideways scroll hid subunit D.
            area_w.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            widget = area_w
        dock.setWidget(widget)
        dock.setMinimumWidth(300 if area == Qt.DockWidgetArea.LeftDockWidgetArea else 480)
        self.addDockWidget(area, dock)
        self._docks = [*getattr(self, "_docks", []), dock]

    def _menus(self) -> None:
        mb = self.menuBar()
        f = mb.addMenu("&File")
        self._action(f, "Open session…", self.sessions.open, "Ctrl+O")
        self._action(f, "Save session…", self.sessions.save, "Ctrl+Shift+S")
        f.addSeparator()
        self._action(f, "Fetch all registry structures", self._fetch_all)
        self._action(f, "Save screenshot…", self._screenshot, "Ctrl+S")
        f.addSeparator()
        self._action(f, "Quit", self.close, QKeySequence.StandardKey.Quit)
        v = mb.addMenu("&View")
        self._action(v, "Side view (cytosol up)", self.scene.side_view, "Ctrl+1")
        self._action(v, "Top view (down the pore)", self.scene.top_view, "Ctrl+2")
        self._action(v, "IP3 site", self._site_view, "Ctrl+3")
        self._action(v, "Fit to view", self.scene.fit_view, "Ctrl+0")
        self._action(v, "Toggle spin", lambda: self.viewport.set_spin(
            0.0 if self.viewport._spin_speed else 20.0), "Space")
        h = mb.addMenu("&Help")
        self._action(h, "Parameters…", self.edit_parameters, "Ctrl+Shift+P")
        self._action(h, "About", self._about)

    def edit_parameters(self) -> None:
        """The registry editor; the banner follows the registry by itself."""
        self.params_dialog = ParametersDialog(self)
        self.params_dialog.exec()

    def _action(self, menu, text, slot, shortcut=None) -> QAction:
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        a.triggered.connect(slot)
        menu.addAction(a)
        return a

    # ------------------------------------------------------------- loading

    def load_structure(self, pdb_id: str) -> None:
        self.statusBar().showMessage(f"loading {pdb_id}…")
        self.modes.clear()
        loader.ALLOW_FETCH = True

        def work():
            st = loader.load(pdb_id)
            return st, measure_channel(st)
        run_async(work, on_done=self._loaded,
                  on_error=lambda e: QMessageBox.warning(self, "Load failed", e))

    def _loaded(self, result) -> None:
        st, summary = result
        self.structure_panel.refresh_list()
        self.structure_panel.set_chains(st.chains)
        entry = get_entry(st.name)
        self.structure_panel.set_info(
            f"<b>{st.name}</b>: {entry.title if entry else ''}<br>"
            f"{st.n_atoms:,} atoms, {st.n_residues:,} residues, chains "
            f"{', '.join(st.chains)}")
        self.channel.show_summary(summary)
        self.channel.set_paralog(summary.numbering.paralog if summary.numbering else None)
        self.morph.reset()
        self.scene.set_structure(st, summary, self._style_kwargs())
        self.transition.set_start(st.name, summary.numbering.paralog
                                  if summary.numbering else None)
        self._apply_sites()
        self.fills.loaded(st)
        self.variants.follow(summary.numbering.paralog if summary.numbering else None)
        self.statusBar().showMessage(
            f"{st.name} loaded — numbering "
            f"{summary.numbering.paralog if summary.numbering else 'none'}; "
            "left-drag rotate, wheel zoom, shift-drag pan, click to identify")
        if self._pending_check:
            self._apply_check(self._pending_check)
            self._pending_check = None
        if self._pending_transition and self._pending_transition[0] == st.name:
            self.build_transition(*self._pending_transition[1:])
        self._pending_transition = None
        self.sessions.loaded(st)

    def _style_kwargs(self) -> dict:
        sp = self.structure_panel
        return {"style": sp.current_style(), "color_by": sp.current_color(),
                "layer": sp.current_layer(), "show_ligands": sp.ligands.isChecked(),
                "visible_chains": sp.visible_chains()}

    def restyle(self) -> None:
        self.scene.restyle(self._style_kwargs())

    def _apply_sites(self) -> None:
        on = [k for k, b in self.structure_panel.site_boxes.items() if b.isChecked()]
        self.scene.highlight_sites(on)

    # ---------------------------------------------------------- analysis

    def compute_modes(self) -> None:
        if self.scene.structure is None:
            self.statusBar().showMessage("load a structure first")
            return
        self.modes.set_busy("solving…")
        run_async(self.scene.compute_modes,
                  on_done=lambda out: self.modes.show_modes(*out),
                  on_error=lambda e: (self.modes.set_busy(""), QMessageBox.warning(
                      self, "Modes failed", e), self.modes.clear()))

    def _animate_mode(self, index: int, amplitude: float) -> None:
        if self.morph.result is not None:
            self.morph.show_frame(0)
            self.transition.follow(0)
        self.scene.animate_mode(index, amplitude)

    def build_transition(self, end_id: str, fit: str = "pore",
                         method: str = "restrained") -> None:
        st = self.scene.structure
        if st is None or not end_id:
            return
        loader.ALLOW_FETCH = True
        self.morph.reset()
        self.transition.set_busy(f"building {st.name} → {end_id} (morph and "
                                 "elastic network)…")
        run_async(build_transition, st, end_id, fit, method,
                  on_done=self._transition_built, on_error=self.transition.failed)

    def _transition_built(self, result) -> None:
        try:
            self.morph.install(result)
        except ValueError as exc:
            return self.transition.failed(str(exc))
        self.transition.show_result(result)
        if self.transition.paint.isChecked():
            self._paint_displacement(True)
        self.sessions.transition_built(result)

    def _transition_preset(self, start: str, end: str) -> None:
        tp = self.transition
        args = (end, tp.fit.currentData(), tp.method.currentData())
        if self.scene.structure is not None and self.scene.structure.name == start:
            self.build_transition(*args)
        else:
            self._pending_transition = (start, *args)
            self.structure_panel.select(start)

    def _paint_displacement(self, on: bool) -> None:
        from ..render.representations import ColorBy
        sp = self.structure_panel
        target = ColorBy.DISPLACEMENT if on else ColorBy.ELEMENT_DOMAIN
        sp.color.setCurrentIndex(sp.color.findData(target))

    def _compare_states(self) -> None:
        from ..structure.states import state_panel
        loader.ALLOW_FETCH = True
        self.channel.states_btn.setEnabled(False)
        paralog = self.channel.panel_paralog
        self.statusBar().showMessage(f"measuring the {paralog} state panel…")
        run_async(state_panel, paralog, on_done=self.channel.show_states,
                  on_error=lambda e: (self.channel.states_btn.setEnabled(True),
                                      QMessageBox.warning(self, "States failed", e)))

    def _unitary_states(self) -> None:
        from ..physics.unitary import unitary_panel
        loader.ALLOW_FETCH = True
        self.channel.unitary_btn.setEnabled(False)
        paralog = self.channel.panel_paralog
        self.statusBar().showMessage(f"solving drift-diffusion over each {paralog} state…")
        run_async(unitary_panel, paralog, sweep=True,
                  on_done=lambda rows: (self.statusBar().clearMessage(),
                                        self.channel.show_unitary(rows)),
                  on_error=lambda e: (self.channel.unitary_btn.setEnabled(True),
                                      QMessageBox.warning(self, "Conductance failed", e)))

    def _ryr_mutants(self) -> None:
        from ..physics.ryr_mutants import open_mutant_panel
        loader.ALLOW_FETCH = True
        self.channel.mutants_btn.setEnabled(False)
        self.statusBar().showMessage("modelling the RyR1 charge mutants…")
        run_async(open_mutant_panel,
                  on_done=lambda res: (self.statusBar().clearMessage(),
                                       self.channel.show_mutants(res)),
                  on_error=lambda e: (self.channel.mutants_btn.setEnabled(True),
                                      QMessageBox.warning(self, "Mutants failed", e)))

    def _tab_shown(self, i: int) -> None:
        w = self.tabs.widget(i)
        if w is self.tree or w is self.genomes or w is self.range:
            w.ensure_loaded()

    def _show_check(self, pdb_id: str, check_id: str) -> None:
        if check_id in CHECK_TREE:
            self.tabs.setCurrentWidget(self.tree)
            return
        if check_id in CHECK_RANGE:
            self.range.ensure_loaded()
            self.tabs.setCurrentWidget(self.range)
            return
        if check_id in CHECK_GENOMES:
            self.genomes.show_layer(CHECK_GENOMES[check_id])
            self.tabs.setCurrentWidget(self.genomes)
            return
        # "" = the check has no structure of its own: use the one on screen
        # (painted only if it is in human numbering, as every site is).
        shown = self.scene.structure.name if self.scene.structure is not None else None
        pdb_id = pdb_id or shown or DEFAULT_STRUCTURE
        if shown == pdb_id:
            self._apply_check(check_id)
        else:
            self._pending_check = check_id
            self.structure_panel.select(pdb_id)

    def _apply_check(self, check_id: str) -> None:
        site = CHECK_SITES.get(check_id)
        for k, b in self.structure_panel.site_boxes.items():
            b.setChecked(k == site)
        if check_id == "S0.pore_profile":
            self.channel.show_pore.setChecked(True)
        if check_id in CHECK_COLOURS:
            from ..render.representations import ColorBy
            sp = self.structure_panel
            sp.color.setCurrentIndex(sp.color.findData(ColorBy(CHECK_COLOURS[check_id])))
        if check_id in CHECK_SITE_VIEW:
            self._site_view()

    def _site_view(self) -> None:
        msg = self.scene.site_view()
        if msg:
            self.statusBar().showMessage(
                f"camera: {msg}" + ("; front clipped at the pocket (Ctrl+0 fits all)"
                                    if self.scene.fit_target == "site" else ""))

    def _highlight_variant(self, paralog: str, resi: int) -> None:
        msg = self.scene.highlight_residue(paralog, resi)
        self.variants.report(msg)

    def _draw_variants(self, paralog: str, classes: tuple, layer) -> None:
        msg = self.scene.show_variants(paralog, classes, layer)
        if msg:
            self.variants.report(msg)

    def _picked(self, index: int) -> None:
        self.statusBar().showMessage(self.scene.describe_atom(index))

    # --------------------------------------------------------------- misc

    def _fetch_all(self) -> None:
        self.statusBar().showMessage("fetching registry structures…")
        run_async(fetch_all, log=lambda *_: None, on_done=lambda s: (
            self.structure_panel.refresh_list(),
            self.statusBar().showMessage(f"fetched: {s}")))

    def _screenshot(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save screenshot", "ip3r.png",
                                              "PNG (*.png)")
        if path:
            self.viewport.grabFramebuffer().save(path)

    def _about(self) -> None:
        QMessageBox.about(self, "About", (
            f"<b>IP3R Structural Simulator {__version__}</b><p>Physics-driven "
            "3-D model of the IP3 receptor, and a re-derivation of the "
            f"ip3r_genes results found at<br><code>{genes_results()}</code>.</p>"
            "<p>Ported from the PIEZO1 simulator. See README.md and "
            "INTERFACE.md.</p>"))
