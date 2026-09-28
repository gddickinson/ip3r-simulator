"""The main window's Round 7.18 additions, as a mixin (split from
``main_window.py`` to keep it under the size limit): the right-click menu,
the sequence window, the selection's clearing, the HUD switches, the guide,
and the Analyses menu's result windows."""

from __future__ import annotations

from .context_menu import build_context_menu

__all__ = ["WindowExtras"]


class WindowExtras:
    """Needs ``scene``, ``selection``, ``hud``, ``viewport``, ``statusBar``,
    ``sequence_window``, ``help_dialog``, ``result_windows``,
    ``measure_action`` on the window."""

    def context_menu(self, pos, index: int) -> None:
        build_context_menu(self, index).exec(self.viewport.mapToGlobal(pos))

    # ------------------------------------------------ Round 7.18 windows

    def show_sequence(self, chain: str | None = None, residue: int | None = None) -> None:
        from .sequence_window import SequenceWindow
        if self.sequence_window is None:
            self.sequence_window = SequenceWindow(self.selection, self)
            st, s = self.scene.structure, self.scene.summary
            if st is not None:
                self.sequence_window.set_structure(
                    st, s.numbering.paralog if s.numbering else None)
        self.sequence_window.show_residue(chain, residue)

    def escape(self) -> None:
        """Esc: stop a recording; else out of full screen if in it; else
        clear the selection."""
        if self.movies.busy:
            self.movies.cancel()
        elif self.presentation.active:
            self.presentation.leave()
        else:
            self.clear_selection()

    def clear_selection(self) -> None:
        self.selection.clear()
        self.selection.clear_distances()
        if self.selection.measuring:
            self.measure_action.setChecked(False)
            self.selection.arm(False)

    def set_hud(self, key: str, on: bool) -> None:
        setattr(self.hud.settings, key, bool(on))
        self.hud.update()

    def show_help(self, topic: str | None = None) -> None:
        from .help_dialog import HelpDialog
        if self.help_dialog is None:
            self.help_dialog = HelpDialog(self)
        if topic:
            self.help_dialog.show_topic(topic)
        self.help_dialog.show()
        self.help_dialog.raise_()

    def open_analysis(self, analysis):
        """An Analyses-menu entry: its window, run at once."""
        from .analyses import command
        st, s = self.scene.structure, self.scene.summary
        paralog = s.numbering.paralog if s is not None and s.numbering else None
        try:
            args = command(analysis, st.name if st is not None else None, paralog)
        except ValueError as exc:
            self.statusBar().showMessage(str(exc))
            return None
        return self.open_command(args, analysis.label, analysis.about, analysis.duration)

    def open_command(self, args, title: str, about: str = "", duration: str = ""):
        from .result_window import ResultWindow
        st = self.scene.structure
        w = ResultWindow(title, list(args), about, duration,
                         pdb=st.name if st is not None else None, parent=self)
        self.result_windows = [r for r in self.result_windows if r.isVisible()] + [w]
        w.show()
        w.run()
        return w
