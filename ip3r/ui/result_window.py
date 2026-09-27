"""A window that runs one CLI analysis and shows its output (Round 7.18).

The command runs as ``python -m ip3r ...`` in a :class:`QProcess` with the
GUI's own interpreter, so it never touches the Qt main thread and can be
stopped. The child is handed the GUI's parameter set explicitly: the
overrides are written to a file and passed as ``IP3R_PARAMETERS``, or the
variable is removed when there are none, so a child cannot inherit a
parameter file the GUI has since moved away from. The window's first lines
are the stamp (:func:`~ip3r.ui.analyses.stamp`): the command, the deposit,
the parameter set.

The command line is editable: an entry is a starting point, and the
flags each command takes are in ``python -m ip3r <command> --help``.
"""

from __future__ import annotations

import os
import shlex
import sys
import tempfile
import time

from PyQt6.QtCore import QProcess, QProcessEnvironment, Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (QFileDialog, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
                             QPushButton, QVBoxLayout, QWidget)

from ..config import DATA_DIR, PROJECT_ROOT
from ..parameters import PARAMETERS
from .analyses import stamp

__all__ = ["ResultWindow"]


class ResultWindow(QWidget):
    def __init__(self, title: str, args: list[str], about: str = "",
                 duration: str = "", pdb: str | None = None, parent=None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setWindowTitle(f"Analysis — {title}")
        self.resize(900, 620)
        self.pdb = pdb
        self.proc: QProcess | None = None
        self._started = 0.0
        self._param_file: str | None = None
        lay = QVBoxLayout(self)
        note = QLabel(f"<b>{title}</b> — {about}"
                      + (f" <i>Takes {duration}.</i>" if duration else ""))
        note.setWordWrap(True)
        lay.addWidget(note)
        row = QHBoxLayout()
        row.addWidget(QLabel("python -m ip3r"))
        self.line = QLineEdit(" ".join(shlex.quote(a) for a in args))
        self.line.returnPressed.connect(self.run)
        row.addWidget(self.line, 1)
        self.run_btn = QPushButton("Run")
        self.run_btn.clicked.connect(self.run)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.clicked.connect(self.stop)
        self.save_btn = QPushButton("Save…")
        self.save_btn.clicked.connect(self._save)
        for b in (self.run_btn, self.stop_btn, self.save_btn):
            row.addWidget(b)
        lay.addLayout(row)
        self.out = QPlainTextEdit()
        self.out.setReadOnly(True)
        font = QFont("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        self.out.setFont(font)
        lay.addWidget(self.out, 1)
        self.state = QLabel("")
        lay.addWidget(self.state)
        self._buttons(False)

    # ---------------------------------------------------------------- run

    @property
    def args(self) -> list[str]:
        return shlex.split(self.line.text())

    @property
    def running(self) -> bool:
        return self.proc is not None and self.proc.state() != QProcess.ProcessState.NotRunning

    def run(self) -> None:
        if self.running:
            return
        args = self.args
        overrides = PARAMETERS.overrides()
        self.out.setPlainText(stamp(args, self.pdb, overrides))
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "1")
        env.remove("IP3R_PARAMETERS")
        if overrides:
            DATA_DIR.mkdir(parents=True, exist_ok=True)
            fd, self._param_file = tempfile.mkstemp(prefix="gui_params_", suffix=".json",
                                                    dir=DATA_DIR)
            os.close(fd)
            PARAMETERS.write_overrides(self._param_file)
            env.insert("IP3R_PARAMETERS", self._param_file)
        self.proc = QProcess(self)
        self.proc.setProcessEnvironment(env)
        self.proc.setWorkingDirectory(str(PROJECT_ROOT))
        self.proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self._read)
        self.proc.finished.connect(self._finished)
        self.proc.errorOccurred.connect(self._error)
        self._started = time.monotonic()
        self.proc.start(sys.executable, ["-m", "ip3r", *args])
        self._buttons(True)
        self.state.setText("running…")

    def stop(self) -> None:
        if self.running:
            self.proc.kill()
            self.state.setText("stopped")

    def wait(self, ms: int) -> bool:
        """Block until the process ends (the smoke test's use)."""
        return self.proc is not None and self.proc.waitForFinished(ms)

    def text(self) -> str:
        return self.out.toPlainText()

    # ------------------------------------------------------------ internals

    def _read(self) -> None:
        data = bytes(self.proc.readAllStandardOutput()).decode(errors="replace")
        self.out.moveCursor(self.out.textCursor().MoveOperation.End)
        self.out.insertPlainText(data)
        self.out.ensureCursorVisible()

    def _finished(self, code: int, status) -> None:
        self._read()
        secs = time.monotonic() - self._started
        how = "stopped" if status == QProcess.ExitStatus.CrashExit else f"exit {code}"
        self.state.setText(f"finished ({how}) in {secs:.1f} s")
        self._buttons(False)
        self._cleanup()

    def _error(self, err) -> None:
        if err == QProcess.ProcessError.FailedToStart:
            self.state.setText(f"could not start {sys.executable}")
            self._buttons(False)
            self._cleanup()

    def _cleanup(self) -> None:
        if self._param_file:
            try:
                os.unlink(self._param_file)
            except OSError:
                pass
            self._param_file = None

    def _buttons(self, running: bool) -> None:
        self.run_btn.setEnabled(not running)
        self.stop_btn.setEnabled(running)

    def _save(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save output", "analysis.txt",
                                              "Text (*.txt)")
        if path:
            with open(path, "w") as fh:
                fh.write(self.text())

    def closeEvent(self, event) -> None:              # noqa: N802
        self.stop()
        if self.proc is not None:
            self.proc.waitForFinished(2000)
        self._cleanup()
        super().closeEvent(event)
