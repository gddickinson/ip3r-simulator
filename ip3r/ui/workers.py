"""Run slow work off the Qt main thread.

Nothing that parses a structure, solves for normal modes, integrates a
model or runs a check may block the event loop. :func:`run_async` runs a
callable on a ``QThreadPool`` worker and delivers its result (or the
exception) back on the main thread through a signal.
"""

from __future__ import annotations

import traceback

from PyQt6.QtCore import QObject, QRunnable, QThreadPool, pyqtSignal

__all__ = ["run_async"]


class _Signals(QObject):
    done = pyqtSignal(object)
    failed = pyqtSignal(str)


class _Job(QRunnable):
    def __init__(self, fn, args, kwargs):
        super().__init__()
        self.fn, self.args, self.kwargs = fn, args, kwargs
        self.signals = _Signals()

    def run(self):
        try:
            out = self.fn(*self.args, **self.kwargs)
        except Exception as exc:  # delivered to the UI, never swallowed
            self.signals.failed.emit(f"{type(exc).__name__}: {exc}\n"
                                     + traceback.format_exc(limit=3))
            return
        self.signals.done.emit(out)


_JOBS: set = set()


def run_async(fn, *args, on_done=None, on_error=None, **kwargs) -> None:
    job = _Job(fn, args, kwargs)
    _JOBS.add(job)                        # keep the signals object alive

    def _finish(_=None):
        _JOBS.discard(job)
    if on_done:
        job.signals.done.connect(on_done)
    if on_error:
        job.signals.failed.connect(on_error)
    job.signals.done.connect(_finish)
    job.signals.failed.connect(_finish)
    QThreadPool.globalInstance().start(job)
