#!/usr/bin/env python
"""Scripted GUI smoke test: launch, load a structure, drive panels, screenshot.

Mechanical refactors of Qt code break the GUI silently (the PIEZO1 project
learned this twice), so this runs the real application with a real OpenGL
context and exits non-zero if any step fails. It also produces the README
screenshots.

    python scripts/screenshot_app.py [--steps lumen,extras] [--checks] [--list]

The steps are in named groups (``screenshot_groups.py``); ``--steps`` runs
the groups a change touches, each set up by itself, and ``--checks`` adds
the findings checks (slow, and verified headlessly by ``make checks``).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import screenshot_groups as groups  # noqa: E402
from screenshot_awake import keep_awake  # noqa: E402

#: A hang detector, not a speed test. The longest waits are the findings
#: checks under --checks (~13 min, one step) and 8TKF's dielectric lumen
#: (Round 7.14, ~5 min).
STEP_LIMIT_S = 1200


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="step groups:\n" + groups.describe(),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--structure", default="6DQN")
    ap.add_argument("--out", default="docs/img")
    ap.add_argument("--checks", action="store_true",
                    help="run every findings check in the GUI first (~13 min); "
                         "the findings screenshot is saved only then")
    ap.add_argument("--steps", default="all",
                    help="comma-separated step groups (default all; see below)")
    ap.add_argument("--list", action="store_true", help="list the step groups")
    args = ap.parse_args()
    if args.list:
        print(groups.describe())
        return 0
    plan = groups.schedule(groups.select(args.steps))

    from types import SimpleNamespace

    from PyQt6.QtCore import QTimer
    from PyQt6.QtWidgets import QApplication

    from ip3r.config import SETTINGS
    from ip3r.ui.gl_widget import configure_surface_format
    from ip3r.ui.main_window import MainWindow
    from ip3r.ui.theme import apply_dark_theme

    import tempfile

    awake = keep_awake()
    tmp = tempfile.TemporaryDirectory()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    configure_surface_format(SETTINGS.render)
    app = QApplication(sys.argv[:1])
    apply_dark_theme(app)
    win = MainWindow()
    win.resize(1600, 980)
    win.show()
    win.raise_()
    ctx = SimpleNamespace(args=args, state={}, tmp=tmp)
    run = {"i": 0, "errors": [], "began": None, "t0": time.monotonic(),
           "groups": {}}
    print(f"{len(plan)} steps in {', '.join(dict.fromkeys(g for g, _, _ in plan))}"
          f"{'' if awake else ' (App Nap not held)'}", file=sys.stderr, flush=True)

    def fail(msg):
        run["errors"].append(msg)
        print("FAIL:", msg)

    def step():
        i = run["i"]
        if i >= len(plan):
            summary = ", ".join(f"{g} {t:.0f} s" for g, t in run["groups"].items())
            print(f"screenshots written to {out}; {summary}")
            return app.quit()
        group, label, fn = plan[i]
        now = time.monotonic()
        if run["began"] is None:
            run["began"] = now
            print(f"{group}: {label} at {now - run['t0']:.0f} s", file=sys.stderr, flush=True)
        try:
            again = fn(win, app, out, ctx)
        except Exception as exc:                 # report and stop
            fail(f"{group}, {label}: {type(exc).__name__}: {exc}")
            return app.quit()
        if again:
            if time.monotonic() - run["began"] > STEP_LIMIT_S:
                fail(f"{group}, {label}: still waiting after {STEP_LIMIT_S} s")
                return app.quit()
            return QTimer.singleShot(500, step)
        run["groups"][group] = run["groups"].get(group, 0.0) + time.monotonic() - run["began"]
        run["i"], run["began"] = i + 1, None
        QTimer.singleShot(1500, step)

    QTimer.singleShot(800, step)
    app.exec()
    return 1 if run["errors"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
