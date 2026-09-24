"""The RyR1 spark steps of the GUI smoke test (``screenshot_app.py``).

Each step checks what the previous one started, saves its screenshot and
starts the next simulation. ``spark_step`` returns True while a worker is
still running (the caller retries the same step) and raises on a wrong
result. The Puffs view is already showing RyR1 sparks when step 0 runs.
"""

from __future__ import annotations

__all__ = ["SPARK_STEPS", "spark_step"]


def _waiting(pz, key: str = "coupled") -> bool:
    """True while the Puffs worker has not delivered ``key``."""
    if pz.result is None or key not in pz.result:
        if pz.text.text().startswith(("Simulating", "Triggering")):
            return True
        raise RuntimeError(pz.text.text())
    return False


def _rerun(pz, model: str, mg: float = 0.0) -> None:
    pz.model.setCurrentIndex(pz.model.findData(model))
    pz.mg.setValue(mg)
    pz.result = None
    pz.run()


def _mean_field(win, pz, out):
    c, u = pz.result["coupled"], pz.result["uncoupled"]
    if pz.result["model"] != "ryr1" or c["large"] == 0 or u["multi"]:
        raise RuntimeError(f"RyR1 sparks drew {pz.result}")
    if pz.mg.value() != 0 or not pz.mg.isEnabled() or pz.mg_scan_btn.isEnabled():
        raise RuntimeError("Mg2+ controls wrong on the mean-field receptor")
    win.grab().save(str(out / "gui_sparks.png"))
    _rerun(pz, "ryr1-cleft")


def _cleft(win, pz, out):
    c, u = pz.result["coupled"], pz.result["uncoupled"]
    if pz.result["model"] != "ryr1-cleft" or c["large"] == 0 or u["multi"]:
        raise RuntimeError(f"cleft sparks drew {pz.result}")
    win.grab().save(str(out / "gui_sparks_cleft.png"))
    _rerun(pz, "ryr1-cleft-fit")


def _fitted(win, pz, out):                # fitted to the measured bell: no end
    if pz.result["coupled"]["open_fraction"] < 0.3:
        raise RuntimeError(f"fitted cleft drew {pz.result}")
    win.grab().save(str(out / "gui_sparks_fitted.png"))
    _rerun(pz, "ryr1-cleft-fit", mg=1000.0)


def _fitted_mg(win, pz, out):             # the fibre's Mg2+: the array is quiet
    if pz.result["coupled"]["open_fraction"] > 0.05 or "Mg²⁺" not in pz.text.text():
        raise RuntimeError(f"fitted cleft under Mg2+ drew {pz.result}")
    pz.result = None
    pz.mg_scan_btn.click()


def _mg_scan(win, pz, out):               # triggered: Mg2+ shuts every spark
    rows = pz.result["mg_scan"]
    if rows[0].ended or rows[-1].ended != rows[-1].trials:
        raise RuntimeError(f"Mg2+ scan drew {[r.row() for r in rows]}")
    win.grab().save(str(out / "gui_sparks_mg.png"))


_STEPS = ((_mean_field, "coupled"), (_cleft, "coupled"), (_fitted, "coupled"),
          (_fitted_mg, "coupled"), (_mg_scan, "mg_scan"))
SPARK_STEPS = len(_STEPS)


def spark_step(k: int, win, app, out) -> bool:
    """Run spark step ``k``; True means "not ready yet, call again"."""
    check, key = _STEPS[k]
    pz = win.dynamics.puffs
    if _waiting(pz, key):
        return True
    app.processEvents()
    check(win, pz, out)
    return False
