"""Write tests/fixtures/anm_reduced.json and tests/fixtures/gating_curves.json.

The reduced elastic network (``ip3r.web.reduced_anm``) and the three gating
models' steady-state curves, computed by the repo's Python on exactly the
inputs the browser uses. ``ip3r/web/static/*.test.mjs`` hold the JavaScript
engines to these numbers; ``tests/test_web_data.py`` checks they are still
what the Python gives.

Run: python scripts/make_anm_fixture.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ip3r.physics import gating, gating_mak, park_drive  # noqa: E402
from ip3r.web.app import PRECOMPUTED  # noqa: E402
from ip3r.web.gating_api import MODELS, model_constants, reference_bells  # noqa: E402
from ip3r.web.reduced_anm import FIXTURE_CASES, reduced_modes  # noqa: E402

FIXTURES = ROOT / "tests" / "fixtures"
#: Ca2+ grid (µM) and IP3 levels (µM) of the gating fixture.
CA_GRID = np.logspace(-3, 3, 61)
IP3_LEVELS = (0.033, 0.3, 1.0, 10.0)


def anm_fixture() -> dict:
    cases = []
    for pdb, stride, cutoff, n in FIXTURE_CASES:
        rec = json.loads((PRECOMPUTED / f"structure_{pdb}.json").read_text())
        cases.append(reduced_modes(rec, stride, cutoff, n))
    return {"cases": cases}


def gating_fixture() -> dict:
    curves = {"dyk": {}, "mak": {}, "pd": {}}
    for p in IP3_LEVELS:
        curves["dyk"][str(p)] = gating.open_probability(CA_GRID, p).tolist()
        curves["mak"][str(p)] = gating_mak.open_probability(CA_GRID, p).tolist()
        curves["pd"][str(p)] = park_drive.open_probability(CA_GRID, p).tolist()
    park = [park_drive.park_fraction(float(c), 1.0) for c in CA_GRID]
    return {"constants": {m: model_constants(m)["values"] for m in MODELS},
            "ca": CA_GRID.tolist(), "ip3": list(IP3_LEVELS), "curves": curves,
            "park_fraction_ip3_1": park, "bells": reference_bells(IP3_LEVELS)}


def main() -> None:
    FIXTURES.mkdir(exist_ok=True)
    for name, obj in (("anm_reduced.json", anm_fixture()), ("gating_curves.json", gating_fixture())):
        path = FIXTURES / name
        path.write_text(json.dumps(obj, separators=(",", ":")))
        print(f"{path.relative_to(ROOT)}: {path.stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
