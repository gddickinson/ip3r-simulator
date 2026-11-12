"""Rebuild the web site's precomputed results in ip3r/web/precomputed/.

Every file is produced by the repo's own code (``ip3r.web.export_*`` call
``measure_channel``, the ANM, the transition and the state panel); the web
app only serves them. Needs the deposits in ``ref/`` (``python -m ip3r fetch``
or ``--fetch`` here). Then refresh the reduced-model fixture:
``python scripts/make_anm_fixture.py``.

Run: python scripts/make_web_data.py [--fetch] [--only structures|modes|transition|states]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ip3r.config import RESOURCE_DIR  # noqa: E402
from ip3r.io import loader  # noqa: E402
from ip3r.web.app import PRECOMPUTED  # noqa: E402
from ip3r.web.export_dynamics import (MODE_DEPOSITS, modes_record,  # noqa: E402
                                      states_record, transition_record)
from ip3r.web.export_structures import (VIEWER_DEPOSITS, palette_record,  # noqa: E402
                                        structure_record)

PARTS = ("structures", "modes", "transition", "states")


def write(name: str, obj) -> None:
    path = PRECOMPUTED / name
    path.write_text(json.dumps(obj, ensure_ascii=False, separators=(",", ":")))
    print(f"  {name}: {path.stat().st_size / 1e3:.0f} kB")


def provenance() -> dict:
    params = (RESOURCE_DIR / "parameters.json").read_bytes()
    return {"generated_by": "scripts/make_web_data.py",
            "date": time.strftime("%Y-%m-%d"),
            "parameters_sha256": hashlib.sha256(params).hexdigest()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true", help="download missing deposits from RCSB")
    ap.add_argument("--only", choices=PARTS, action="append")
    args = ap.parse_args()
    loader.ALLOW_FETCH = args.fetch
    parts = args.only or PARTS
    PRECOMPUTED.mkdir(parents=True, exist_ok=True)
    if "structures" in parts:
        write("palette.json", palette_record())
        for pdb in VIEWER_DEPOSITS:
            write(f"structure_{pdb}.json", structure_record(pdb))
    if "modes" in parts:
        for pdb in MODE_DEPOSITS:
            write(f"modes_{pdb}.json", modes_record(pdb))
    if "transition" in parts:
        write("transition_8TKG_8TKF.json", transition_record("8TKG", "8TKF"))
    if "states" in parts:
        write("states_ITPR3.json", states_record())
    write("provenance.json", provenance())


if __name__ == "__main__":
    main()
