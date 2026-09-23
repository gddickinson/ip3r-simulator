import subprocess
import sys
from pathlib import Path

from ip3r.parameters import PARAMETERS

ROOT = Path(__file__).resolve().parent.parent


def test_registry_builds_and_validates():
    r = subprocess.run([sys.executable, str(ROOT / "scripts/build_parameters.py"), "--check"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout


def test_json_matches_the_table():
    sys.path.insert(0, str(ROOT / "scripts"))
    from parameter_table import P
    assert {p["key"]: p["value"] for p in P} == {k: p.default for k, p in PARAMETERS.parameters.items()}


def test_overrides_are_tracked():
    try:
        PARAMETERS.set_value("gating.d1", 0.2)
        assert PARAMETERS.modified
    finally:
        PARAMETERS.reset()
    assert not PARAMETERS.modified
