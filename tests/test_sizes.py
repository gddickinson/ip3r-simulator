"""The project rule: no Python file over 500 lines."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_no_file_over_500_lines():
    long = [(p, n) for p in ROOT.glob("**/*.py")
            if ".git" not in p.parts
            and (n := len(p.read_text().splitlines())) > 500]
    assert not long, long
