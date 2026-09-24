from ip3r.cli import main


def test_gating(capsys):
    assert main(["gating", "--ip3", "1.0"]) == 0
    assert "peak P_open" in capsys.readouterr().out


def test_params(capsys):
    assert main(["params"]) == 0
    assert "gating.d1" in capsys.readouterr().out


from conftest import needs_structure  # noqa: E402


@needs_structure("8TKG", "8TKF")
def test_transition(capsys):
    assert main(["transition", "8TKG", "8TKF"]) == 0
    out = capsys.readouterr().out
    assert "lowest collective A mode" in out and "random direction" in out


def test_puffs_park_drive(capsys):
    assert main(["puffs", "--model", "park-drive", "--duration", "1"]) == 0
    assert '"large_share"' in capsys.readouterr().out
