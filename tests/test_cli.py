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


def test_microdomain(capsys):
    assert main(["microdomain", "--duration", "2", "--clamp", "none"]) == 0
    out = capsys.readouterr().out
    assert "store free" in out and "blip dF/F0" in out


@needs_structure("8TKG")
def test_lumen_of_a_shut_pore(capsys):
    assert main(["lumen", "8TKG"]) == 0
    assert "no K+ path joins the two baths" in capsys.readouterr().out
