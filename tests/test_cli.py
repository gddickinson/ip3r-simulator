from ip3r.cli import main


def test_gating(capsys):
    assert main(["gating", "--ip3", "1.0"]) == 0
    assert "peak P_open" in capsys.readouterr().out


def test_params(capsys):
    assert main(["params"]) == 0
    assert "gating.d1" in capsys.readouterr().out
