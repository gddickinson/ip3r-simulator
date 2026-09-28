"""The two checks ip3r_genes' S29 answered: each still reads the
pre-S29 publication as a discrepancy, so the confirmation now is the
correction's, not a softened check."""

from pathlib import Path

from ip3r.analysis import checks as C
from conftest import needs_genes, needs_structure
from test_checks_calibration import _clear_caches, _copy_sources


def _setup(check_id, tmp_path, monkeypatch):
    check = {c.id: c for c in C.all_checks()}[check_id]
    (tmp_path / "results").mkdir()
    _copy_sources(check, tmp_path)
    monkeypatch.setenv("IP3R_GENES_DIR", str(tmp_path))
    _clear_caches()
    return check


@needs_genes
def test_old_both_metrics_sentence_is_a_discrepancy(tmp_path, monkeypatch):
    check = _setup("P5.report_both_metrics", tmp_path, monkeypatch)
    assert C.run_check(check).status == "confirmed"
    report = Path(tmp_path, "results/constraint/report.md")
    text = report.read_text()
    report.write_text(text + "\nThe gate and the selectivity filter are the most "
                      "constrained elements of the protein, on both metrics and "
                      "in all three paralogs.\n")
    _clear_caches()
    out = C.run_check(check)
    assert out.status == "discrepancy" and "pre-S29" in out.outcome.detail
    report.write_text("neither sentence")
    _clear_caches()
    assert C.run_check(check).status == "not_run"
    _clear_caches()


@needs_genes
@needs_structure("6DQN", "8TKG", "8TKF", "8TKH", "7T3P", "8TLA")
def test_without_s29s_table_the_old_claim_is_a_discrepancy(tmp_path, monkeypatch):
    check = _setup("P6.contacts_heavy_atom", tmp_path, monkeypatch)
    assert C.run_check(check).status == "confirmed"
    Path(tmp_path, "results/ligand_site/contact_rule.tsv").unlink()
    _clear_caches()
    out = C.run_check(check)
    assert out.status == "discrepancy" and "503 at 4.78" in out.outcome.found
    _clear_caches()
