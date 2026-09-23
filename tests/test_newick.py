import pytest

from ip3r.analysis.newick import leaves, mrca, parse, smallest_clade_containing

TREE = "((A:1,B:2)90/99:0.5,('C x':1,(D,E)x:1):0.3,F);"


def test_parse_keeps_labels_and_lengths():
    root = parse(TREE)
    assert sorted(lf.label for lf in leaves(root)) == ["A", "B", "C x", "D", "E", "F"]
    assert root.children[0].label == "90/99"
    assert root.children[0].length == 0.5


def test_mrca_and_clade():
    root = parse(TREE)
    assert smallest_clade_containing(root, {"D", "C x"}) == {"C x", "D", "E"}
    assert smallest_clade_containing(root, {"A", "D"}) == {"A", "B", "C x", "D", "E", "F"}


def test_missing_label_raises():
    with pytest.raises(KeyError):
        mrca(parse(TREE), {"A", "Z"})


def test_requires_semicolon():
    with pytest.raises(ValueError):
        parse("(A,B)")
