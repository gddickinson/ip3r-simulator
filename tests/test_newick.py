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


def _split_set(root):
    """Every bipartition as a frozenset of the side without the first tip."""
    tips = {lf.label for lf in leaves(root)}
    anchor = sorted(tips)[0]
    out = set()
    for n in _walk(root):
        side = {lf.label for lf in leaves(n)}
        if 1 < len(side) < len(tips) - 1:
            out.add(frozenset(side if anchor not in side else tips - side))
    return out


def _walk(n):
    yield n
    for c in n.children:
        yield from _walk(c)


def test_reroot_keeps_every_split_its_support_and_the_tree_length():
    from ip3r.analysis.newick import clade_node, reroot
    root = parse(TREE)
    new = reroot(root, {"D", "E"})
    assert {lf.label for lf in leaves(new.children[0])} == {"D", "E"}
    assert _split_set(new) == _split_set(root)
    # the (C,(D,E)) edge was labelled "" and the (A,B) edge 90/99: both kept
    assert clade_node(new, {"A", "B"}).label == "90/99"
    total = lambda r: sum(n.length or 0 for n in _walk(r))
    assert total(new) == pytest.approx(total(root))
    # an outgroup given as the complement of a clade works; a non-clade fails
    assert {lf.label for lf in leaves(reroot(root, {"C x", "D", "E", "F"}).children[0])} \
        in ({"A", "B"}, {"C x", "D", "E", "F"})
    with pytest.raises(ValueError):
        reroot(root, {"A", "D"})


def test_reroot_moves_a_label_to_the_flipped_side_of_its_edge():
    """Rooted at A, the (A,B) edge's 90/99 now sits on the node holding the
    other four tips: the same bipartition, seen from the new root."""
    from ip3r.analysis.newick import clade_node, reroot
    new = reroot(parse(TREE), {"A"})
    assert clade_node(new, {"C x", "D", "E", "F"}).label == "90/99"
    assert clade_node(new, {"B", "C x", "D", "E", "F"}) is not None
