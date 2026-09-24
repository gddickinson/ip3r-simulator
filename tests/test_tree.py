"""The tree questions on small trees whose answers are known by construction."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from ip3r.analysis.newick import parse  # noqa: E402
from ip3r.analysis.tree import (bipartitions, cyclostome_clades, genus_of,  # noqa: E402
                                group_of, is_supported, paralog_clades, support)
from ip3r.analysis.tree_figure import draw_tree, layout  # noqa: E402

HAG, LAM = "vertebrate_basal_Myxine_glutinosa_", "vertebrate_basal_Petromyzon_marinus_"
# RyR outgroup | (invertebrate, (vertebrates)); in the vertebrates a
# cyclostome pair branches first, and a second pair is sister to ITPR2+ITPR3.
_RYR = "(RYR_a:1,RYR_b:1)100/100:1"
_CYC1 = f"({HAG}1:1,{LAM}1:1)99/100:1"
_CYC2 = f"({HAG}2:1,{LAM}2:1)100/100:1"
_P1 = "(ITPR1_h:1,ITPR1_m:1,vertebrate_basal_Callorhinchus_s:1)100/100:1"
_P23 = "((ITPR2_h:1,ITPR2_m:1)100/100:1,(ITPR3_h:1,ITPR3_m:1)97/96:1)70/90:1"
_JOIN = f"({_CYC2},{_P23})83.1/77:1"
_VERT = f"({_CYC1},({_P1},{_JOIN})90/99:1)60/70:1"
TOY = f"({_RYR},(invert_metazoa_x:1,{_VERT})100/100:1);"


def test_groups_and_genera():
    assert group_of("vertebrate_basal_Myxine_x") == "vertebrate_basal"
    assert group_of("invert_metazoa_Dugesia") == "invert_metazoa"
    assert group_of("ITPR2_Homo") == "ITPR2" and group_of("weird") == "other"
    assert genus_of("vertebrate_basal_Petromyzon_marinus_A") == "Petromyzon"
    assert genus_of("ITPR1_Homo_sapiens") is None


def test_support_parsing_and_bar():
    root = parse("((a:1,b:1)83.1/77:1,c:1);")
    n = root.children[0]
    assert support(n) == (83.1, 77.0) and not is_supported(n)
    assert support(parse("((a,b)x,c);").children[0]) is None


def test_root_edge_is_one_bipartition():
    # both root children are labelled: they are the same split, counted once
    root = parse(TOY)
    nodes = bipartitions(root)
    labelled = [n for n in _internal(root) if support(n)]
    assert len(nodes) == len(labelled) - 1
    assert root.children[1] not in nodes


def test_paralog_clades_take_in_unnamed_vertebrates():
    cl = paralog_clades(parse(TOY))
    assert len(cl["ITPR1"].tips) == 3 and cl["ITPR1"].unnamed and not cl["ITPR1"].foreign
    assert cl["ITPR2"].support == (100.0, 100.0) and len(cl["ITPR3"].tips) == 2


def test_a_misplaced_tip_is_foreign():
    bad = TOY.replace("ITPR1_m", "ITPR2_x").replace("ITPR2_m", "ITPR1_m")
    cl = paralog_clades(parse(bad))
    assert cl["ITPR1"].foreign and cl["ITPR2"].foreign


def test_cyclostome_clades():
    cc = cyclostome_clades(parse(TOY))
    assert [len(c.tips) for c in cc] == [2, 2]
    assert all(c.genera == ["Myxine", "Petromyzon"] for c in cc)
    first = [c for c in cc if c.first_among_vertebrates]
    assert len(first) == 1 and first[0].support == (99.0, 100.0)
    join = next(c for c in cc if not c.first_among_vertebrates)
    assert set(join.sister) == {"ITPR2", "ITPR3"} and join.parent_support == (83.1, 77.0)


def test_split_cyclostomes_are_separate_clades():
    # swap a lamprey into ITPR3: two lone tips, neither holding both genera
    split = TOY.replace(f"{LAM}2", "TMP").replace("ITPR3_m", f"{LAM}2").replace("TMP", "ITPR3_m")
    cc = cyclostome_clades(parse(split))
    assert sorted(len(c.tips) for c in cc) == [1, 1, 2]


def test_layout_and_drawing():
    root = parse(TOY)
    lay = layout(root)
    assert sorted(lay.y[id(t)] for t in lay.tips) == list(range(len(lay.tips)))
    assert lay.x[id(root)] == 0.0 and max(lay.x.values()) == pytest.approx(7.0)
    fig, ax = plt.subplots()
    info = draw_tree(ax, root, labels=True)
    plt.close(fig)
    assert info["clades"] == {"ITPR1": 3, "ITPR2": 2, "ITPR3": 2}
    assert info["cyclostome_clades"] == [2, 2]
    assert info["supported_nodes"] == sum(map(is_supported, bipartitions(root)))


def _internal(root):
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        if n.children and n is not root:
            out.append(n)
        stack.extend(n.children)
    return out
