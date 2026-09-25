"""Fills through an alignment (Round 7.9), and fills of long stretches.

Synthetic first. The prediction is the deposit's own chain as another isoform:
five residues spliced out and the rest renumbered. So the map is known, and a
fill must land exactly and carry the deposit's numbers. A stretch that needs
the missing residues must be refused for that reason, and so must one where
the model has extra residues.

Then real data. Rat 7LHF is filled from rat isoform 8, and every filled
residue is the construct's own. Hidden long stretches are right where
AlphaFold is confident. 8TKH's 926-943, the only resolved run below pLDDT 50,
is the case that says a very-low fill is not a position.
"""

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.core.structure import AA3TO1
from ip3r.structure.graft import Stretch, fill_stretches
from ip3r.structure.graft_numbering import by_alignment, by_number, runs
from ip3r.structure.numbering import chain_residues
from ip3r.structure.symmetry import rotation_matrix
from test_graft import _chain, _structure

AA1TO3 = {v: k for k, v in AA3TO1.items()}
SPLICE = (40, 44)


def _isoform(hidden=(25, 34), extra=False):
    """(deposit with a gap, its construct, an isoform lacking SPLICE, truth)."""
    ca, names = _chain(80, seed=3)
    keep = np.ones(len(ca), bool)
    keep[hidden[0] - 1:hidden[1]] = False
    idx = np.flatnonzero(keep)
    blocks = [("A", ca[p], names[p], int(p[0]) + 1, np.zeros(len(p), bool))
              for p in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)]
    dep = _structure(blocks, "DEP")
    construct = {i + 1: str(n) for i, n in enumerate(names)}
    rot = rotation_matrix(np.array([0.3, -1.0, 2.0]), 0.7)
    iso = np.ones(len(ca), bool)
    iso[SPLICE[0] - 1:SPLICE[1]] = False
    pca, pnames = ca[iso], names[iso]
    if extra:                                   # two residues inserted at 30
        pca = np.insert(pca, 29, pca[29] + [0.5, 0.5, 0.0], axis=0)
        pca = np.insert(pca, 29, pca[29] + [0.4, 0.2, 0.0], axis=0)
        pnames = np.insert(pnames, 29, ["GLY", "GLY"])
    pred = _structure([("A", pca @ rot.T + [3, 1, -4], pnames, 1,
                        np.zeros(len(pca), bool))], "AF-ISO", plddt=80.0)
    return dep, construct, pred, ca


def test_the_alignment_map_skips_the_splice_segment():
    _, construct, pred, _ = _isoform()
    nm = by_alignment(construct, chain_residues(pred, "A"))
    assert nm.route == "alignment" and nm.identity == 1.0
    assert nm.unmapped == (SPLICE,)
    assert nm(39) == 39 and nm(45) == 40 and nm(80) == 75 and nm(42) is None
    # By number the same isoform agrees only up to the splice.
    dep, _, _, _ = _isoform()
    assert by_number(chain_residues(dep, "A"), chain_residues(pred, "A")).identity < 0.9


def test_a_fill_through_the_map_lands_in_deposit_numbers():
    dep, construct, pred, truth = _isoform(hidden=(55, 64))   # after the splice
    nm = by_alignment(construct, chain_residues(pred, "A"))
    model = fill_stretches(dep, pred, [Stretch("A", 55, 64, "gap")], numbering=nm)
    assert len(model.fills) == 1 and not model.skipped
    f = model.fills[0]
    assert f.anchor_rmsd < 1e-3 and f.joined
    a = model.atoms
    ca = a.atom_name == "CA"
    assert a.res_seq[ca].tolist() == list(range(55, 65))
    assert np.allclose(a.xyz[ca], truth[54:64], atol=1e-3)
    assert all(construct[int(r)] == n for r, n in zip(a.res_seq[ca], a.res_name[ca]))


def test_by_number_the_same_isoform_would_fill_the_wrong_residues():
    dep, construct, pred, truth = _isoform(hidden=(55, 64))
    model = fill_stretches(dep, pred, [Stretch("A", 55, 64, "gap")])
    # Anchors must match in sequence, so a shifted numbering is refused, not misfilled.
    assert not model.fills and "anchors" in model.skipped[0].reason


def test_a_stretch_needing_the_spliced_residues_is_refused_by_name():
    dep, construct, pred, _ = _isoform(hidden=(38, 47))
    nm = by_alignment(construct, chain_residues(pred, "A"))
    model = fill_stretches(dep, pred, [Stretch("A", 38, 47, "gap")], numbering=nm)
    assert not model.fills
    assert model.skipped[0].reason.startswith("the model lacks residues 40-44")


def test_an_insertion_in_the_model_is_refused():
    dep, construct, pred, _ = _isoform(hidden=(25, 34), extra=True)
    nm = by_alignment(construct, chain_residues(pred, "A"))
    model = fill_stretches(dep, pred, [Stretch("A", 25, 34, "gap")], numbering=nm)
    assert not model.fills and "insertion" in model.skipped[0].reason


def test_runs():
    assert runs([5, 1, 2, 3, 7, 8]) == [(1, 3), (5, 5), (7, 8)]


# ---------------------------------------------------------------- real data

def _models():
    from ip3r.io.predictions import local_predictions
    if not local_predictions():
        pytest.skip("no AlphaFold model downloaded")


@needs_structure("7LHF")
def test_the_construct_includes_what_was_not_built():
    from ip3r.io.loader import local_path
    from ip3r.io.poly_seq import construct_sequence
    c = construct_sequence(str(local_path("7LHF")), "A")
    assert len(c) == 2736 and min(c) == 6 and max(c) == 2741
    assert c[324] == "GLU"                        # unmodelled, still named


@needs_structure("7LHF")
def test_7lhf_is_filled_from_rat_isoform_8_through_the_alignment():
    from ip3r.io.loader import load
    from ip3r.structure.graft import fill_structure, prediction_for
    from ip3r.structure.graft_numbering import construct_of
    _models()
    st = load("7LHF")
    pred, nm = prediction_for(st)
    assert pred.name == "AF-P29994-8-F1" and nm.route == "alignment"
    assert nm.identity == 1.0 and nm.unmapped == ((318, 332), (1693, 1732))
    model = fill_structure(st, "gaps", prediction=pred, numbering=nm)
    assert len(model.fills) == 40 and model.n_residues == 1264
    construct = construct_of(st, "A")
    a = model.atoms
    ca = (a.atom_name == "CA") & (a.chain == "A")
    assert all(construct[int(r)] == n for r, n in zip(a.res_seq[ca], a.res_name[ca]))
    lacks = {sk.stretch.label(): sk.reason for sk in model.skipped if sk.stretch.chain == "A"}
    assert "318-332" in lacks["A:324-351"] and "1693-1732" in lacks["A:1691-1723"]


@needs_structure("8TKG")
def test_long_windows_fill_where_alphafold_is_confident():
    from ip3r.io.loader import load
    from ip3r.parameters import PARAMETERS as P
    from ip3r.structure.graft import prediction_for
    from ip3r.structure.graft_calibration import window_trials
    _models()
    st = load("8TKG")
    pred, nm = prediction_for(st)
    trials = window_trials(st, pred, lengths=(30, 60), numbering=nm)
    assert len(trials) >= 30
    for n in (30, 60):
        ts = [t for t in trials if t.stretch.n_residues == n]
        assert np.median([t.rmsd_fill for t in ts]) < 2.0 < 10.0 < \
            np.median([t.rmsd_line for t in ts])
        assert np.median([t.plddt_scored for t in ts]) > 70   # ordered, so easy
    # True seams of right fills close within the join tolerance at any length.
    good = [max(t.joins) for t in trials if t.rmsd_fill < 2.0]
    assert max(good) <= P.value("graft.join_tolerance")


@needs_structure("8TKH", "7T3P")
def test_islands_are_right_when_confident_and_not_when_very_low():
    from ip3r.io.loader import load
    from ip3r.structure.graft import prediction_for
    from ip3r.structure.graft_calibration import island_trials
    _models()
    st = load("8TKH")
    pred, nm = prediction_for(st)
    low = dict(island_trials(st, pred, numbering=nm))[(926, 943)]
    assert low.stretch.n_residues > 50 and low.plddt_scored < 50
    assert not low.beats_line and low.rmsd_fill > 20.0
    st = load("7T3P")
    pred, nm = prediction_for(st)
    sure = [t for _, t in island_trials(st, pred, numbering=nm) if t.plddt_scored >= 70]
    assert len(sure) >= 10 and all(t.beats_line for t in sure)
    assert np.median([t.rmsd_fill for t in sure]) < 2.0
