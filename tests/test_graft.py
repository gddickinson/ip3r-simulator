"""AlphaFold fills: a synthetic chain with a known answer first, then real data.

The synthetic prediction is the deposit's own chain moved rigidly, so a fill
must land on the hidden residues exactly, close its seams at the chain's own
C-alpha spacing, and follow the deposit when it moves. Each rule that keeps a
fill honest (anchors must match in sequence, stubs are not anchors, clashes
are counted) is shown to fire on a case built to trip it.
"""

import numpy as np
import pytest

from conftest import needs_structure
from ip3r.core.structure import Structure
from ip3r.structure.graft import (GraftRefusal, Stretch, fill_stretches,
                                  fill_structure, placed_ca, prediction_for)
from ip3r.structure.graft_calibration import trial
from ip3r.structure.symmetry import rotation_matrix

NAMES = np.array(["LEU", "LYS", "GLU", "SER", "TRP", "ASP", "PHE", "ARG", "VAL", "ILE"])


def _chain(n=60, seed=0):
    """A helix-like C-alpha trace plus one side-chain atom per residue."""
    rng = np.random.default_rng(seed)
    t = np.arange(n) * 100 * np.pi / 180
    ca = np.column_stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)])
    names = NAMES[rng.integers(0, len(NAMES), n)]
    return ca, names


def _structure(blocks, name="SYN", plddt=0.0):
    """``blocks``: list of (chain, ca, names, first_resnum, stub_mask)."""
    rows = []
    for chain, ca, names, first, stub in blocks:
        for k, (x, nm) in enumerate(zip(ca, names)):
            r = first + k
            rows.append((chain, r, nm, "CA", "C", x))
            if not stub[k]:
                rows.append((chain, r, nm, "CG", "C", x + [1.5, 0.0, 0.0] * (x / (np.linalg.norm(x[:2]) + 1e-9))))
    n = len(rows)
    st = Structure(
        xyz=np.array([r[5] for r in rows], np.float32),
        element=np.array([r[4] for r in rows], "U2"),
        atom_name=np.array([r[3] for r in rows], "U6"),
        res_name=np.array([r[2] for r in rows], "U5"),
        res_seq=np.array([r[1] for r in rows], np.int32),
        chain=np.array([r[0] for r in rows], "U6"),
        hetero=np.zeros(n, bool), b_factor=np.full(n, plddt, np.float32),
        occupancy=np.ones(n, np.float32), alt_loc=np.full(n, ".", "U2"),
        entity=np.full(n, "1", "U6"), name=name)
    st._build_residue_index()
    return st


def _pair(hidden=(25, 34), stub=()):
    """(deposit with a gap, prediction = whole chain moved rigidly, truth)."""
    ca, names = _chain()
    keep = np.ones(len(ca), bool)
    keep[hidden[0] - 1:hidden[1]] = False
    stub_mask = np.isin(np.arange(1, len(ca) + 1), stub)
    offset = np.array([40.0, 0.0, 0.0])
    blocks = []
    for ch, shift in (("A", 0 * offset), ("B", offset)):
        c = ca + shift
        idx = np.flatnonzero(keep)
        for part in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
            blocks.append((ch, c[part], names[part], int(part[0]) + 1, stub_mask[part]))
    dep = _structure(blocks, "DEP")
    rot = rotation_matrix(np.array([1.0, 2.0, 3.0]), 1.1)
    pred = _structure([("A", ca @ rot.T + [5, -7, 11], names, 1,
                        np.zeros(len(ca), bool))], "AF-SYN", plddt=80.0)
    return dep, pred, ca


def test_a_fill_lands_on_the_hidden_residues():
    dep, pred, truth = _pair()
    model = fill_structure(dep, "gaps", prediction=pred)
    assert len(model.fills) == 2 and not model.skipped and model.identity == 1.0
    fa = next(f for f in model.fills if f.stretch.chain == "A")
    assert fa.anchor_rmsd < 1e-4 and fa.plddt == 80.0 and fa.confident == 1.0
    got = placed_ca(model, fa)
    assert sorted(got) == list(range(25, 35))
    assert max(np.linalg.norm(got[r] - truth[r - 1]) for r in got) < 1e-3
    step = np.linalg.norm(truth[1] - truth[0])            # the chain's own spacing
    assert np.allclose(fa.joins, step, atol=1e-3) and fa.joined


def test_the_fill_follows_the_deposit():
    dep, pred, _ = _pair()
    model = fill_structure(dep, "gaps", prediction=pred)
    rot = rotation_matrix(np.array([0.0, 1.0, 0.0]), 0.7)
    moved = dep.xyz.astype(float) @ rot.T + [3.0, 1.0, -2.0]
    assert np.allclose(model.place(moved), model.atoms.xyz.astype(float) @ rot.T
                       + [3.0, 1.0, -2.0], atol=1e-3)


def test_a_terminus_is_filled_only_when_asked():
    dep, pred, truth = _pair()
    dep = dep.subset(dep.res_seq <= 50)                    # C-terminal 51-60 unresolved
    assert {f.stretch.kind for f in fill_structure(dep, "gaps", prediction=pred).fills} == {"gap"}
    full = fill_structure(dep, "full", prediction=pred)
    ct = [f for f in full.fills if f.stretch.kind == "c_term"]
    assert len(ct) == 2 and all(len(f.joins) == 1 for f in ct)
    got = placed_ca(full, ct[0])
    assert max(np.linalg.norm(got[r] - truth[r - 1]) for r in got) < 1e-3


def test_anchors_must_match_in_sequence():
    """A prediction numbered one off is refused, not fitted to the wrong residues."""
    dep, pred, _ = _pair()
    shifted = pred.copy_with_coords(pred.xyz)
    shifted.res_seq = pred.res_seq + 1
    model = fill_stretches(dep, shifted, [Stretch("A", 25, 34, "gap")])
    assert not model.fills and "anchors" in model.skipped[0].reason


def test_prediction_for_refuses_a_foreign_numbering(tmp_path):
    dep, pred, _ = _pair()
    good, bad = tmp_path / "AF-GOOD-F1-model_v1.pdb", tmp_path / "AF-BAD-F1-model_v1.pdb"
    pred.to_pdb(good)
    shifted = pred.copy_with_coords(pred.xyz)
    shifted.res_seq = pred.res_seq + 1
    shifted.to_pdb(bad)
    chosen, ident = prediction_for(dep, [bad, good])
    assert chosen.name == "AF-GOOD-F1" and ident == 1.0
    with pytest.raises(GraftRefusal, match="AF-BAD-F1"):
        prediction_for(dep, [bad])


def test_stubs_are_not_anchors():
    dep, pred, _ = _pair(stub=range(19, 25))               # every anchor before the gap
    model = fill_stretches(dep, pred, [Stretch("A", 25, 34, "gap")])
    assert not model.fills and model.skipped[0].reason.startswith("anchors 0 before")


def test_clashes_count_other_subunits():
    dep, pred, truth = _pair()
    clean = fill_structure(dep, "gaps", prediction=pred)
    assert all(f.clashes == 0 for f in clean.fills)
    # Chain B moved so that one of its atoms sits on hidden residue 30 of chain A.
    xyz = dep.xyz.copy()
    b = np.flatnonzero((dep.chain == "B") & (dep.res_seq == 5) & (dep.atom_name == "CA"))[0]
    xyz[b] = truth[29] + [0.5, 0.0, 0.0]
    hit = fill_stretches(dep.copy_with_coords(xyz), pred, [Stretch("A", 25, 34, "gap")])
    assert hit.fills[0].clashes == 1


def test_a_broken_seam_is_reported():
    dep, pred, _ = _pair()
    bent = pred.xyz.copy()
    bent[(pred.res_seq >= 25) & (pred.res_seq <= 34)] += [0.0, 0.0, 9.0]
    model = fill_stretches(dep, pred.copy_with_coords(bent), [Stretch("A", 25, 34, "gap")])
    assert not model.fills[0].joined


def test_calibration_scores_a_right_and_a_wrong_prediction():
    dep, pred, _ = _pair(hidden=(0, -1))                   # nothing hidden: the host
    host = dep
    good = trial(host, pred, "A", 25, 34)
    assert good.rmsd_fill < 1e-3 < good.rmsd_line and good.beats_line
    wrong = pred.xyz.copy()
    inside = (pred.res_seq >= 25) & (pred.res_seq <= 34)
    wrong[inside] = wrong[inside][::-1]                     # the loop run backwards
    bad = trial(host, pred.copy_with_coords(wrong), "A", 25, 34)
    assert bad.rmsd_fill > good.rmsd_fill + 1.0


# ---------------------------------------------------------------- real data

@needs_structure("8TKG")
def test_8tkg_fills_every_internal_gap():
    from ip3r.io.loader import load
    from ip3r.io.predictions import local_predictions
    if not local_predictions():
        pytest.skip("no AlphaFold model downloaded")
    model = fill_structure(load("8TKG"), "gaps")
    assert model.prediction == "AF-Q14573-F1" and model.identity > 0.99
    assert len(model.fills) == 56 and not model.skipped
    # AlphaFold is least sure exactly where the map is empty.
    assert np.mean([f.plddt for f in model.fills]) < 50


@needs_structure("9YKK", "7LHF")
def test_deposits_in_no_downloaded_numbering_are_refused():
    from ip3r.io.loader import load
    from ip3r.io.predictions import local_predictions
    if not local_predictions():
        pytest.skip("no AlphaFold model downloaded")
    for pdb in ("9YKK", "7LHF"):
        with pytest.raises(GraftRefusal, match="no AlphaFold model is in its numbering"):
            fill_structure(load(pdb), "gaps")


@needs_structure("8TKG", "6DQN", "7T3P", "8TKF", "6DQJ", "8TKH", "8TLA")
def test_fills_beat_a_straight_line_on_stretches_8tkg_resolves():
    from ip3r.io.loader import load
    from ip3r.io.predictions import local_predictions
    from ip3r.structure.graft_calibration import calibrate
    if not local_predictions():
        pytest.skip("no AlphaFold model downloaded")
    host = load("8TKG")
    pred, _ = prediction_for(host)
    others = [load(i) for i in ("6DQN", "6DQJ", "7T3P", "8TKF", "8TKH", "8TLA")]
    trials = calibrate(host, others, pred)
    assert len(trials) == 16
    fill = np.array([t.rmsd_fill for t in trials])
    assert np.median(fill) < 1.5 < 5.0 < np.median([t.rmsd_line for t in trials])
    assert sum(t.beats_line for t in trials) >= 15
    # The fill is local: the global superposition is far worse.
    assert np.median([t.rmsd_global for t in trials]) > 5.0
    # pLDDT says which fills to believe.
    from ip3r.analysis.stats import spearman
    assert spearman([t.plddt for t in trials], fill)["rho"] < -0.5
