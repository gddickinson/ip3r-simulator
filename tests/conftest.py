"""Shared fixtures. Real-data tests skip (never fail) when the data is absent."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ip3r.config import genes_results  # noqa: E402
from ip3r.io.loader import is_local  # noqa: E402

needs_genes = pytest.mark.skipif(not genes_results().is_dir(),
                                 reason="ip3r_genes checkout not found")


def needs_structure(*ids):
    missing = [i for i in ids if not is_local(i)]
    return pytest.mark.skipif(bool(missing), reason=f"not downloaded: {missing}")


def c4_tetramer(n_per=60, seed=0, noise=0.0):
    """A synthetic C4 tetramer: one random chain turned by k x 90° about z."""
    rng = np.random.default_rng(seed)
    # A chain-like random walk well off the axis, so the four copies do not overlap.
    steps = rng.normal(size=(n_per, 3)) * 1.5
    chain = np.cumsum(steps, axis=0) + np.array([9.0, 0.0, 0.0])
    blocks = []
    for k in range(4):
        t = np.pi / 2 * k
        r = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
        blocks.append(chain @ r.T + rng.normal(size=chain.shape) * noise)
    return blocks


def make_structure(blocks, chains="ABCD", resname="ALA", atom="CA", element="C"):
    """A Structure of CA atoms, one chain per block, residues numbered from 1."""
    from ip3r.core.structure import Structure
    xyz = np.vstack(blocks).astype(np.float32)
    n = len(xyz)
    per = len(blocks[0])
    st = Structure(
        xyz=xyz, element=np.full(n, element, "U2"), atom_name=np.full(n, atom, "U6"),
        res_name=np.full(n, resname, "U5"),
        res_seq=np.tile(np.arange(1, per + 1, dtype=np.int32), len(blocks)),
        chain=np.repeat(np.array(list(chains[:len(blocks)]), "U6"), per),
        hetero=np.zeros(n, bool), b_factor=np.zeros(n, np.float32),
        occupancy=np.ones(n, np.float32), alt_loc=np.full(n, ".", "U2"),
        entity=np.full(n, "1", "U6"), name="SYN")
    st._build_residue_index()
    return st
