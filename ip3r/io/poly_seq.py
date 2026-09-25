"""The deposited construct's full sequence, observed or not.

An mmCIF's ``_pdbx_poly_seq_scheme`` lists every residue of every polymer
strand in the construct, with its author number (``pdb_seq_num``), including
the residues no atom was built for. The atoms alone cannot give this: an
unresolved stretch has no residue names. Filling a deposit through an
alignment needs it, because the alignment must see the residues that are to
be filled.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from .cif_reader import parse_cif_categories

__all__ = ["construct_sequences", "construct_sequence"]

_SCHEME = "pdbx_poly_seq_scheme"


@lru_cache(maxsize=8)
def construct_sequences(path: str) -> dict[str, dict[int, str]]:
    """``strand -> {author residue number: three-letter name}`` for one file.

    Empty when the file has no poly-seq scheme (a PDB-format file, or a
    predicted model).
    """
    p = Path(path)
    if not p.exists() or ".cif" not in p.name:
        return {}
    cat = parse_cif_categories(p, {_SCHEME}).get(_SCHEME)
    if not cat or "pdb_strand_id" not in cat:
        return {}
    out: dict[str, dict[int, str]] = {}
    for strand, num, name in zip(cat["pdb_strand_id"], cat["pdb_seq_num"], cat["mon_id"]):
        try:
            out.setdefault(strand, {})[int(num)] = name
        except ValueError:
            continue
    return out


def construct_sequence(path: str, strand: str) -> dict[int, str]:
    """One strand's construct, ``{author number: three-letter name}`` (empty if absent)."""
    return construct_sequences(str(path)).get(strand, {})
