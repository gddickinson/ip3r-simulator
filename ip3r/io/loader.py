"""Load a registry structure by PDB id, from ``ref/`` or (optionally) RCSB.

Parsed structures are memoised per process, because a 139,000-atom tetramer
takes a few seconds to read and several checks and panels ask for the same
one. Whether a missing file may be downloaded is a process-wide switch
(:data:`ALLOW_FETCH`), off by default so that running the checks never
touches the network unless the user asked it to (``--fetch``).
"""

from __future__ import annotations

from functools import lru_cache

from ..core.structure import Structure
from .fetch import fetch_structure
from .registry import local_path

__all__ = ["load", "is_local", "ALLOW_FETCH", "StructureUnavailable"]

#: Set True (by ``python -m ip3r ... --fetch`` or the GUI) to allow downloads.
ALLOW_FETCH = False


class StructureUnavailable(FileNotFoundError):
    """The structure is not in ``ref/`` and fetching is not allowed."""


def is_local(pdb_id: str) -> bool:
    return local_path(pdb_id).exists()


@lru_cache(maxsize=8)
def load(pdb_id: str) -> Structure:
    pdb_id = pdb_id.upper()
    path = local_path(pdb_id)
    if not path.exists():
        if not ALLOW_FETCH:
            raise StructureUnavailable(
                f"{pdb_id} is not downloaded — run `python -m ip3r fetch` "
                f"or pass --fetch")
        path = fetch_structure(pdb_id)
    st = Structure.from_file(path, name=pdb_id)
    st.meta["pdb_id"] = pdb_id
    return st
