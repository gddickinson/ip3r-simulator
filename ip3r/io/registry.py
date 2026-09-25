"""The structure registry: which depositions the application knows about.

Read from ``resources/structures.json``, which ``scripts/sync_genes.py``
builds from ``ip3r_genes`` — 6DQN from S0's measurement record, the ITPR
references and state panel S11 selected by seven stated rules, and the six
IP3-bound entries S22 measured ligand shells in. None of it is typed here.

The RyR1 state panel is appended from ``resources/ryr1.json``, which
``scripts/curate_ryr.py`` selects from the PDB by stated rules (no RyR
structure is in ``ip3r_genes``). Open-state controls not in ip3r_genes'
panel (7T3T, Round 7.6) come from ``resources/ip3r_controls.json`` with the
role ``open_control``; the state panel leaves them out.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from ..config import RESOURCE_DIR, STRUCTURE_DIR

__all__ = ["StructureEntry", "load_registry", "get_entry", "local_path"]


@dataclass(frozen=True)
class StructureEntry:
    pdb_id: str
    paralog: str
    organism: str
    uniprot: str
    resolution: float
    state: str
    title: str
    roles: tuple = field(default_factory=tuple)
    ip3_bound: bool = False

    @property
    def label(self) -> str:
        return (f"{self.pdb_id} — {self.paralog} "
                f"({self.organism.split()[0][0]}. {self.organism.split()[-1]}), "
                f"{self.state}, {self.resolution:g} Å")

    @property
    def human(self) -> bool:
        return self.organism == "Homo sapiens"

    @property
    def is_control(self) -> bool:
        """Added here as a control, not part of ip3r_genes' selection."""
        return "open_control" in self.roles

    @property
    def family(self) -> str:
        return "RyR" if self.paralog.startswith("RYR") else "IP3R"


@lru_cache(maxsize=1)
def load_registry() -> tuple[StructureEntry, ...]:
    raw = json.loads((RESOURCE_DIR / "structures.json").read_text())["structures"]
    for extra in ("ryr1.json", "ip3r_controls.json"):
        path = RESOURCE_DIR / extra
        if path.exists():
            raw = raw + json.loads(path.read_text())["structures"]
    return tuple(StructureEntry(
        pdb_id=e["pdb_id"], paralog=e["paralog"], organism=e["organism"],
        uniprot=e.get("uniprot", ""), resolution=float(e["resolution"]),
        state=e["state"], title=e["title"], roles=tuple(e.get("roles", ())),
        ip3_bound=bool(e.get("ip3_bound"))) for e in raw)


def get_entry(pdb_id: str) -> StructureEntry | None:
    pdb_id = pdb_id.upper()
    return next((e for e in load_registry() if e.pdb_id == pdb_id), None)


def local_path(pdb_id: str, directory: Path | None = None) -> Path:
    """Where a fetched mmCIF for ``pdb_id`` lives (it may not exist yet)."""
    return (directory or STRUCTURE_DIR) / f"{pdb_id.upper()}.cif.gz"
