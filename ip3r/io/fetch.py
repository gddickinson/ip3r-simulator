"""Download deposited structures into ``ref/structures``.

Files are fetched from RCSB as gzipped mmCIF and kept gzipped; the reader
handles ``.cif.gz`` directly. A file that already exists and decompresses is
never re-downloaded, so ``python -m ip3r fetch`` is idempotent and a second
run is offline.

If the ``ip3r_genes`` data root already holds a copy (its S11 stage
downloaded the reference panel), set ``IP3R_STRUCTURE_MIRROR`` to that
directory and matching files are copied rather than downloaded.
"""

from __future__ import annotations

import gzip
import os
import shutil
import urllib.request
from pathlib import Path

from ..config import STRUCTURE_DIR, ensure_dirs
from .registry import load_registry, local_path

__all__ = ["fetch_structure", "fetch_all", "is_valid"]

RCSB_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"
_TIMEOUT_S = 120


def is_valid(path: Path) -> bool:
    """A gzip that decompresses and contains an ``_atom_site`` loop.

    Scans to the loop however far into the file it is. The first version gave
    up after 20,000 lines and rejected 7LHF and 8TKG, whose headers (struct_conn
    and friends for a 10,000-residue tetramer) are longer than that — a
    validity check that fails on valid files, found by its first real run.
    """
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith("_atom_site."):
                    return True
    except (OSError, EOFError):
        return False
    return False


def _from_mirror(pdb: str, dest: Path) -> bool:
    mirror = os.environ.get("IP3R_STRUCTURE_MIRROR")
    if not mirror:
        return False
    for name in (f"{pdb.lower()}.cif", f"{pdb.upper()}.cif"):
        src = Path(mirror) / name
        if src.exists():
            with open(src, "rb") as fin, gzip.open(dest, "wb") as fout:
                shutil.copyfileobj(fin, fout)
            return is_valid(dest)
    return False


def fetch_structure(pdb_id: str, directory: Path | None = None,
                    force: bool = False, log=print) -> Path:
    """Return the local path of ``pdb_id``, downloading it if needed."""
    ensure_dirs()
    dest = local_path(pdb_id, directory)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force and is_valid(dest):
        return dest
    if _from_mirror(pdb_id, dest):
        log(f"  {pdb_id}: copied from IP3R_STRUCTURE_MIRROR")
        return dest
    url = RCSB_URL.format(pdb=pdb_id.upper())
    tmp = dest.with_suffix(".part")
    log(f"  {pdb_id}: downloading {url}")
    with urllib.request.urlopen(url, timeout=_TIMEOUT_S) as resp, \
            open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    if not is_valid(tmp):
        tmp.unlink(missing_ok=True)
        raise OSError(f"{pdb_id}: downloaded file is not a valid mmCIF")
    tmp.replace(dest)
    return dest


def fetch_all(pdb_ids: list[str] | None = None, log=print) -> dict[str, str]:
    """Fetch every registry entry (or the listed ones). Returns id -> status."""
    ids = pdb_ids or [e.pdb_id for e in load_registry()]
    status = {}
    for pdb in ids:
        try:
            fetch_structure(pdb, STRUCTURE_DIR, log=log)
            status[pdb] = "ok"
        except OSError as exc:
            status[pdb] = f"failed: {exc}"
            log(f"  {pdb}: {exc}")
    return status
