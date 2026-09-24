"""AlphaFold DB models, downloaded into ``ref/alphafold``.

The model for each accession is **discovered from the AlphaFold DB API**, and
its entry id and version are never guessed. What the database holds for this
family is not what one would assume (queried 2026-09-24):

* ITPR3 (Q14573): the canonical sequence, ``AF-Q14573-F1`` (2,671 aa).
* ITPR1 (Q14643): only isoform 4 (``AF-Q14643-4-F1``, 2,695 aa), not the
  canonical 2,758-residue sequence the S17 tables number by.
* ITPR2 (Q14571): only isoform 2, 181 aa. It cannot fill a receptor, so it
  is not requested (:data:`UNAVAILABLE`).
* rat ITPR1 (P29994): only isoform 8 (2,695 aa).

So a model is never matched to a deposit by name. :func:`ip3r.structure.graft.
prediction_for` measures which downloaded model is in the deposit's own
numbering, by the same identity rule the variant painting uses.
"""

from __future__ import annotations

import gzip
import json
import shutil
import urllib.request
from functools import lru_cache
from pathlib import Path

from ..config import REF_DIR
from ..core.structure import Structure

__all__ = ["PREDICTION_DIR", "ACCESSIONS", "UNAVAILABLE", "fetch_prediction",
           "fetch_predictions", "local_predictions", "load_prediction"]

PREDICTION_DIR = REF_DIR / "alphafold"

#: Accessions whose AlphaFold DB model is downloaded, with what it is.
ACCESSIONS = {
    "Q14573": "human ITPR3, canonical",
    "Q14643": "human ITPR1 (AlphaFold DB holds isoform 4 only)",
    "P29994": "rat ITPR1 (AlphaFold DB holds isoform 8 only)",
}

#: What AlphaFold DB does not usefully have, and why. Kept rather than left as
#: a silent gap in the list above.
UNAVAILABLE = {
    "Q14571": "human ITPR2: AlphaFold DB models only isoform 2, 181 residues",
}

API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
_TIMEOUT_S = 120


def _api(acc: str) -> list[dict]:
    with urllib.request.urlopen(API_URL.format(acc=acc), timeout=_TIMEOUT_S) as r:
        data = json.load(r)
    return data if isinstance(data, list) else []


def _is_model(path: Path) -> bool:
    """A gzip that decompresses to an mmCIF with an ``_atom_site`` loop."""
    try:
        with gzip.open(path, "rt", errors="replace") as fh:
            return any(line.startswith("_atom_site.") for line in fh)
    except (OSError, EOFError):
        return False


def local_predictions(directory: Path | None = None) -> list[Path]:
    """Every downloaded model, sorted by name."""
    d = directory or PREDICTION_DIR
    return sorted(d.glob("AF-*-model_v*.cif.gz")) if d.exists() else []


def fetch_prediction(acc: str, directory: Path | None = None, log=print) -> Path:
    """Download the model AlphaFold DB currently serves for ``acc``."""
    d = directory or PREDICTION_DIR
    d.mkdir(parents=True, exist_ok=True)
    entries = _api(acc)
    if not entries:
        raise OSError(f"{acc}: AlphaFold DB has no entry")
    e = entries[0]
    entry_id = e.get("modelEntityId") or e.get("entryId")
    version = e["latestVersion"]
    dest = d / f"{entry_id}-model_v{version}.cif.gz"
    if dest.exists() and _is_model(dest):
        return dest
    url = e.get("cifUrl") or (f"https://alphafold.ebi.ac.uk/files/"
                              f"{entry_id}-model_v{version}.cif")
    log(f"  {acc}: downloading {url}")
    tmp = dest.with_suffix(".part")
    with urllib.request.urlopen(url, timeout=_TIMEOUT_S) as resp, \
            gzip.open(tmp, "wb") as fh:
        shutil.copyfileobj(resp, fh)
    if not _is_model(tmp):
        tmp.unlink(missing_ok=True)
        raise OSError(f"{acc}: downloaded file is not a valid mmCIF")
    tmp.replace(dest)
    return dest


def fetch_predictions(directory: Path | None = None, log=print) -> dict[str, str]:
    """Fetch every model in :data:`ACCESSIONS`. Returns accession -> status."""
    status = {}
    for acc in ACCESSIONS:
        try:
            status[acc] = fetch_prediction(acc, directory, log=log).name
        except OSError as exc:
            status[acc] = f"failed: {exc}"
            log(f"  {acc}: {exc}")
    return status


@lru_cache(maxsize=4)
def load_prediction(path: str) -> Structure:
    """Parse one model (memoised); pLDDT is in ``b_factor``."""
    p = Path(path)
    st = Structure.from_file(p, name=p.name.split("-model")[0])
    st.meta["prediction"] = p.name
    return st
