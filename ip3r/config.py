"""Project-wide paths, reference identities and runtime settings.

Every other module imports paths from here rather than computing its own, so
relocating a data directory is a one-line change (or one environment
variable).

Three places data comes from, kept apart on purpose:

* ``ip3r/resources/`` — curated JSON **committed** with the package. Most of it
  is imported from the ``ip3r_genes`` publication project by
  ``scripts/sync_genes.py``, which records the SHA-256 of every source table so
  a resource that has drifted from its source is visible.
* ``ref/`` — downloads (mmCIF files). Git-ignored and regenerable with
  ``python -m ip3r fetch``.
* ``GENES_DIR`` — the ``ip3r_genes`` checkout itself. The findings checks read
  its committed ``results/`` tables directly, because confirming a result means
  reading the table the result was written to, not a copy of it.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Filesystem layout
# --------------------------------------------------------------------------

PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PACKAGE_DIR.parent

#: Curated, authored annotation shipped with the package.
RESOURCE_DIR = PACKAGE_DIR / "resources"

#: Downloads. Git-ignored, regenerable.
REF_DIR = Path(os.environ.get("IP3R_REF_DIR", PROJECT_ROOT / "ref"))
STRUCTURE_DIR = REF_DIR / "structures"

#: Derived artefacts (mode caches, exported figures). Git-ignored.
DATA_DIR = Path(os.environ.get("IP3R_DATA_DIR", PROJECT_ROOT / "data"))
CACHE_DIR = DATA_DIR / "cache"
DERIVED_DIR = DATA_DIR / "derived"

#: The publication project whose results this application illustrates and
#: re-checks. A sibling checkout by default.
GENES_DIR = Path(os.environ.get("IP3R_GENES_DIR",
                                PROJECT_ROOT.parent / "ip3r_genes"))


def genes_results() -> Path:
    """``<GENES_DIR>/results`` — resolved at call time so tests can redirect."""
    return Path(os.environ.get("IP3R_GENES_DIR", GENES_DIR)) / "results"


def ensure_dirs() -> None:
    """Create every directory the application writes to."""
    for d in (STRUCTURE_DIR, CACHE_DIR, DERIVED_DIR):
        d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Reference molecules
# --------------------------------------------------------------------------

#: Canonical human UniProt accessions. The ``ip3r_genes`` constraint tables are
#: keyed on exactly these, so a residue number here means the same residue as
#: a residue number there.
PARALOG_ACC = {"ITPR1": "Q14643", "ITPR2": "Q14571", "ITPR3": "Q14573"}
PARALOGS = tuple(PARALOG_ACC)

#: Ryanodine receptors, curated here (``scripts/curate_ryr.py``), not from
#: ``ip3r_genes``. Rabbit RyR1, because the deposits are rabbit: a residue
#: number on an RyR1 deposit is a P11716 number.
RYR_ACC = {"RYR1": "P11716"}

#: Every reference numbering a deposit can be in. ``PARALOGS`` stays the
#: three IP3Rs the publication's tables are keyed on; RyR1 has a sequence
#: and domains but no conservation, sites or variants, so those paint grey.
NUMBERINGS = PARALOGS + tuple(RYR_ACC)

#: IP3 receptors are homotetramers with C4 symmetry about the pore axis.
N_SUBUNITS = 4

#: The chemical component id of inositol 1,4,5-trisphosphate in the PDB.
IP3_COMP_ID = "I3P"

#: The structure the application opens with: human ITPR3 with IP3 bound, the
#: one ``ip3r_genes`` S0 measured the channel on.
DEFAULT_STRUCTURE = "6DQN"


@dataclass
class RenderSettings:
    """Defaults for the OpenGL viewport."""

    gl_major: int = 4
    gl_minor: int = 1
    samples: int = 4
    background: tuple[float, float, float, float] = (0.055, 0.063, 0.086, 1.0)
    ambient_occlusion: bool = True
    depth_cue: bool = True
    fov_degrees: float = 35.0
    target_fps: int = 60


@dataclass
class AppSettings:
    """Top-level runtime settings, mutable from the GUI."""

    render: RenderSettings = field(default_factory=RenderSettings)
    default_structure: str = DEFAULT_STRUCTURE


SETTINGS = AppSettings()
