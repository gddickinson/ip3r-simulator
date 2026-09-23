#!/usr/bin/env python
"""Import the curated annotation this application needs from ``ip3r_genes``.

Four resources are written to ``ip3r/resources/``, each carrying the path and
SHA-256 of every source table it was built from and the ``ip3r_genes`` commit
at import time:

* ``sequences.json`` — the three human reference sequences, read residue by
  residue out of the S17 constraint tables (so a residue number here is the
  number those tables use, by construction).
* ``domains.json`` — the S17 domain map: Pfam elements measured per accession
  and the 6DQN structural elements transferred with an anchor test.
* ``sites.json`` — the measured functional residues (ten IP3 contacts, two
  filter- and two gate-lining residues) in each paralog's own numbering.
* ``constraint.json`` — per-residue conservation (JSD) of each human paralog
  in the four S17 layers, ``null`` where the layer is unreliable or below
  S17's occupancy floor, so the viewer can colour without ip3r_genes.
* ``variants.json`` — the S17 variant harvest, one row per variant, with
  its class bucket, so variants can be placed on structures.
* ``structures.json`` — the structure registry: 6DQN from S0's measurement
  record, every ITPR reference and state-panel entry S11 selected, and every
  IP3-bound deposition S22 measured its ligand shells in. Nothing is typed.

Why import rather than read live: the application must run, and its tests
must pass, on a machine without the publication project. The findings checks
are the exception — they read ``ip3r_genes/results`` live, because confirming a
result means reading the table it was written to.

Usage::

    python scripts/sync_genes.py            # rebuild the resources
    python scripts/sync_genes.py --check    # exit 1 if any source has drifted
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ip3r.config import GENES_DIR, PARALOG_ACC, RESOURCE_DIR  # noqa: E402

RESULTS = GENES_DIR / "results"
SOURCES = {
    "domain_map": "constraint/domain_map.tsv",
    "functional_sites": "constraint/functional_sites.tsv",
    "structure_meta": "s0_baseline/review_figures/structure_meta.json",
    "reference_selection": "structures/reference_selection.tsv",
    "shell_agreement": "ligand_site/shell_agreement.tsv",
    "variants": "constraint/variants.tsv",
    **{f"constraint_{g}": f"constraint/constraint_{g}_{a}.tsv"
       for g, a in PARALOG_ACC.items()},
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(rel: str) -> list[dict]:
    with open(RESULTS / rel, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genes_commit() -> str:
    try:
        return subprocess.run(["git", "-C", str(GENES_DIR), "rev-parse",
                               "--short", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def provenance(keys: list[str]) -> dict:
    return {"genes_commit": genes_commit(),
            "sources": {SOURCES[k]: sha256(RESULTS / SOURCES[k]) for k in keys}}


# ----------------------------------------------------------------- builders

def build_sequences() -> dict:
    out = {}
    for gene, acc in PARALOG_ACC.items():
        rows = read_tsv(SOURCES[f"constraint_{gene}"])
        resi = [int(r["resi"]) for r in rows]
        if resi != list(range(1, len(rows) + 1)):
            raise SystemExit(f"{gene}: constraint table is not residues 1..N")
        out[gene] = {"accession": acc, "length": len(rows),
                     "sequence": "".join(r["aa"] for r in rows)}
    return {"provenance": provenance([f"constraint_{g}" for g in PARALOG_ACC]),
            "paralogs": out}


def build_domains() -> dict:
    out: dict[str, list] = {g: [] for g in PARALOG_ACC}
    for r in read_tsv(SOURCES["domain_map"]):
        out[r["paralog"]].append({
            "element": r["element"], "kind": r["kind"], "source": r["source"],
            "start": int(r["start"]), "end": int(r["end"])})
    return {"provenance": provenance(["domain_map"]), "paralogs": out}


def build_sites() -> dict:
    out: dict[str, list] = {g: [] for g in PARALOG_ACC}
    for r in read_tsv(SOURCES["functional_sites"]):
        out[r["paralog"]].append({
            "site_class": r["site_class"], "resi": int(r["resi"]),
            "aa": r["aa"], "source_resi": int(r["source_resi"]),
            "transferred": r["transferred"] == "True"})
    return {"provenance": provenance(["functional_sites"]), "paralogs": out}


def build_structures() -> dict:
    meta = json.loads((RESULTS / SOURCES["structure_meta"]).read_text())
    ip3_bound = {r["pdb_id"] for r in read_tsv(SOURCES["shell_agreement"])}
    entries = {meta["pdb_id"]: {
        "pdb_id": meta["pdb_id"], "paralog": "ITPR3",
        "organism": "Homo sapiens", "uniprot": meta["uniprot"],
        "resolution": meta["resolution_A"], "state": "IP3-bound",
        "title": meta["description"], "roles": ["s0_measurement"],
        "ip3_bound": meta["pdb_id"] in ip3_bound}}
    for r in read_tsv(SOURCES["reference_selection"]):
        if r["call"] != "ITPR":
            continue
        e = entries.setdefault(r["entry_id"], {
            "pdb_id": r["entry_id"], "paralog": r["paralog"],
            "organism": r["organism"],
            "uniprot": PARALOG_ACC.get(r["paralog"], "")
            if r["organism"] == "Homo sapiens" else "",
            "resolution": float(r["resolution"]), "state": r["state"],
            "title": r["title"], "roles": [],
            "ip3_bound": r["entry_id"] in ip3_bound})
        e["roles"].append(r["role"])
    missing = ip3_bound - set(entries)
    if missing:
        raise SystemExit(f"IP3-bound entries with no registry row: {missing}")
    return {"provenance": provenance(["structure_meta", "reference_selection",
                                      "shell_agreement"]),
            "structures": sorted(entries.values(),
                                 key=lambda e: (e["paralog"], e["pdb_id"]))}


LAYERS = ("deep", "shallow", "vert", "family")
MIN_OCCUPANCY = 0.5          # S17's floor; mirrored by constraint.min_occupancy


def _score(row: dict, layer: str):
    try:
        v = float(row[f"{layer}_jsd"])
        occ = float(row[f"{layer}_occupancy"] or "nan")
    except ValueError:
        return None
    if occ == occ and occ < MIN_OCCUPANCY:
        return None
    if layer == "deep" and row["deep_reliable"] != "True":
        return None
    return round(v, 4)


def build_constraint() -> dict:
    out = {}
    for gene in PARALOG_ACC:
        rows = read_tsv(SOURCES[f"constraint_{gene}"])
        out[gene] = {la: [_score(r, la) for r in rows] for la in LAYERS}
    return {"provenance": provenance([f"constraint_{g}" for g in PARALOG_ACC]),
            "layers": list(LAYERS), "min_occupancy": MIN_OCCUPANCY,
            "paralogs": out}


def build_variants() -> dict:
    keep = ("gene", "source", "variation_id", "resi", "ref_aa", "alt_aa",
            "classification", "class_bucket", "condition", "review_status")
    rows = [{k: r[k] for k in keep} for r in read_tsv(SOURCES["variants"])]
    for r in rows:
        r["resi"] = int(r["resi"])
    return {"provenance": provenance(["variants"]), "variants": rows}


BUILDERS = {"constraint.json": build_constraint,
            "variants.json": build_variants,"sequences.json": build_sequences, "domains.json": build_domains,
            "sites.json": build_sites, "structures.json": build_structures}


def check() -> int:
    """Report every resource whose recorded source hash no longer matches."""
    drift = 0
    for name in BUILDERS:
        path = RESOURCE_DIR / name
        if not path.exists():
            print(f"  {name}: missing — run scripts/sync_genes.py")
            drift += 1
            continue
        for rel, digest in json.loads(path.read_text())["provenance"][
                "sources"].items():
            src = RESULTS / rel
            if not src.exists():
                print(f"  {name}: source {rel} no longer exists")
                drift += 1
            elif sha256(src) != digest:
                print(f"  {name}: source {rel} has changed since import")
                drift += 1
    print("resources in step with ip3r_genes" if not drift
          else f"{drift} drifted source(s)")
    return 1 if drift else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    if not RESULTS.exists():
        raise SystemExit(f"ip3r_genes results not found at {RESULTS}; set "
                         f"IP3R_GENES_DIR")
    if args.check:
        return check()
    for name, build in BUILDERS.items():
        data = build()
        (RESOURCE_DIR / name).write_text(json.dumps(data, indent=1))
        print(f"wrote {name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
