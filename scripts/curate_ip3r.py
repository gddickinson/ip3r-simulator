#!/usr/bin/env python
"""Curate the extended IP3R deposit set: every full-length EM tetramer.

``ip3r_genes`` S11 chose one deposit per ITPR3 state, and S22 the six it
measured ligand shells in; that selection stays the publication's and is
imported by ``sync_genes.py``. This script adds the rest of the PDB's
full-length IP3R structures so they can be viewed, measured, superposed
and morphed, by stated rules and with nothing typed: it searches, filters
and writes ``ip3r/resources/ip3r_extended.json`` with the URL and SHA-256 of
every response it used and the date it ran. Entries get the role
``extended``: the state panel, the shortfall panel and every findings check
leave them out, so no publication number moves.

**Rules** (every PDB entry mapped to a human ITPR, or rat ITPR1 P29994,
the species of the reference 7LHF):

1. single-particle EM, four chains of one ITPR accession, each sample at
   least ``MIN_CHAIN_LENGTH`` residues (full-length, not a deletion
   construct), and at least ``MIN_MODELLED`` residues modelled over the
   entry (not a partial model or a crystal of a domain);
2. wild type: no ``pdbx_mutation``, no "mutant" in the title;
3. resolution at most ``MAX_RESOLUTION`` A, so side chains are placed
   (the RyR1 panel's bar);
4. not already registered (``structures.json`` or ``ip3r_controls.json``
   keep their roles).

The state label is read from the title (``state_of``), never assigned.

Usage::

    python scripts/curate_ip3r.py            # rebuild the resource
    python scripts/curate_ip3r.py --dry-run  # print the selection only
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ip3r.config import PARALOG_ACC, RESOURCE_DIR  # noqa: E402
from scripts.curate_ryr import GRAPHQL, SEARCH, _SOURCES, _get  # noqa: E402

#: accession -> (paralog, organism). Rat ITPR1 because 7LHF is rat.
ACCESSIONS = {acc: (name, "Homo sapiens") for name, acc in PARALOG_ACC.items()}
ACCESSIONS["P29994"] = ("ITPR1", "Rattus norvegicus")

MIN_CHAIN_LENGTH = 2600       # ITPR3 2671, ITPR1 2749, ITPR2 2701; 7T3* 2633
MIN_MODELLED = 7000           # a tetramer's channel and most of its cytosol
MAX_RESOLUTION = 4.0

#: Title fragment -> state. Checked in order: "inactive" before "active",
#: "preactivated+Ca" before "preactivated", "labile resting" before "resting".
STATE_WORDS = (("inactive-like", "inactive-like"), ("inactive", "inactive"),
               ("inhibited", "inhibited"), ("preactivated+ca", "preactivated + Ca2+"),
               ("pre-active", "preactivated"), ("preactivated", "preactivated"),
               ("activated", "activated"), ("active state", "active"),
               ("labile resting", "labile resting"), ("resting", "resting"),
               ("apo", "apo"), ("ligand-free", "apo"),
               ("calcium/ip3/atp", "Ca2+/IP3/ATP"), ("high ca", "high Ca2+"),
               ("low ip3 ca", "low IP3 + Ca2+"), ("high ip3 ca", "high IP3 + Ca2+"),
               ("ca2+-bound", "Ca2+-bound"), ("ip3-bound", "IP3-bound"))
#: What tells two deposits of one state apart, kept beside it.
QUALIFIERS = (r"\b(class \d)", r"pre-active ([a-c]) state", r"state (\d)\b",
              r"(nanodisc)")

GQL = """{ entries(entry_ids: %s) { rcsb_id struct { title } exptl { method }
 rcsb_entry_info { resolution_combined deposited_modeled_polymer_monomer_count }
 rcsb_primary_citation { pdbx_database_id_DOI }
 nonpolymer_entities { nonpolymer_comp { chem_comp { id } } }
 polymer_entities { rcsb_polymer_entity { pdbx_mutation pdbx_number_of_molecules }
  rcsb_polymer_entity_container_identifiers { reference_sequence_identifiers {
   database_accession } }
  entity_poly { rcsb_sample_sequence_length } } } }"""


def entry_ids() -> list[str]:
    ids: set[str] = set()
    for acc in ACCESSIONS:
        q = {"query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers."
                         "reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": acc}},
            "return_type": "entry", "request_options": {"return_all_hits": True}}
        body = _get(SEARCH, q)        # an accession with no entries answers empty
        if body:
            ids |= {h["identifier"] for h in json.loads(body)["result_set"]}
    return sorted(ids)


def state_of(title: str) -> str | None:
    t = title.lower()
    base = next((label for word, label in STATE_WORDS if word in t), None)
    if base is None:
        return None
    for pattern in QUALIFIERS:
        m = re.search(pattern, t)
        if m:
            q = m.group(1)
            return f"{base} ({q.upper() if len(q) == 1 and q.isalpha() else q})"
    return base


def accession_of(e: dict) -> tuple[str | None, list[dict]]:
    """The one ITPR accession of an entry and its polymer entities."""
    by_acc: dict[str, list[dict]] = {}
    for p in e["polymer_entities"] or []:
        ids = p["rcsb_polymer_entity_container_identifiers"]["reference_sequence_identifiers"] or []
        for x in ids:
            if x["database_accession"] in ACCESSIONS:
                by_acc.setdefault(x["database_accession"], []).append(p)
    if len(by_acc) != 1:
        return None, []
    acc, ents = next(iter(by_acc.items()))
    return acc, ents


def verdict(e: dict) -> str | None:
    """``None`` if the entry passes rules 1-3, else the first rule it fails."""
    info = e["rcsb_entry_info"]
    if e["exptl"][0]["method"] != "ELECTRON MICROSCOPY":
        return "1: not EM"
    acc, ents = accession_of(e)
    if acc is None:
        return "1: not one ITPR accession"
    n = sum(p["rcsb_polymer_entity"]["pdbx_number_of_molecules"] for p in ents)
    if n != 4 or any(p["entity_poly"]["rcsb_sample_sequence_length"] < MIN_CHAIN_LENGTH
                     for p in ents):
        return "1: not four full-length chains"
    if (info["deposited_modeled_polymer_monomer_count"] or 0) < MIN_MODELLED:
        return "1: too little modelled"
    if "mutant" in e["struct"]["title"].lower() or any(
            p["rcsb_polymer_entity"]["pdbx_mutation"] for p in ents):
        return "2: mutant"
    if (info["resolution_combined"] or [99.0])[0] > MAX_RESOLUTION:
        return "3: resolution"
    return None


def registered() -> set[str]:
    out = set()
    for name in ("structures.json", "ip3r_controls.json"):
        out |= {s["pdb_id"] for s in json.loads((RESOURCE_DIR / name).read_text())["structures"]}
    return out


def extended() -> tuple[list[dict], dict[str, str]]:
    known = registered()
    body = json.loads(_get(GRAPHQL, {"query": GQL % json.dumps(entry_ids())}))
    chosen, rejected = [], {}
    for e in body["data"]["entries"]:
        why = verdict(e)
        if why is None and e["rcsb_id"] in known:
            why = "4: already registered"
        state = state_of(e["struct"]["title"])
        if why is None and state is None:
            why = "state: none named in the title"
        if why:
            rejected[e["rcsb_id"]] = why
            continue
        acc, _ = accession_of(e)
        paralog, organism = ACCESSIONS[acc]
        ligands = {n["nonpolymer_comp"]["chem_comp"]["id"]
                   for n in e["nonpolymer_entities"] or []}
        chosen.append({
            "pdb_id": e["rcsb_id"], "paralog": paralog, "organism": organism,
            "uniprot": acc if organism == "Homo sapiens" else "",
            "resolution": e["rcsb_entry_info"]["resolution_combined"][0],
            "state": state, "title": e["struct"]["title"],
            "roles": ["extended"], "ip3_bound": "I3P" in ligands,
            "doi": (e["rcsb_primary_citation"] or {}).get("pdbx_database_id_DOI", "")})
    return sorted(chosen, key=lambda r: (r["paralog"], r["organism"], r["pdb_id"])), rejected


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    structures, rejected = extended()
    for s in structures:
        print(f"{s['pdb_id']}  {s['paralog']} {s['organism'].split()[0][:3]}  "
              f"{s['state']:22s} {s['resolution']:.2f} A  {s['title']}")
    counts: dict[str, int] = {}
    for why in rejected.values():
        counts[why.split(":")[0]] = counts.get(why.split(":")[0], 0) + 1
    print(f"{len(structures)} chosen; {len(rejected)} rejected, by first failed rule: {counts}")
    if args.dry_run:
        return 0
    out = {"provenance": {"retrieved": datetime.date.today().isoformat(),
                          "sources": _SOURCES,
                          "rules": {"accessions": sorted(ACCESSIONS),
                                    "min_chain_length": MIN_CHAIN_LENGTH,
                                    "min_modelled": MIN_MODELLED,
                                    "max_resolution_A": MAX_RESOLUTION}},
           "structures": structures, "rejected": rejected}
    path = RESOURCE_DIR / "ip3r_extended.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
