#!/usr/bin/env python
"""Curate the RyR1 resource: reference sequence, domains and a state panel.

``ip3r_genes`` has no ryanodine-receptor structures, so this resource is
curated here, from public databases, by stated rules. Nothing is typed:
the script searches, filters and writes ``ip3r/resources/ryr1.json`` with
the URL and SHA-256 of every response it used and the date it ran.

**Reference.** Rabbit RyR1, UniProt P11716. The deposits are rabbit, so
their residue numbers are P11716's; the numbering check confirms it per
deposit (the same S24 rule as the ITPRs). Domains are InterPro's Pfam
matches on that accession. PF00520 is named ``channel``, as for the ITPRs,
because the pore-domain span is read from it.

**The state panel, by rule** (every entry of the PDB mapped to P11716):

1. single-particle EM, four P11716 chains of at least
   ``MIN_CHAIN_LENGTH`` residues in the sample, and at least
   ``MIN_MODELLED`` residues modelled: full-length tetramers, not domain
   crystals or local refinements (titles naming a local, focused or
   signal-subtracted refinement are also excluded);
2. wild type: no ``pdbx_mutation``, no "mutant" or "chimeric" in the title;
3. only the physiological or activating ligands in ``ALLOWED_LIGANDS`` (no
   drug, toxin or ryanodine: the panel asks what gating does to the pore);
4. resolution at most ``MAX_RESOLUTION`` A, so side chains are placed and
   can carry the wall charge in the conductance model;
5. the state is the one the depositors' title names (``STATE_WORDS``);
   per state, the best resolution wins;
6. the morph needs a start from the open deposit's own preparation, so
   the best-resolved passing primed deposit of the same paper (DOI) is
   added with the role ``morph_start``, if rule 5 did not already pick it.
   Gating seen across two preparations would mix in the preparation.

Usage::

    python scripts/curate_ryr.py            # rebuild the resource
    python scripts/curate_ryr.py --dry-run  # print the selection only
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ip3r.config import RESOURCE_DIR, RYR_ACC  # noqa: E402

ACC = RYR_ACC["RYR1"]
UNIPROT = f"https://rest.uniprot.org/uniprotkb/{ACC}.fasta"
INTERPRO = (f"https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/"
            f"{ACC}/?page_size=100")
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"

MIN_CHAIN_LENGTH = 5000       # the sample, not the model: full-length construct
MIN_MODELLED = 12000          # residues modelled over the whole entry
MAX_RESOLUTION = 4.0
ALLOWED_LIGANDS = {"CA", "ZN", "MG", "K", "NA", "CL", "ATP", "ADP", "ACP",
                   "CFF", "POV", "PCW", "PC1", "CLR", "LBN"}
EXCLUDED_TITLE = ("local refinement", "focus", "signal subtracted", "mutant",
                  "chimeric")
#: Title word -> state label. Checked in order, so "closed-inactivated"
#: is not read as "inactivated", nor "close" inside "closed" twice.
STATE_WORDS = (("closed-inactivated", "closed-inactivated"),
               ("inactivated", "inactivated"), ("primed", "primed"),
               ("open", "open"), ("closed", "closed"), ("close state", "closed"))

#: Pfam accession -> element name. Names shared with the ITPRs where the
#: domain is shared.
PFAM_NAMES = {"PF08709": "nterm_trefoil", "PF02815": "MIR", "PF01365": "RIH",
              "PF08454": "RIH_assoc", "PF00520": "channel", "PF00622": "SPRY",
              "PF02026": "RyR_repeat", "PF06459": "RyR_TM4_6",
              "PF21119": "junctional_solenoid"}

GQL = """{ entries(entry_ids: %s) { rcsb_id struct { title } exptl { method }
 rcsb_entry_info { resolution_combined deposited_modeled_polymer_monomer_count }
 rcsb_primary_citation { pdbx_database_id_DOI }
 nonpolymer_entities { nonpolymer_comp { chem_comp { id } } }
 polymer_entities { rcsb_polymer_entity { pdbx_mutation pdbx_number_of_molecules }
  rcsb_polymer_entity_container_identifiers { reference_sequence_identifiers {
   database_accession } }
  entity_poly { rcsb_sample_sequence_length } } } }"""

_SOURCES: dict[str, str] = {}


def _get(url: str, payload: dict | None = None) -> bytes:
    data = None if payload is None else json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data,
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        body = r.read()
    key = url if payload is None else f"{url} {hashlib.sha256(data).hexdigest()[:12]}"
    _SOURCES[key] = hashlib.sha256(body).hexdigest()
    return body


def sequence() -> str:
    lines = _get(UNIPROT).decode().splitlines()
    return "".join(x.strip() for x in lines if not x.startswith(">"))


def domains() -> list[dict]:
    out = []
    for r in json.loads(_get(INTERPRO))["results"]:
        acc = r["metadata"]["accession"]
        for p in r["proteins"]:
            for loc in p["entry_protein_locations"]:
                for f in loc["fragments"]:
                    out.append({"element": PFAM_NAMES.get(acc, acc), "kind": "pfam",
                                "source": f"InterPro {acc}",
                                "start": int(f["start"]), "end": int(f["end"])})
    return sorted(out, key=lambda e: e["start"])


def entry_ids() -> list[str]:
    q = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers."
                     "reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": ACC}},
        "return_type": "entry", "request_options": {"return_all_hits": True}}
    return sorted(h["identifier"] for h in json.loads(_get(SEARCH, q))["result_set"])


def state_of(title: str) -> str | None:
    t = title.lower()
    return next((label for word, label in STATE_WORDS if word in t), None)


def verdict(e: dict) -> str | None:
    """``None`` if the entry passes rules 1-4, else the first rule it fails."""
    title = e["struct"]["title"]
    info = e["rcsb_entry_info"]
    if e["exptl"][0]["method"] != "ELECTRON MICROSCOPY":
        return "1: not EM"
    ryr = [p for p in e["polymer_entities"] or []
           if ACC in [x["database_accession"] for x in
                      p["rcsb_polymer_entity_container_identifiers"]
                      ["reference_sequence_identifiers"] or []]]
    n = sum(p["rcsb_polymer_entity"]["pdbx_number_of_molecules"] for p in ryr)
    if n != 4 or any(p["entity_poly"]["rcsb_sample_sequence_length"] < MIN_CHAIN_LENGTH
                     for p in ryr):
        return "1: not four full-length chains"
    if (info["deposited_modeled_polymer_monomer_count"] or 0) < MIN_MODELLED:
        return "1: too little modelled"
    if any(w in title.lower() for w in EXCLUDED_TITLE):
        return "1/2: title excluded"
    if any(p["rcsb_polymer_entity"]["pdbx_mutation"] for p in ryr):
        return "2: mutant"
    ligands = {n["nonpolymer_comp"]["chem_comp"]["id"]
               for n in e["nonpolymer_entities"] or []}
    if ligands - ALLOWED_LIGANDS:
        return f"3: ligand {','.join(sorted(ligands - ALLOWED_LIGANDS))}"
    if (info["resolution_combined"] or [99.0])[0] > MAX_RESOLUTION:
        return "4: resolution"
    return None


def panel() -> tuple[list[dict], dict[str, str]]:
    ids = entry_ids()
    body = json.loads(_get(GRAPHQL, {"query": GQL % json.dumps(ids)}))
    best: dict[str, dict] = {}
    passing: list[dict] = []
    rejected = {}
    for e in body["data"]["entries"]:
        why = verdict(e)
        state = state_of(e["struct"]["title"])
        if why is None and state is None:
            why = "5: no state named in the title"
        if why:
            rejected[e["rcsb_id"]] = why
            continue
        res = e["rcsb_entry_info"]["resolution_combined"][0]
        row = {
                "pdb_id": e["rcsb_id"], "paralog": "RYR1",
                "organism": "Oryctolagus cuniculus", "uniprot": ACC,
                "resolution": res, "state": state, "title": e["struct"]["title"],
                "roles": ["ryr_panel"], "ip3_bound": False,
                "doi": (e["rcsb_primary_citation"] or {}).get("pdbx_database_id_DOI", "")}
        passing.append(row)
        if state not in best or res < best[state]["resolution"]:
            best[state] = row
    chosen = {r["pdb_id"]: r for r in best.values()}
    if "open" in best:
        partners = [r for r in passing if r["state"] == "primed"
                    and r["doi"] and r["doi"] == best["open"]["doi"]]
        if partners:
            start = min(partners, key=lambda r: r["resolution"])
            chosen.setdefault(start["pdb_id"], start)
            chosen[start["pdb_id"]]["roles"] = chosen[start["pdb_id"]]["roles"] + ["morph_start"]
            best["open"]["roles"] = best["open"]["roles"] + ["morph_end"]
    return sorted(chosen.values(), key=lambda r: r["pdb_id"]), rejected


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    seq, doms = sequence(), domains()
    structures, rejected = panel()
    for s in structures:
        print(f"{s['pdb_id']}  {s['state']:20s} {s['resolution']:.2f} A  "
              f"{','.join(s['roles'])}  {s['title']}")
    counts: dict[str, int] = {}
    for why in rejected.values():
        counts[why.split(":")[0]] = counts.get(why.split(":")[0], 0) + 1
    print(f"{len(rejected)} rejected, by first failed rule: {counts}")
    if args.dry_run:
        return 0
    out = {"provenance": {"retrieved": datetime.date.today().isoformat(),
                          "sources": _SOURCES,
                          "rules": {"min_chain_length": MIN_CHAIN_LENGTH,
                                    "min_modelled": MIN_MODELLED,
                                    "max_resolution_A": MAX_RESOLUTION,
                                    "allowed_ligands": sorted(ALLOWED_LIGANDS)}},
           "paralogs": {"RYR1": {"accession": ACC, "length": len(seq),
                                 "sequence": seq, "domains": doms}},
           "structures": structures, "rejected": rejected}
    path = RESOURCE_DIR / "ryr1.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
