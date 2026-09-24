"""Structural findings, recomputed from the deposited coordinates.

S0 measured the channel on 6DQN with a centroid-plane axis; these checks
find the axis by superposition instead (``structure.symmetry``) and re-measure
everything downstream with this project's code. S22's ligand shells are
recomputed on all six IP3-bound depositions.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from ..core import genes_data as G
from ..io.loader import load
from ..parameters import PARAMETERS as _P
from ..structure.channel import measure_channel
from ..structure.ligand import ligand_sites, residue_distances
from ..structure.pore import pore_profile
from .checks import agree, not_run, register

META = "s0_baseline/review_figures/structure_meta.json"
PORE = "s0_baseline/review_figures/structure_pore.tsv"
SHELLS = "ligand_site/shell_agreement.tsv"


@lru_cache(maxsize=4)
def _summary(pdb: str):
    return measure_channel(load(pdb), include_hetero=True)


_P.subscribe(_summary.cache_clear)      # measured with registered constants

def _tol() -> float:
    return _P.value("check.length_tol")


@register("S0.c4_symmetry", "structure",
          "6DQN is four-fold symmetric: subunit A turned by 90° lands on its "
          "neighbours to 0.058 Å RMSD.",
          "Axis from Kabsch superposition of each subunit on the next (not "
          "the centroid plane S0 used); the residual re-measured about it.",
          "recomputed", (META,), ("6DQN",))
def c4_symmetry():
    meta = G.read_json(META)
    s = _summary("6DQN")
    axis0 = np.asarray(meta["four_fold_axis_in_deposit"], float)
    tilt = float(np.degrees(np.arccos(min(1.0, abs(s.frame.axis @ axis0 /
                                                   np.linalg.norm(axis0))))))
    ok = (abs(s.c4_residual - meta["c4_residual_rmsd_A"]) <= _tol()
          and abs(s.superposition_angle - 90.0) < 0.5 and tilt < 1.0)
    return agree(ok, f"{meta['c4_residual_rmsd_A']} Å",
                 f"{s.c4_residual:.3f} Å; subunit-to-subunit rotation "
                 f"{s.superposition_angle:.2f}°; axis within {tilt:.3f}° of S0's",
                 "Two different axis constructions agree; the residual is "
                 "compared within the registered geometry tolerance.",
                 residual=s.c4_residual, angle=s.superposition_angle, tilt=tilt)


def _constriction_check(name: str, meta_key: str, lining_key: str):
    meta = G.read_json(META)
    s = _summary("6DQN")
    c = s.constrictions.get(name)
    if c is None:
        return not_run(f"no {name} constriction found")
    pub_z, pub_r = meta[f"{meta_key}_z_A"], meta[f"{meta_key}_min_radius_A"]
    pub_res = meta[lining_key]
    ok = (abs(c.radius - pub_r) <= _tol() and abs(c.z - pub_z) <= 1.0
          and c.residues == pub_res)
    return agree(ok, f"r = {pub_r} Å at z = {pub_z} Å, lined by {', '.join(pub_res)}",
                 f"r = {c.radius:.2f} Å at z = {c.z:.1f} Å, lined by "
                 f"{', '.join(c.residues)}", z=c.z, radius=c.radius,
                 residues=c.residues)


@register("S0.selectivity_filter", "structure",
          "The narrowest luminal point of the pore is the GGGVGD selectivity "
          "filter, lined by Asn2472 and Gly2473.",
          "Pore profile re-measured about the independently found axis; the "
          "luminal minimum and the residues lining it re-identified.",
          "recomputed", (META,), ("6DQN",))
def selectivity_filter():
    out = _constriction_check("filter", "filter", "filter_lining_residues")
    meta = G.read_json(META)
    motif = meta["filter_motif_start"]
    nums = [int("".join(ch for ch in r if ch.isdigit()))
            for r in out.data.get("residues", [])]
    on_motif = any(motif - 2 <= n <= motif + 6 for n in nums)
    if not on_motif and out.status == "confirmed":
        out.status = "discrepancy"
    out.detail = (f"Lining residues {'fall' if on_motif else 'do NOT fall'} on "
                  f"the GGGVGD motif at {motif}-{motif + 5}. The motif is not "
                  "used to find the constriction, so landing on it is a test.")
    return out


@register("S0.gate", "structure",
          "The narrowest cytosolic point is the gate, 2.52 Å, lined by "
          "Phe2513 and Ile2517.",
          "As for the filter, on the cytosolic half of the pore domain.",
          "recomputed", (META,), ("6DQN",))
def gate():
    return _constriction_check("gate", "gate", "gate_lining_residues")


@register("S0.pore_profile", "structure",
          "The committed 6DQN pore-radius profile (minimum heavy-atom "
          "distance to the axis).",
          "Profile re-measured on S0's own z grid about the independently "
          "fitted axis; every point compared.",
          "recomputed", (PORE,), ("6DQN",))
def pore_profile_check():
    rows = G.read_tsv(PORE)
    z = np.array([float(r["z_along_axis"]) for r in rows])
    r_pub = np.array([float(r["min_heavy_atom_radius"]) for r in rows])
    s = _summary("6DQN")
    step = float(np.median(np.diff(z)))
    prof = pore_profile(load("6DQN"), s.frame, z[0], z[-1], step=step,
                        include_hetero=True)
    r_new = np.interp(z, prof.z, prof.r_min)
    diff = np.abs(r_new - r_pub)
    tol = _tol()
    frac = float((diff <= tol).mean())
    ok = frac >= 0.95
    return agree(ok, f"{len(z)} points from z = {z[0]:.1f} to {z[-1]:.1f} Å",
                 f"{frac:.1%} of points within {tol} Å (median |Δ| "
                 f"{np.median(diff):.3f} Å, max {diff.max():.3f} Å)",
                 "Rule: at least 95 % of the grid within the geometry tolerance "
                 "— a slightly different axis moves the atoms nearest the axis "
                 "by hundredths of an Å, and a real disagreement moves many "
                 "points.", z=z.tolist(), r_published=r_pub.tolist(),
                 r_recomputed=r_new.tolist())


@register("S0.ip3_contacts", "structure",
          "IP3 is contacted (≤ 4.5 Å) by ten residues of its own subunit and "
          "none of a neighbour's: R266, T268, L269, R270, R503, K507, R510, "
          "Y567, R568, K569.",
          "Contacts re-measured at every one of the four sites, each site "
          "assigned to the subunit its atoms are nearest.",
          "recomputed", (META,), ("6DQN",))
def ip3_contacts():
    meta = G.read_json(META)
    pub = meta["ip3_contacts"]["same_subunit"]
    s = _summary("6DQN")
    per_site = {k: [f"{a}{b}" for a, b in v["same_subunit"]]
                for k, v in s.ip3_contacts.items()}
    cross = {k: v["other_subunit"] for k, v in s.ip3_contacts.items()
             if v["other_subunit"]}
    ok = (len(per_site) == 4 and all(v == pub for v in per_site.values())
          and not cross)
    found = (f"{len(per_site)} sites; " + "; ".join(
        f"{k}: {'identical' if v == pub else ', '.join(v)}" for k, v in per_site.items()))
    return agree(ok, ", ".join(pub), found,
                 "Cross-subunit contacts: " + (str(cross) if cross else "none"),
                 contacts=per_site)


def _shell_contacts(pdb: str, heavy_only: bool) -> tuple[list[int], dict]:
    """Residues within the cutoff, best over subunits (S22's rule)."""
    st = load(pdb)
    cutoff = _P.value("ligand.contact_cutoff")
    best: dict[int, float] = {}
    for site in ligand_sites(st):
        for resi, d in residue_distances(st, site, heavy_only=heavy_only).items():
            best[resi] = min(best.get(resi, np.inf), d)
    return sorted(r for r, d in best.items() if d <= cutoff), best


def _s0_numbers() -> list[int]:
    return [int("".join(c for c in r if c.isdigit()))
            for r in G.read_json(META)["ip3_contacts"]["same_subunit"]]


@register("P6.shell_agreement", "ligand",
          "S0's ten IP3 contacts are recovered in all six IP3-bound "
          "depositions (6DQN, 8TKG, 8TKF, 8TKH, 7T3P, 8TLA), with the "
          "contact counts and extra contacts in shell_agreement.tsv.",
          "Minimum distances re-measured per residue over all atoms, "
          "hydrogens included, best over the subunits whose own site is "
          "occupied — S22's rule exactly.",
          "recomputed", (SHELLS, META), ("6DQN", "8TKG", "8TKF", "8TKH", "7T3P", "8TLA"))
def shell_agreement():
    s0 = _s0_numbers()
    lines, ok = [], True
    for row in G.read_tsv(SHELLS):
        got, _ = _shell_contacts(row["pdb_id"], heavy_only=False)
        extra = [r for r in got if r not in s0]
        pub_extra = [int(x) for x in row["extra_contacts"].split(";") if x]
        n_s0 = sum(r in got for r in s0)
        same = (len(got) == int(row["n_contact_le_4.5A"])
                and n_s0 == int(row["n_s0_contacts_recovered"]) and extra == pub_extra)
        ok &= same
        lines.append(f"{row['pdb_id']}: {len(got)} contacts, {n_s0}/10 S0, "
                     f"extra {extra or '—'}" + ("" if same else " ≠ published"))
    return agree(ok, "six depositions, all 10/10", "; ".join(lines))


@register("P6.contacts_heavy_atom", "ligand",
          "The positive control holds under S0's own contact definition: "
          "the ten contacts are within 4.5 Å by heavy atoms in all six "
          "IP3-bound depositions.",
          "The S22 recomputation repeated with hydrogens excluded, the "
          "convention S0 defined the ten contacts with.",
          "recomputed", (SHELLS, META), ("6DQN", "8TKG", "8TKF", "8TKH", "7T3P", "8TLA"))
def contacts_heavy_atom():
    s0 = _s0_numbers()
    lines, ok = [], True
    for row in G.read_tsv(SHELLS):
        got, best = _shell_contacts(row["pdb_id"], heavy_only=True)
        miss = {r: round(best.get(r, np.inf), 2) for r in s0 if r not in got}
        ok &= not miss
        lines.append(f"{row['pdb_id']}: {10 - len(miss)}/10"
                     + (f" (missing {', '.join(f'{r} at {d} Å' for r, d in miss.items())})"
                        if miss else ""))
    return agree(ok, "10/10 in all six", "; ".join(lines),
                 "S22's all-atom rule counts hydrogens; this is the heavy-atom "
                 "reading of the same structures. A miss here means the "
                 "positive control passes only through a hydrogen contact.")
