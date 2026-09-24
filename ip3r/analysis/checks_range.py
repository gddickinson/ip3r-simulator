"""Paper 1 re-derived: where the receptor is, and where it has been lost.

Every count is rebuilt from per-proteome, per-record or per-genome rows with
this project's own code (:mod:`.range_table`) and compared with the number
the paper states. Presence is taken two ways: from the presence table's call
count, and record by record from the assignment tables. The genome-level
absences are rebuilt from the per-genome ledgers and compared with
``absence_at_genome.tsv`` row by row. The absence target list is rebuilt
from the one-sentence rule S23 states. Each outcome carries what its exhibit
draws.
"""

from __future__ import annotations

from collections import Counter

from ..core import genes_data as G
from ..parameters import PARAMETERS as _P
from . import range_table as RT
from .checks import agree, register

PROTEOME_SOURCES = (RT.PRESENCE, RT.TAXONOMY)
GENOME_SOURCES = (RT.MANIFEST, RT.CONTROLS, RT.COPIES, RT.COPY_LEDGER)

#: Clade counts the paper states: clade → (proteomes, with a call).
STATED_PRESENT = {"Arthropoda": (311, 292), "Nematoda": (110, 106),
                  "Mollusca": (32, 32), "Cnidaria": (14, 14), "Porifera": (3, 3)}
STATED_ABSENT = {"Streptophyta": (384, 0), "Chlorophyta": (48, 15),
                 "Ascomycota": (1034, 0), "Basidiomycota": (319, 0),
                 "Mucoromycota": (34, 18), "Chytridiomycota": (16, 6),
                 "Apicomplexa": (60, 0), "Microsporidia": (29, 0),
                 "Glomeromycota": (27, 0), "Mortierellomycota": (18, 0),
                 "Kickxellomycota": (35, 0), "Bacillariophyta": (16, 0),
                 "Rhodophyta": (12, 0)}


def _fmt(d: dict) -> str:
    return ", ".join(f"{c} {k}/{n}" for c, (n, k) in d.items())


def _clade_data(ps, clades) -> dict:
    return {"clades": [(c.clade, c.supergroup, c.n, c.present) for c in clades]}


@register("P1.presence_range", "range",
          "662 of the 6,928 non-vertebrate proteomes carry an IP3 receptor "
          "call, and 45 of the 135 clades swept do: 292 of 311 arthropod, 106 "
          "of 110 nematode and all molluscan, cnidarian and sponge proteomes; "
          "no archaeal (634) or bacterial (3,537) proteome.",
          "Proteomes counted per taxonomy clade from proteome_presence.tsv's "
          "call-count column (not its status label); presence re-derived "
          "record by record from the six assignment tables and compared "
          "taxon by taxon.",
          "rederived", PROTEOME_SOURCES + RT.ASSIGNMENTS)
def presence_range():
    ps = RT.load_proteomes()
    clades = RT.clade_rows(ps)
    present = sum(p.present for p in ps)
    with_call = sum(c.present > 0 for c in clades)
    by_clade = {c.clade: (c.n, c.present) for c in clades}
    named = {c: by_clade.get(c, (0, 0)) for c in STATED_PRESENT}
    prok = {d: (sum(p.domain == d for p in ps), sum(p.present for p in ps if p.domain == d))
            for d in ("Archaea", "Bacteria")}
    calls = RT.taxon_calls()
    per_taxon: dict[str, int] = Counter()
    for p in ps:
        per_taxon[p.taxid] += p.n_itpr
    disagree = sorted(t for t in set(per_taxon) | set(calls)
                      if (calls.get(t, 0) > 0) != (per_taxon.get(t, 0) > 0))
    ok = ((present, len(ps), with_call, len(clades)) == (662, 6928, 45, 135)
          and named == STATED_PRESENT
          and prok == {"Archaea": (634, 0), "Bacteria": (3537, 0)} and not disagree)
    return agree(ok, f"662/6,928 proteomes; 45/135 clades; {_fmt(STATED_PRESENT)}; "
                 "Archaea 0/634, Bacteria 0/3,537",
                 f"{present}/{len(ps):,} proteomes; {with_call}/{len(clades)} clades; "
                 f"{_fmt(named)}; Archaea {prok['Archaea'][1]}/{prok['Archaea'][0]}, "
                 f"Bacteria {prok['Bacteria'][1]}/{prok['Bacteria'][0]:,}",
                 f"presence from the assignment records agrees with the presence "
                 f"table in {len(per_taxon) - len(disagree)} of {len(per_taxon)} taxa",
                 **_clade_data(ps, clades))


@register("P1.kingdom_absences", "range",
          "Land plants and Dikarya have lost the receptor their relatives "
          "retain: none of 384 Streptophyta proteomes against 15 of 48 "
          "Chlorophyta; none of 1,034 Ascomycota or 319 Basidiomycota against "
          "18 of 34 Mucoromycota and 6 of 16 Chytridiomycota; also none of 60 "
          "apicomplexan, 29 microsporidian, 27 glomeromycote, 18 "
          "mortierellomycote, 35 kickxellomycote, 16 diatom and 12 red algal "
          "proteomes.",
          "Each named phylum counted from the per-proteome presence table "
          "joined to the taxonomy.",
          "rederived", PROTEOME_SOURCES)
def kingdom_absences():
    ps = RT.load_proteomes()
    got = {c: RT.count_at(ps, "phylum", c) for c in STATED_ABSENT}
    ok = got == STATED_ABSENT
    return agree(ok, _fmt(STATED_ABSENT), _fmt(got),
                 **_clade_data(ps, RT.clade_rows(ps)))


@register("P1.relaxed_controls", "range",
          "At the relaxed sensitivity PF08709 returns no substantial match in "
          "land plants against 26 in Chlorophyta, and none in Dikarya against "
          "16 in Mucoromycota, while the MIR domain PF02815 returns 633 in land "
          "plants and 4,376 in Dikarya; neither full-length profile has a "
          "substantial match in either absence, and the one non-MIR exception "
          "is a short accessory domain, not a receptor.",
          "relaxed_hits.tsv filtered by this project at the registered "
          "substantial-match bar (range.substantial_evalue, "
          "range.substantial_coverage), each target counted once per "
          "lineage, lineages placed by the taxonomy table.",
          "rederived", (RT.RELAXED, RT.TAXONOMY))
def relaxed_controls():
    res = RT.substantial_table()
    t = res["table"]

    def sub(lin, prof):
        return t.get(lin, {}).get(prof, (0, 0))[1]

    want = {("land plants", "PF08709"): 0, ("Chlorophyta", "PF08709"): 26,
            ("Dikarya", "PF08709"): 0, ("Mucoromycota", "PF08709"): 16,
            ("land plants", "PF02815"): 633, ("Dikarya", "PF02815"): 4376}
    got = {k: sub(*k) for k in want}
    full = {lin: (sub(lin, "itpr"), sub(lin, "ryr")) for lin in ("land plants", "Dikarya")}
    exc = res["exceptions"]
    receptor_like = [a for a, ps in exc.items()
                     if {"PF08709", "itpr", "ryr"} & set(ps)]
    ok = got == want and all(v == (0, 0) for v in full.values()) \
        and not receptor_like and len(exc) == 1
    fmt = "; ".join(f"{p} {lin} {v:,}" for (lin, p), v in got.items())
    return agree(ok, "; ".join(f"{p} {lin} {v:,}" for (lin, p), v in want.items())
                 + "; full-length profiles 0; 1 non-MIR exception",
                 f"{fmt}; full-length (itpr, ryr) {full}; non-MIR exceptions "
                 f"{len(exc)} ({', '.join(f'{a}: {p}' for a, p in exc.items())})",
                 f"receptor-like exceptions: {len(receptor_like)}",
                 table=t, lineages=list(RT.LINEAGES))


@register("P1.absence_targets", "range",
          "The genome sweep's absence targets are re-derived from the "
          "presence table: every clade with at least ten proteomes swept and "
          "no call, 35 of them, including three no summary had named: diatoms "
          "(0 of 16), red algae (0 of 12) and Cestoda (0 of 11).",
          "S23's rule G3 re-implemented from its statement (a eukaryotic "
          "phylum or class, range.absence_min_proteomes swept, none called) "
          "on the per-proteome table, and compared with the clades "
          "absence_at_genome.tsv reports, including each one's proteome count.",
          "rederived", PROTEOME_SOURCES + (RT.ABSENCE,))
def absence_targets():
    ps = RT.load_proteomes()
    ours = RT.absence_targets(ps)
    table = {(r["rank"], r["clade"]): int(r["swept_proteomes"]) for r in G.read_tsv(RT.ABSENCE)}
    swept = {k: RT.count_at(ps, *k)[0] for k in ours}
    three = {c: RT.count_at(ps, r, c) for r, c in
             (("phylum", "Bacillariophyta"), ("phylum", "Rhodophyta"), ("class", "Cestoda"))}
    ok = (set(ours) == set(table) and len(ours) == 35
          and all(swept[k] == table[k] for k in ours)
          and three == {"Bacillariophyta": (16, 0), "Rhodophyta": (12, 0), "Cestoda": (11, 0)})
    only_ours = sorted(set(ours) - set(table))
    only_table = sorted(set(table) - set(ours))
    return agree(ok, "35 targets; Bacillariophyta 0/16, Rhodophyta 0/12, Cestoda 0/11",
                 f"{len(ours)} targets by the rule at ≥ "
                 f"{_P.value('range.absence_min_proteomes'):.0f} proteomes; {_fmt(three)}",
                 f"only by the rule: {only_ours or 'none'}; only in the table: "
                 f"{only_table or 'none'}", **_clade_data(ps, RT.clade_rows(ps)))


@register("P1.absences", "range",
          "35 clade-level absences are confirmed in genome assemblies, each "
          "genome carrying a positive control chosen for its clade: Ascomycota "
          "has no receptor in 31 controlled assemblies, Streptophyta none in 25 "
          "and Basidiomycota none in 17; two land-plant assemblies carry only "
          "translated-search traces; and no absence rests on the two ciliate "
          "genomes whose control is partial or absent.",
          "Each target's genomes found in the S23 manifest by rank and name; "
          "controlled, complete-gene and trace-only genomes counted from the "
          "control ledger, copies.tsv and the copy-number ledger; each "
          "rebuilt row compared with absence_at_genome.tsv.",
          "rederived", PROTEOME_SOURCES + GENOME_SOURCES + (RT.ABSENCE,))
def absences():
    ps = RT.load_proteomes()
    table = {(r["rank"], r["clade"]): r for r in G.read_tsv(RT.ABSENCE)}
    rows = RT.genome_absences(ps, sorted(table))
    cols = ("swept_proteomes", "genomes_swept", "genomes_controlled",
            "genomes_with_full_itpr", "genomes_with_trace_only")
    differ = [r.clade for r in rows if tuple(int(table[(r.rank, r.clade)][c]) for c in cols)
              != (r.proteomes, r.genomes, r.controlled, r.with_full, r.trace_only)]
    held = [r for r in rows if r.holds]
    big = {r.clade: (r.controlled, r.with_full) for r in rows
           if r.clade in ("Ascomycota", "Streptophyta", "Basidiomycota")}
    traces = sum(r.trace_only for r in rows if r.clade == "Streptophyta")
    control = {r["accession"]: r["verdict"] for r in G.read_tsv(RT.CONTROLS)}
    weak = {a for a, v in control.items() if v not in RT.CONTROLLED}
    target_genomes = {m["accession"] for m in G.read_tsv(RT.MANIFEST)
                      if (("phylum", m["phylum"]) in table or ("class", m["class"]) in table)}
    resting = weak & target_genomes
    ok = (len(held) == len(rows) == 35 and not differ and traces == 2 and len(weak) == 2
          and not resting and big == {"Ascomycota": (31, 0), "Streptophyta": (25, 0),
                                      "Basidiomycota": (17, 0)})
    return agree(ok, "35 hold; Ascomycota 31, Streptophyta 25, Basidiomycota 17 controlled, "
                 "none with a gene; 2 trace-only land plants; 2 weakly controlled genomes, "
                 "no absence resting on either",
                 f"{len(held)}/{len(rows)} hold; " + ", ".join(
                     f"{c} {n} controlled, {k} with a gene" for c, (n, k) in big.items())
                 + f"; {traces} trace-only land plants; {len(weak)} weakly controlled, "
                 f"{len(resting)} in an absence target",
                 f"rebuilt rows differ from absence_at_genome.tsv in "
                 f"{len(differ)} of {len(rows)} ({', '.join(differ) or 'none'})",
                 absences=[(r.clade, r.genomes, r.controlled, r.with_full) for r in rows])


@register("P1.copy_number", "range",
          "Of the 194 non-vertebrate genomes, 117 carry no IP3 receptor gene "
          "and 43 carry one; the flatworm Macrostomum lignano carries 18 "
          "complete gene models, the ciliate Stentor coeruleus 13, the sponge "
          "Dysidea avara 8, and the green alga Cymbomonas tetramitiformis 3.",
          "Complete gene models counted per genome from the rows of copies.tsv, "
          "zeros filled from the manifest, and the named species found by "
          "their counts rather than their counts by name.",
          "rederived", (RT.COPIES, RT.MANIFEST))
def copy_number():
    cn = RT.copy_numbers()
    dist = Counter(n for _, n in cn.values())
    top = sorted(cn.values(), key=lambda x: -x[1])[:3]
    cym = [n for org, n in cn.values() if org.startswith("Cymbomonas")]
    want_top = [("Macrostomum lignano", 18), ("Stentor coeruleus", 13), ("Dysidea avara", 8)]
    ok = (len(cn), dist[0], dist[1]) == (194, 117, 43) and cym == [3] and \
        [(o.split(" (")[0], n) for o, n in top] == want_top
    return agree(ok, "194 genomes; 117 none, 43 one; Macrostomum lignano 18, Stentor "
                 "coeruleus 13, Dysidea avara 8; Cymbomonas 3",
                 f"{len(cn)} genomes; {dist[0]} none, {dist[1]} one; "
                 + ", ".join(f"{o} {n}" for o, n in top) + f"; Cymbomonas {cym}",
                 copies=sorted(n for _, n in cn.values()))


@register("P1.record_chase", "range",
          "All 99 plant and fungal IP3 receptor records (64 plant, 35 fungal) "
          "are real genes (47) or fragments of one (52); none is a "
          "contamination suspect and none lacks genome backing. Surviving "
          "plant records sit at 19.9–39.9 % identity to anything outside "
          "their kingdom and fungal ones at 20.5–33.5 %; no record of the 99 "
          "exceeds 45.8 %.",
          "Each record's verdict re-applied from its own columns with the "
          "registered chase thresholds: a fragment below "
          "range.family_floor_aa or flagged a fragment by UniProt, a contaminant at range.contaminant_pident "
          "over range.contaminant_qcov to one relative outside its kingdom, "
          "and genome backing from its EMBL and proteome cross-references.",
          "rederived", (RT.CHASE,))
def record_chase():
    rows = G.read_tsv(RT.CHASE)
    floor = _P.value("range.family_floor_aa")
    pid, qcov = _P.value("range.contaminant_pident"), _P.value("range.contaminant_qcov")
    ident = [float(r["out_kingdom_pident"]) for r in rows]
    contam = [r["accession"] for r, i in zip(rows, ident)
              if i >= pid and float(r["out_kingdom_qcov"]) >= qcov]
    unbacked = [r["accession"] for r in rows
                if int(r["n_embl"]) == 0 and int(r["n_proteomes"]) == 0]
    # R5: short of the floor, *or* flagged a fragment by UniProt at any length
    ours = ["fragment" if int(r["length"]) < floor or r["fragment"] == "fragment"
            else "real_gene" for r in rows]
    label_diff = sum(o != r["verdict"] for o, r in zip(ours, rows))
    kingdoms = Counter(r["kingdom"] for r in rows)
    real = Counter(ours)
    spans = {}
    for k in ("Viridiplantae", "Fungi"):
        v = [i for i, o, r in zip(ident, ours, rows) if o == "real_gene" and r["kingdom"] == k]
        spans[k] = (min(v), max(v))
    ok = ((real["real_gene"], real["fragment"]) == (47, 52) and not contam and not unbacked
          and (kingdoms["Viridiplantae"], kingdoms["Fungi"]) == (64, 35)
          and spans == {"Viridiplantae": (19.9, 39.9), "Fungi": (20.5, 33.5)}
          and max(ident) == 45.8)
    return agree(ok, "64 plant + 35 fungal; 47 real, 52 fragments; 0 contaminants; "
                 "plants 19.9–39.9 %, fungi 20.5–33.5 %; max 45.8 %",
                 f"{kingdoms['Viridiplantae']} plant + {kingdoms['Fungi']} fungal; "
                 f"{real['real_gene']} real, {real['fragment']} fragments; "
                 f"{len(contam)} contaminants, {len(unbacked)} unbacked; "
                 + ", ".join(f"{k} {a}–{b} %" for k, (a, b) in spans.items())
                 + f"; max {max(ident)} %",
                 f"the re-applied verdict differs from the table's in {label_diff} of "
                 f"{len(rows)} records",
                 points=[(int(r["length"]), i, o, r["kingdom"])
                         for r, i, o in zip(rows, ident, ours)],
                 floor=floor, contaminant=pid)
