# SCIENCE_CHECKS — the findings checks, paper by paper

Moved out of `SCIENCE.md` in Round 7.7 to keep both files short. Every check
here is calibrated by a planted flip (`tests/test_checks_calibration.py`).

## Paper 6: the module contrast

The **ligand core** is the smallest span that holds all ten IP3 contacts. The
**pore module** is PF00520 less the luminal loop. Both are rebuilt in
`core.modules` from the imported sites and domain map, and validated against
what they must contain. For each orthologue in S17's deep alignment, the
identity to the human reference is computed in each module. Only columns the
tip covers are counted, and a tip must cover at least 50 % of each module
(`ligand.module_min_coverage`). The test is on the paired difference core −
pore: an exact sign test, and a signed-rank test that drops zeros, corrects
for ties and applies no continuity correction. The residue → column map comes
from walking the alignment's reference row, and is refused unless that row's
ungapped sequence is the UniProt sequence. S17's `deep_col` column is not
used.

What agreement establishes: from the alignments, the published spans, tip
counts, means, sign counts and p-values follow by independent code, and the
reversal with the loop counted in is real. It says nothing about the
alignments themselves.

## Paper 6: ligand shells

For each residue, the distance is the minimum all-atom distance to the IP3 on
its own subunit, with hydrogens included as in S22's reader. Per deposit it
is the best over the four subunits, and the consensus is the median over the
six IP3-bound ITPR3 depositions. Shells are half-open: [0, 4.5), [4.5, 8),
[8, 11.5) and [11.5, 15) Å. Residues beyond 15 Å belong to no shell.

S22 carried the pocket to ITPR1/ITPR2 with a MAFFT pairwise alignment. Here
it is carried with a separate alignment: Gotoh's affine-gap global
alignment, BLOSUM62, gap open 10 and extend 0.5 (EMBOSS needle's defaults),
end gaps free (`core/pairwise.py`). All 250 transferred positions agree with
S22's. Conservation is S17's per-residue deep JSD. Each shell is tested
against every scored residue with a one-sided Mann–Whitney (normal
approximation, tie and continuity corrections). The trend is Spearman's ρ
with a t-distribution p.

"No step at 4.5 Å" is prose, not a table, so it is tested two ways:

1. The contact shell does not beat the second (one-sided Mann–Whitney).
   With 12 and 14 residues this test is weak on its own.
2. The contact boundary does not carry the largest drop in mean JSD
   between adjacent shells.

It holds both ways in all three paralogs. The largest drop is at the 11.5 Å
boundary (third → fourth) everywhere. Note that the FEL purifying fraction
(the paper's §8) falls after the *second* shell instead, so the two
instruments put the drop in different places.

## Paper 6 §8: the same questions asked of the substitution rate

S17's FEL run (`constraint/fel_sites.tsv`) gives each codon site α, β, q and
a verdict. `analysis/fel_rates.py` joins it to structure built here: the
pocket recomputed from six ITPR3 deposits and carried by this project's own
alignment, and the modules rebuilt from the annotation. As in S22, β is
compared by rank, because it is exactly zero at most sites. A site is
purifying at FEL's own q ≤ `check.fel_q`.

- **The literature core** (`ibc_literature`) is Bosanac 2002's crystallised
  IP3-binding core, ITPR1 224–604 (`ligand.ibc_start`/`ligand.ibc_end`).
  It is carried to ITPR2 and ITPR3 by the Gotoh alignment, which lands on
  S17's transfer exactly (224–604, 225–604). With it, `P6.module_map` now
  compares all twelve spans.
- **The q-values are corrected together**, as S22 corrects them: the twelve
  shell rows in one Benjamini–Hochberg family, and the nine module rows
  (three definitions × three paralogues) in another. The primary row's q
  can only be reproduced with all nine rows present.
- All 21 rows reproduce: counts, means, the common-language effect, p and
  q. The shell share falls from 1.000 at contact to a floor, and the fourth
  shell sits above the third in all three paralogues. On β, the pore
  evolves faster than the core in ITPR1 (q 6.2 × 10⁻⁸), against the
  conservation result.


The P2 tree checks ask their questions of the committed `rooted.nwk` with this
project's Newick reader (`analysis/newick.py`, `analysis/tree.py`). No
ip3r_genes clade table goes into the answer. A tip's group comes from the
census prefix on its label. Cyclostomes are recognised by genus (*Myxine*,
*Petromyzon* and the other hagfish and lamprey genera).

- **A paralog's whole clade** is the MRCA of every tip whose record names it.
  The check requires that clade to contain nothing foreign, meaning no other
  paralog, no cyclostome and no invertebrate. A second rule, "expand the
  named core while the clade stays pure", gives the same 19/13/19.
- **A cyclostome clade** is a maximal clade made only of cyclostome tips. A
  clade "branches first among the vertebrates" when its parent holds exactly
  the vertebrate tips.
- **Support is counted per bipartition.** In a rooted tree the root's two
  edges are a single bipartition of the unrooted tree, and IQ-TREE labels
  both of them. Counting nodes gave 92 of 132; counting bipartitions gives
  the published 91 of 131. The bars are the registered parameters
  `tree.alrt_min` (80, Guindon 2010) and `tree.ufboot_min` (95, Hoang 2018).

### The `--bnni` re-search (S7 §5.7)

`analysis/tree_robustness.py` re-asks every clade claim of UFBoot's
model-violation guard. The claims' tip sets are rebuilt from the census
prefixes: a paralogue's **core** is its largest *pure* clade. That rule
matters: five shark and coelacanth records name a paralogue but sit outside
its pure clade, so "every tip naming it" gives 15/11/18 and not S7's
13/10/16. The pairs and the triple are unions of the MRCA clades. Both
trees are rooted on RyR with `newick.reroot`, which carries each support
label with its edge. Each set is then read as a clade in each tree. All 13
sets equal `claim_members.tsv`. Nine of the ten clades hold. The ITPR1
core (47.8/95 → 47.5/73) was never well supported. The Tree tab's
"Beside --bnni" draws the two trees side by side.

## Papers 3 and 4: the genome × paralog grid

The retention sweep asked 309 vertebrate assemblies for four cells each:
ITPR1, ITPR2, ITPR3 and a ryanodine-receptor control. `analysis/genome_grid.py`
joins three per-cell tables into one grid: what the search found
(`contiguity_cells.tsv`), the S15a evidence state (`character_matrix.tsv`)
and protein-record recovery (`gene_recovery.tsv`).

- **The contiguity bar** is the registered `genomes.contiguity_bar_bp`
  (142,212 bp, D4's median measured ITPR span). This project splits genomes
  on contig N50 ≥ bar with its own comparison, and the check tests that the
  split agrees with the table's `spans_gene` column in every cell.
- **The method control** is recomputed with `analysis/stats.py`: Wilson
  intervals, a two-sided Fisher exact test (hypergeometric sum, R's 1e-7
  tie tolerance), and a logistic regression of found on log10 N50 (IRLS,
  Wald p). The Mann–Whitney test of found against missed N50 is
  **one-sided**. The first version doubled it and was off by exactly a
  factor of 2, which was the checker's error.
- **The recovery channel** is rebuilt from the table's count columns in
  this order: gene present? reference proteome? any family record? any
  full-length record resolving to a cell? one resolving to *this* cell?
  It agrees with the table's own label in all 1,236 cells.

## Paper 1: presence and absence across eukaryotes

`analysis/range_table.py` joins the S20 proteome sweep (per-proteome
presence, the taxonomy, the per-record assignments, the relaxed search) and
the S23 genome sweep (manifest, control ledger, copies, copy-number ledger).
A clade is the taxonomy table's own `clade`. Seven checks re-derive the
paper's numbers from those rows:

- **Presence two ways.** Presence is counted from the presence table's call
  count, and again record by record from the six assignment tables. The two
  agree in all 6,854 taxa. They can only be compared per taxon, because 72
  taxa have more than one proteome.
- **Absence targets from the rule.** S23's rule G3 is re-implemented from its
  one-sentence statement: a eukaryotic phylum or class with ≥
  `range.absence_min_proteomes` (10) swept and no call. It gives exactly the
  35 clades `absence_at_genome.tsv` reports.
- **Genome absences from the ledgers.** Each target's genomes, controlled
  genomes (`controlled_cross_kingdom` or `controlled_by_target` only),
  complete-gene genomes and trace-only genomes are counted from the
  per-genome ledgers. All 35 rows match the published table.
- **Substantial matches** are full E ≤ `range.substantial_evalue` (1e-5)
  and model coverage ≥ `range.substantial_coverage` (0.5), with each target
  counted once per lineage. All 24 cells of S20's four-lineage table
  reproduce.
- **The record chase.** Rule R5 calls a record a fragment if it is below
  `range.family_floor_aa` (2,000) *or* UniProt flags it as a fragment. The
  first checker applied only the floor and got 50 real genes instead of 47:
  three records of 2,366–2,858 aa are UniProt-flagged fragments. That was
  the checker's error.

### The family-call benchmark (S1)

`analysis/family_benchmark.py` measures the labelled-bait margin again from
the 56 UniProt sequences S1 committed. Each control is aligned pairwise to
each human ITPR and RyR bait by this project's aligner, not S1's single
MAFFT alignment. The margin is the identity to the nearest ITPR bait minus
the identity to the nearest RyR bait.

- **Full-alignment identity** (identical columns over columns where either
  sequence has a residue) is S1's scorer metric. It calls all 25 positives
  ITPR and all 6 RyRs RyR. The gap is 0.582 (S1: 0.601), and iplA is again
  the narrowest true member, at +0.048 (+0.065), inside the 0.10 no-call
  band. The largest margin difference is 0.029.
- **Covered-column identity is not re-measured by this route.** A pairwise
  optimum aligns unrelated long sequences at a covered identity of
  0.25–0.28. So iplA sits at that floor against both families (−0.002
  against S1's +0.042), and at gap costs 11/1 the iplA–RyR alignment
  collapses to about 54 pairs. An MSA constrains those columns; a pairwise
  alignment does not. That is a limit of this route, not an error in S1.
- **The one miss is within aligner noise.** *Drosophila* Itpr against worm
  itr-1 reads 0.361 pairwise against 0.342 in S1's MSA. The breadth bar is
  0.35, so the paper's "0.008 short" is true of S1's instrument, but it is
  smaller than the difference between two alignment routes.
- **The counts** are rebuilt from each control's fired components with the
  registered `bench.points_*`. The sister cap is decided by the recomputed
  margin, and the evidence gate by the family-specific components. Every
  score and cap matches: recall 24/25, specificity 31/31, RyR 6/6.

## Paper 5 §8: the VUS stratification

`analysis/vus_strata.py` rebuilds S17's rule from its description. For one
gene on one layer, each class (P/LP, B/LB, VUS) is scored once per
**position**. A residue scores only at occupancy ≥
`constraint.min_occupancy` (0.5). On the deep layer that is the same set as
the `deep_reliable` flag: the two agree at every residue of all three
paralogs. A residue carrying alleles of two classes counts in both. The
thresholds are the labelled sets' own medians. A VUS at or above the P/LP
median is pathogenic-like, and one at or below the B/LB median is
benign-like; ties count. `P5.vus_stratification` rebuilds all 12 rows of
`vus_stratification.tsv` from `variants.tsv` and the per-residue tables, and
every field agrees.

The stratification takes every source, so the 13 curated UniProt P/LP
records are in its P/LP median. The classifier test (`P5.variant_auc`) is
ClinVar only. Run ClinVar-only, the table is unchanged: every curated record
sits at a position ClinVar already labels P/LP, or at one no layer scores. So
the difference in rules changes no number here.

The viewer computes the same stratification from the committed resources
(`constraint.json`, `variants.json`), and a test proves that it reproduces
S17's table.

## Every check, as of Round 7.7

| check | kind | verdict | re-derived |
|---|---|---|---|
| `S0.c4_symmetry` | recomputed | confirmed | 0.051 Å (S0: 0.058); rotation 90.00°; axes agree to 0.02° |
| `S0.selectivity_filter` | recomputed | confirmed | 5.08 Å (5.06) at −89.8 Å, Asn2472/Gly2473, on GGGVGD |
| `S0.gate` | recomputed | confirmed | 2.55 Å (2.52), Phe2513/Ile2517 |
| `S0.pore_profile` | recomputed | confirmed | 99.3 % of 144 points within 0.05 Å |
| `S0.ip3_contacts` | recomputed | confirmed | the same ten residues at all four sites |
| `P6.shell_agreement` | recomputed | confirmed | all six depositions, counts and extras identical |
| `P6.contacts_heavy_atom` | recomputed | **discrepancy** | R503 outside 4.5 Å by heavy atoms in 8TKG, 8TKH |
| `P5.element_means` | rederived | confirmed | 27 element means, largest Δ 0.0000 |
| `P5.gate_filter_most_conserved` | rederived | confirmed | as the paper states, both metrics |
| `P5.report_both_metrics` | rederived | **discrepancy** | report sentence overstates the JSD ranking |
| `P5.luminal_loop_least` | rederived | confirmed | last in all three |
| `P5.gate_identical` | rederived | confirmed | FGVII in all three |
| `P5.variant_auc` | rederived | confirmed | all AUCs and position counts |
| `P5.deep_ranks_third` | rederived | confirmed | family 0.872 > vert 0.854 > deep 0.758 > shallow 0.684 |
| `P5.vus_count` | rederived | confirmed | 1,546 |
| `P5.vus_stratification` | rederived | confirmed | all 12 gene × layer rows, every count, median and fraction; ClinVar-only gives the same table |
| `P5.omega_range` | read | confirmed | ω 0.024, 0.043, 0.042 |
| `P6.contacts_vs_core` | rederived | confirmed | +0.069, +0.092, +0.073 |
| `P6.module_map` | rederived | confirmed | all twelve module spans identical (the literature core carried by our alignment) |
| `P6.module_contrast` | rederived | confirmed | core − pore −0.0239 / −0.0005 / −0.0199; ITPR2 p = 0.134 |
| `P6.loop_reverses` | rederived | confirmed | loop counted as pore: +0.055 / +0.044 / +0.056 |
| `P6.shell_distances` | recomputed | confirmed | 125 residues, 12/14/40/59 per shell, medians to 2×10⁻⁵ Å |
| `P6.shell_constraint` | rederived | confirmed | all 12 shell means above the protein; every field identical |
| `P6.shell_trend` | rederived | confirmed | ρ −0.175 / −0.436 / −0.168; significant in ITPR2 only |
| `P6.no_contact_step` | rederived | confirmed | 4.5 Å drop −0.007 / +0.003 / +0.002; largest drop at 11.5 Å in all three |
| `P6.shell_rates` | rederived | confirmed | FEL purifying share 1.000 at contact in all three, fourth above third; all 12 rows, q over the family |
| `P6.module_rates` | rederived | confirmed | β: ITPR1 pore faster (q 6.2e-08), ITPR2/ITPR3 n.s.; all 9 rows |
| `P2.sister_pair` | rederived | confirmed | ITPR2 + ITPR3 |
| `P2.au_test` | read | confirmed | only H3_23 retained (p_AU 0.48) |
| `P2.teleost_itpr1` | rederived | confirmed | ≥ 2 copies: ITPR1 87 %, ITPR2 4 %, ITPR3 2 % |
| `P2.paralog_clades` | rederived | confirmed | 19/13/19 at 100/100, 4/2/1 unnamed; only the 6 cyclostomes outside |
| `P2.cyclostome_lineages` | rederived | confirmed | two cyclostome-only clades, both species in each; 4 first among vertebrates (99.5/100), 2 join ITPR2+ITPR3 (83.1/77) |
| `P2.support_bar` | rederived | confirmed | 91 of 131 bipartitions clear both bars |
| `P2.bnni_robustness` | rederived | confirmed | 13 claim sets rebuilt = claim_members.tsv; 9 of 10 clades held, ITPR1 core 47.8/95 → 47.5/73 |
| `P3.no_absent_cells` | rederived | confirmed | 927 cells, 0 absent |
| `P3.dollo_zero` | read | confirmed | 0 losses |
| `P3.false_negatives` | rederived | confirmed | 140/923, 5/563 contiguous |
| `P4.unreachable` | rederived | confirmed | 744 / 923 |
| `P3.miss_by_contiguity` | rederived | confirmed | median N50 23,460 vs 3,396,515; chromosome 3/512 and 0/172; above bar 189 genomes, 5/563; below 0.3917/0.4333/0.3000 |
| `P3.contiguity_tests` | rederived | confirmed | all 12 tests: OR per 10× 8.10 / 19.99, ITPR vs RyR Fisher p 0.578 |
| `P4.recovery_channels` | rederived | confirmed | 257/309, 260/307, 227/307, RyR 196/309; reasons 289/186/254/15, 179 reachable |
| `P1.presence_range` | rederived | confirmed | 662/6,928 proteomes, 45/135 clades; Archaea 0/634, Bacteria 0/3,537 |
| `P1.kingdom_absences` | rederived | confirmed | Streptophyta 0/384, Ascomycota 0/1034, Basidiomycota 0/319, Apicomplexa 0/60 |
| `P1.relaxed_controls` | rederived | confirmed | all 24 cells; no full-length ITPR or RyR in land plants or Dikarya |
| `P1.absence_targets` | rederived | confirmed | 35 targets by rule G3 at ≥ 10 proteomes |
| `P1.absences` | rederived | confirmed | 35 clades |
| `P1.copy_number` | rederived | confirmed | 194 genomes; Macrostomum 18, Stentor 13, Dysidea 8, Cymbomonas 3 |
| `P1.record_chase` | rederived | confirmed | 64 plant + 35 fungal records: 47 real, 52 fragments |
| `P1.bait_margin` | recomputed | confirmed | pairwise, from the sequences: 31/31 calls, gap 0.582, iplA +0.048 inside ±0.10 |
| `P1.benchmark_counts` | rederived | confirmed | every control's score rebuilt; recall 24/25, specificity 31/31, RyR 6/6 |
| `LEDGER.claims` | read | confirmed | 287 + 245 + 713 ledger rows all ok |

## What the findings checks establish

| kind | what agreement means |
|---|---|
| `recomputed` | the published number follows from the deposited coordinates by a route sharing no code with the publication |
| `rederived` | the published output follows from the published inputs by independent arithmetic; says nothing about the inputs |
| `read` | the prose states what the table says |

A check is only trusted after it has been shown to flip on a planted change
(`tests/test_checks_calibration.py`).
