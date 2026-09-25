# CLAUDE.md — IP3R Structural Simulator

@INTERFACE.md

An interactive, physics-driven 3-D model of the **IP3 receptor** (ITPR1/2/3,
the ER's IP3- and Ca²⁺-gated Ca²⁺-release channel, a C4 homotetramer of
~2,700-residue subunits), and an instrument that **illustrates, demonstrates
and re-derives the results of the `ip3r_genes` publication project**
(`../ip3r_genes`). Ported from the PIEZO1 simulator (`../piezo1_simulation`),
whose `INTERFACE.md` maps a reference implementation for most features still
on the roadmap. Port the method; never port a result.

## Session protocol (every session)

1. **Orient.** Read `INTERFACE.md`, then `ROADMAP.md` (its **Next:** line
   names the task unless the user asks otherwise; Round 6 is in
   `ROADMAP_RYR.md`), then the last entry
   of `SESSION_LOG.md`.
2. **Sync.** `git pull` here and in `../ip3r_genes`. Then `make sync-check`:
   if any `ip3r_genes` source table changed since the last import, run
   `make sync`, then `make checks`, and record in `SESSION_LOG.md` which
   verdicts moved. A moved verdict is a finding about the publication
   project — report it to the user; **never edit `ip3r_genes` from here**.
3. **Environment.** `conda activate ip3r_sim` (every `make` target does this
   itself). Structures: `make fetch` (idempotent, ~27 MB into `ref/`).
4. **Work** the one task. Before finishing: `make test`, `make lint`,
   `make sizes`, and — after any change under `ip3r/ui/` — `make screenshots`
   (the scripted GUI smoke test; it exits non-zero on a broken panel).
5. **Close.** Tick the round in `ROADMAP.md` with the result it measured;
   append to `SESSION_LOG.md` (what, and *why*); update `INTERFACE.md` for any
   structural change and `README.md` for any user-visible one; commit with a
   message explaining why; `git push`.

## Rules that are not optional

- **Files stay under 500 lines** (`tests/test_sizes.py` enforces it). Split
  before you exceed it.
- **Every number a calculation depends on is a registered parameter**:
  declare it in `scripts/parameter_table.py` with unit, bounds, kind,
  citation (a key in `scripts/reference_table.py`, or a sentinel with a
  `source_note` saying why) and rebuild with `make params`. Consume it as
  `_P.value("key")` at call time.
- **A checking instrument is a measuring instrument; calibrate it first.**
  Every findings check needs an entry in `PLANTS` in
  `tests/test_checks_calibration.py` — a planted change to one of its
  *declared* source tables that flips its verdict. The same test proves the
  check's `sources` list is complete (it runs on a copy of only those files).
  When a check disagrees with the publication, **suspect the checker first**:
  on the first run of this project three of five discrepancies were checker
  bugs (variants counted instead of positions, the RyR cell counted as a
  paralog, hydrogens excluded where S22 included them).
- **Check kinds are honest**: `recomputed` (from coordinates, no shared
  code), `rederived` (from the publication's input tables with our code),
  `read` (the table is read and the prose tested against it). Do not label a
  `read` check anything stronger.
- **Residue numbers** are canonical human UniProt numbers of the named
  paralog (Q14643 / Q14571 / Q14573), exactly as in the ip3r_genes S17
  tables. Residue-keyed annotation (elements, conservation, variants) is
  painted on a deposit only when `structure.numbering` says the deposit is in
  that numbering; otherwise grey, never a guess. Rat 7LHF fits no human
  numbering.
- **Missing values are grey, never the low end of a scale.** Quantitative
  colour scales are fixed, never auto-ranged per structure.
- **Physics modules never import `render` or `ui`.** Dependency order:
  `io → core → structure → physics → analysis`; `render` and `ui` consume.
- **Nothing slow on the Qt main thread** — use `ui.workers.run_async`.
- **Nothing downloaded is committed.** `ref/` and `data/` are git-ignored;
  `ip3r/resources/*.json` are committed (curated, with source SHA-256).

## Running

```
python -m ip3r                 # GUI (loads 6DQN)
python -m ip3r checks          # re-derive the ip3r_genes findings
python -m ip3r states          # pore radius of every ITPR3 gating state
python -m ip3r info 8TKF       # measure one deposit
make help                      # everything else
```
