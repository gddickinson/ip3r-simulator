# Contributing to the IP3R Structural Simulator

This guide is for anyone changing the code. Start with the
[README](README.md) for what the project is, and read
[`INTERFACE.md`](INTERFACE.md) before opening any source file: it maps
every module, class and function.

## The repository is organised into layers that depend in one direction.

```
io ──▶ core ──▶ structure ──▶ physics ──▶ analysis
                                  │
                   render ◀───────┴───────▶ ui
```

The science (everything left of `render` and `ui`) never imports the viewer,
so it runs headless in the command line, tests and notebooks.

| Path | What it holds |
|---|---|
| `ip3r/io/` | reading structure files, the structure registry, downloads |
| `ip3r/core/` | the structure container, annotations, pairwise alignment, access to `ip3r_genes` |
| `ip3r/structure/` | measurement: symmetry axis, pore profile, ligand contacts, morphs, AlphaFold fills |
| `ip3r/physics/` | the models: normal modes, gating, puffs and sparks, permeation in 1-D and 3-D |
| `ip3r/analysis/` | the 52 findings checks and the statistics and tree code they use |
| `ip3r/render/` | the OpenGL renderer and colour scales |
| `ip3r/ui/` | the PyQt6 desktop app |
| `ip3r/resources/` | curated data (committed, each with the SHA-256 of its source) |
| `scripts/` | the parameter table, data import and curation, figure scripts, the GUI smoke test |
| `tests/` | the pytest suite |
| `docs/` | the science documents and screenshots |
| `ref/`, `data/` | downloads and derived output (not committed; regenerable) |

## The test suite checks the science, the checks and the interface.

```bash
make test          # the full suite, about 740 tests (tests needing missing data are skipped)
make test-quick    # only the tests that need no structures and no ip3r_genes
make lint          # static checks with ruff
make sizes         # fails if any Python file exceeds 500 lines
make screenshots   # drives the real app through every panel and rewrites docs/img/
```

The tests calibrate each model against a case with a known answer before
using it on a real structure. Examples include textbook conductance of a
cylinder, the Debye–Hückel limit, and a planted flip in each findings check.
The GUI smoke test (`make screenshots`) exits with an error if any panel
fails. Its steps are grouped so that one area can be checked on its own,
for example `make screenshots STEPS=lumen,extras` (`make screenshot-groups`
lists the groups).

## Contributions follow a small set of project rules.

- **Every number is a registered parameter.** Declare it in
  `scripts/parameter_table*.py` with a unit, bounds and a citation, then
  rebuild with `make params`. Code reads it at run time with
  `_P.value("key")`.
- **Every findings check is calibrated.** Add a planted change to
  `PLANTS` in `tests/test_checks_calibration.py` that flips its verdict.
- **Check labels are honest.** A check that only reads a table is labelled
  `read`, never `rederived` or `recomputed`.
- **Missing values are grey.** Colour scales are fixed and never stretched to
  fit one structure.
- **Nothing slow runs on the interface thread.** Use `ui.workers.run_async`.
- **Files stay under 500 lines.** `make sizes` enforces this.
- **Nothing downloaded is committed.** `ref/` and `data/` are ignored by git.
- **The publication project is never edited from here.** Discrepancies are
  reported to it instead.

Before submitting a change, run `make test`, `make lint`, `make sizes` and,
for any change under `ip3r/ui/`, the relevant `make screenshots` groups.

## Each working session follows the same protocol.

The project is developed in rounds of one focused task. At the start of a
session, read `INTERFACE.md`, the **Next:** line of `ROADMAP.md` and the last
entry of `SESSION_LOG.md`, then `git pull` here and in `../ip3r_genes` and
run `make sync-check`. If a publication table changed, run `make sync` and
`make checks` and record which verdicts moved. At the end, tick the round in
`ROADMAP.md` with what it measured, add an entry to `SESSION_LOG.md`
explaining what was done and why, update `INTERFACE.md` and the README for
any structural or user-visible change, and commit with a message that
explains the reason for the change.
