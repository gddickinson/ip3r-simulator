"""The Analyses menu's catalogue, Qt-free (Round 7.18).

Most of Round 7's science (the 3-D shortfall, the wall in 3-D, the salt
bridge, the image cost, charge–space competition, selectivity and
protonation) had no way in from the GUI: it was command-line only. PIEZO1
closed the same gap with a result window per analysis; here each entry is
the CLI command itself, run in a child process (``python -m ip3r ...``),
so the window prints exactly what the command-line prints and nothing is
implemented twice. The child process also keeps the minutes-long solves
off the GUI's interpreter entirely, and a Stop button can end one.

``{pdb}`` is the deposit on screen; ``{paralog}`` its numbering (ITPR3 when
it has none). An entry that needs a deposit is disabled without one.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Analysis", "ANALYSES", "GROUPS", "command", "stamp"]


@dataclass(frozen=True)
class Analysis:
    key: str
    group: str
    label: str
    args: tuple[str, ...]
    #: How long it takes, said before it starts.
    duration: str
    #: What it computes (the menu tooltip and the window's first line).
    about: str

    @property
    def needs_deposit(self) -> bool:
        return any("{pdb}" in a for a in self.args)


GROUPS = ("Structure", "Permeation", "Gating and puffs", "Findings")

ANALYSES = (
    Analysis("info", "Structure", "Measure this deposit", ("info", "{pdb}"), "seconds",
             "C4 axis two ways, numbering, pore profile, constrictions, IP3 contacts."),
    Analysis("states", "Structure", "Pore radius of every state", ("states", "--paralog",
             "{paralog}"), "under a minute", "Every human deposit of the paralog, measured alike."),
    Analysis("modes", "Structure", "Normal modes of this deposit", ("modes", "{pdb}"),
             "a minute", "C4-labelled ANM: irreps and collectivity of the lowest modes."),
    Analysis("cutoff", "Structure", "Transition 8TKG → 8TKF: cutoff scan",
             ("transition", "--cutoff-scan"), "a few minutes",
             "The lowest collective A mode's overlap against the ANM cutoff (Round 2)."),
    Analysis("graft", "Structure", "AlphaFold fill of this deposit", ("graft", "{pdb}"),
             "a minute", "Each unresolved stretch filled from the prediction: anchor RMSD, "
             "seams, pLDDT, clashes; skips named."),
    Analysis("unitary", "Permeation", "K+ conductance of every state",
             ("unitary", "--paralog", "{paralog}"), "a minute",
             "1-D drift-diffusion: neutral, charged and paired walls vs Mak / Vais."),
    Analysis("selectivity", "Permeation", "Selectivity (Vais 2010's protocols)",
             ("selectivity", "{pdb}"), "a minute",
             "P_Cl:P_K, P_Ca:P_K and i_Ca on every wall reading."),
    Analysis("protonation", "Permeation", "Protonation of the lining",
             ("protonation", "{pdb}"), "a few minutes",
             "Lining pKas by the Tanford–Kirkwood network and PROPKA; selectivity under each."),
    Analysis("shortfall", "Permeation", "3-D shortfall of the open deposits",
             ("shortfall",), "several minutes",
             "1-D, lumen-area 1-D and 3-D conductance of 8TKF, 7T3T, 9HEO vs measured."),
    Analysis("lumen", "Permeation", "Where the voltage falls", ("lumen", "{pdb}"),
             "a minute", "3-D vs 1-D drop and each constriction's share (neutral)."),
    Analysis("wall3d", "Permeation", "Wall charge in 3-D", ("wall3d", "{pdb}"),
             "several minutes", "Neutral and each closure (slice / local / pb), paired and not."),
    Analysis("bridge", "Permeation", "The lining salt bridge", ("bridge", "{pdb}"),
             "several minutes", "The bridge's pKas three ways, then its field under the "
             "dielectric closure (Round 7.13)."),
    Analysis("born", "Permeation", "Image cost of the low-ε wall", ("born", "{pdb}"),
             "seconds cached; a quarter hour on every core uncached",
             "W on the axis, K+ g with it, the dipole with and without it (Round 7.15)."),
    Analysis("csc", "Permeation", "Charge–space competition", ("csc", "{pdb}"),
             "a few minutes", "Filter binding split as Gillespie's Fig. 7; P_Ca:P_K under "
             "donnan / csc (Round 7.17)."),
    Analysis("gating", "Gating and puffs", "Gating bells, three models",
             ("gating", "--model", "pd"), "seconds",
             "Park/drive open probability against Ca2+ at four IP3 levels."),
    Analysis("oscillate", "Gating and puffs", "Oscillation window",
             ("oscillate", "--window"), "under a minute",
             "Li–Rinzel cell: the IP3 window that oscillates vs 0.36–0.63 µM measured."),
    Analysis("microdomain", "Gating and puffs", "Puffs in the microdomain",
             ("microdomain",), "a minute",
             "Park/drive cluster with fluo-4 in Cao's microdomain: IPIs, amplitudes."),
    Analysis("checks", "Findings", "Run every findings check", ("checks",), "a few minutes",
             "Re-derive the ip3r_genes results; each verdict and its kind."),
)


def command(a: Analysis, pdb: str | None, paralog: str | None) -> list[str]:
    """The CLI arguments with the deposit filled in (ValueError without one)."""
    if a.needs_deposit and not pdb:
        raise ValueError(f"{a.label} needs a deposit on screen")
    return [s.replace("{pdb}", pdb or "").replace("{paralog}", paralog or "ITPR3")
            for s in a.args]


def stamp(args: list[str], pdb: str | None, overrides: dict) -> str:
    """The header every result carries: the command, the deposit on screen,
    and the parameter set it ran under (a result without it is unreadable
    once a parameter has been edited)."""
    params = ("registered defaults" if not overrides else
              f"{len(overrides)} override(s): " +
              ", ".join(f"{k} = {v:g}" for k, v in sorted(overrides.items())))
    return (f"$ python -m ip3r {' '.join(args)}\n"
            f"# deposit on screen: {pdb or 'none'}; parameters: {params}\n")
