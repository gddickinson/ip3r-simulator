"""The ITPR3 state panel at the pore, as one figure: every deposit of
ip3r_genes S11's panel measured by :func:`ip3r.structure.states.state_panel`
(one axis method, one profile, one constriction rule), plus the open-state
control 7T3T. (a) Each pore profile, z measured from that deposit's filter so
the constrictions line up; (b) the gate and filter radius of each, ordered by
the gate. Drawn in the findings panel's colours.

    python scripts/figure_states.py [OUT.png]
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from ip3r.analysis.exhibits import dark  # noqa: E402
from ip3r.io import loader  # noqa: E402
from ip3r.io.registry import get_entry  # noqa: E402
from ip3r.structure.channel import measure_channel  # noqa: E402
from ip3r.structure.states import StateRow, state_panel  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CONTROL = "7T3T"
_C = ["#5b9ef2", "#f28c4d", "#72cc80", "#cc80e6", "#f2cc4d", "#66d9d9", "#e06c6c"]
_FG = "#d7dbe3"


def rows() -> list[StateRow]:
    panel = state_panel("ITPR3")
    e = get_entry(CONTROL)
    ctrl = StateRow(CONTROL, e.state, e.resolution, e.ip3_bound,
                    measure_channel(loader.load(CONTROL)))
    return panel + [ctrl]


def main(out: Path) -> int:
    data = rows()
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 5.0),
                                 gridspec_kw={"width_ratios": [1.5, 1]})
    for k, r in enumerate(data):
        p = r.summary.profile
        zf = r.summary.constrictions["filter"].z
        control = r.pdb_id == CONTROL
        col = "#8a8f99" if control else _C[k % len(_C)]
        wide = r.state.split()[0] in ("activated", "active", "open")
        ax.plot(np.asarray(p.z) - zf, p.r_min, color=col,
                lw=1.8 if wide else 0.9, ls="--" if control else "-",
                label=f"{r.pdb_id} {r.state.split('(')[0].strip()}"
                      + (" (control)" if control else ""))
    ax.set_xlim(-25, 40)
    ax.set_ylim(0, 14)
    ax.set_xlabel("z from the filter (Å; luminal ← → cytosolic)")
    ax.set_ylabel("min heavy-atom radius (Å)")
    ax.set_title("(a) pore profiles, aligned on the filter", fontsize=9)
    ax.axvline(0, color="#4a505c", lw=0.8, ls=":")
    ax.annotate("filter", (0, 0.4), xytext=(3, 0), textcoords="offset points",
                color=_FG, fontsize=7)
    handles, labels = ax.get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=7,
                     frameon=False)
    for t in leg.get_texts():
        t.set_color(_FG)

    order = sorted(data, key=lambda r: r.radius("gate"))
    y = np.arange(len(order))
    bx.scatter([r.radius("gate") for r in order], y, color="#f28c4d", s=28,
               label="gate", zorder=3)
    bx.scatter([r.radius("filter") for r in order], y, color="#5b9ef2", s=28,
               marker="s", label="filter", zorder=3)
    for yi, r in zip(y, order):
        bx.plot([r.radius("gate"), r.radius("filter")], [yi, yi],
                color="#4a505c", lw=0.8, zorder=1)
    bx.set_yticks(y, [f"{r.pdb_id} {r.state.split('(')[0].strip()}"
                      + (" *" if r.pdb_id == CONTROL else "") for r in order],
                  fontsize=7)
    bx.set_xlabel("radius at the constriction (Å)")
    bx.set_title("(b) gate and filter per state", fontsize=9)
    leg = bx.legend(fontsize=7, frameon=False, loc="lower right")
    for t in leg.get_texts():
        t.set_color(_FG)
    for a in (ax, bx):
        dark(fig, a)
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    for r in order:
        print(f"{r.pdb_id}\t{r.state}\tgate {r.radius('gate'):.2f}\t"
              f"filter {r.radius('filter'):.2f}")
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]) if len(sys.argv) > 1
                          else ROOT / "docs/img/state_panel_pore.png"))
