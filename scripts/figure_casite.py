"""Round 7.27's site search as one figure, from ``python -m ip3r casite``'s
own output (``data/casite/cs_<pdb>.out``): P_Ca:P_K at Vais 2010's reversal
against the site's depth, one line per kind of site, beside the uncharged
pore and the measured 15.2. Drawn in the findings panel's colours.

    python scripts/figure_casite.py [OUT.png]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ip3r.analysis.exhibits import dark  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
KINDS = (("uncompensated site", "", "#8a8f99"),
         ("compensated site", "", "#5b9ef2"),
         ("uncompensated site", ", K+ block 1", "#cc80e6"),
         ("compensated site", ", K+ block 1", "#f28c4d"))
_ROW = re.compile(r"^\s+((?:un)?compensated site) ([\d.]+) kT at span(, K\+ block 1)?"
                  r"\s+([\d.]+)\s+([\d.]+)")
_HEAD = re.compile(r"measured P_Cl:P_K ([\d.]+), P_Ca:P_K ([\d.]+); uncharged "
                   r"pore ([\d.]+), ([\d.]+)")
_ROOT = re.compile(r"P_Ca:P_K = ([\d.]+) at d = ([\d.]+) kT: compensated site .*K\+ block")


def parse(path: Path) -> dict:
    out = {"rows": {}, "root": None}
    for line in path.read_text().splitlines():
        if m := _HEAD.search(line):
            out["measured"], out["neutral"] = float(m[2]), float(m[4])
        elif m := _ROW.match(line):
            key = (m[1], m[3] or "")
            out["rows"].setdefault(key, []).append((float(m[2]), float(m[5])))
        elif m := _ROOT.search(line):
            out["root"] = (float(m[2]), float(m[1]))
    return out


def main(out: Path) -> int:
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True)
    for ax, pdb in zip(axes, ("8TKF", "7T3T")):
        d = parse(ROOT / f"data/casite/cs_{pdb.lower()}.out")
        for kind, block, col in KINDS:
            pts = d["rows"][(kind, block)]
            ax.plot(*zip(*pts), "o-", color=col, ms=3, lw=1.3,
                    label=kind + (" + K+ block" if block else ""))
        ax.axhline(d["measured"], color="#72cc80", lw=1.0, ls="--",
                   label=f"measured (Vais 2010) {d['measured']:g}")
        ax.axhline(d["neutral"], color="#d7dbe3", lw=0.8, ls=":",
                   label="uncharged pore")
        ax.annotate(f"uncharged pore {d['neutral']:.2f}", (12, d["neutral"]),
                    xytext=(0, 3), textcoords="offset points", ha="right",
                    color="#d7dbe3", fontsize=7)
        if d["root"]:
            ax.plot(*d["root"], "*", color="#f2cc4d", ms=11)
            ax.annotate(f"crosses at {d['root'][0]:.2f} kT", d["root"],
                        xytext=(8, -12), textcoords="offset points",
                        color="#f2cc4d", fontsize=7)
        ax.set_yscale("log")
        ax.set_xlabel("site depth d (kT; K_d = s / (e^d - 1))")
        dark(fig, ax)
        ax.set_title(f"{pdb}: P_Ca:P_K at bi-ionic reversal", fontsize=9)
    axes[0].set_ylabel("P_Ca:P_K")
    handles, labels = axes[0].get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7,
                     frameon=False)
    for t in leg.get_texts():
        t.set_color("#d7dbe3")
    fig.tight_layout(rect=(0, 0.14, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]) if len(sys.argv) > 1
                          else ROOT / "docs/img/casite_ratio.png"))
