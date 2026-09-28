"""Round 7.29's mole-fraction prediction as one figure, from
``python -m ip3r molefrac``'s own output (``data/molefrac/mf_<pdb>_*.out``):
against luminal Ca²⁺, P_Ca:P_K at Vais 2010's reversal, the K⁺ current at
the hold relative to the Ca²⁺-free pore's, and the site's occupancy. One
colour per reading, one line style per deposit. Drawn in the findings
panel's colours.

    python scripts/figure_molefrac.py [OUT.png]
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
COLOURS = {"site": "#f28c4d", "site, no block": "#5b9ef2",
           "uncharged": "#8a8f99"}
STYLES = {"8TKF": "-", "7T3T": "--"}
_PDB = re.compile(r"^(\w{4}): band '([^']+)', hold")
_HEAD = re.compile(r"^  (site|site, no block|uncharged)(?::|;)")
_ROW = re.compile(r"^\s+([\d.e+-]+)\s+([+-][\d.]+)\s+([\d.-]+)\s+([\d.e+-]+)"
                  r"\s+([\d.e+-]+)\s+([\d.-]+)\s+([\d.-]+)\s+([\d.]+)\s+([\d.]+)")


def parse(path: Path) -> dict:
    """{(pdb, reading): [(Ca mM, P_Ca:P_K, i_K/0, theta), ...]}"""
    out, pdb, reading = {}, None, None
    for line in path.read_text().splitlines():
        if m := _PDB.match(line):
            pdb = m[1]
        elif m := _HEAD.match(line):
            reading = m[1]
        elif (m := _ROW.match(line)) and pdb and reading:
            out.setdefault((pdb, reading), []).append(
                (float(m[1]), float(m[3]), float(m[6]), float(m[8])))
    return out


def main(out: Path) -> int:
    data = {}
    for f in sorted((ROOT / "data/molefrac").glob("mf_*.out")):
        data.update(parse(f))
    if not data:
        print("no molefrac output under data/molefrac")
        return 1
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.9))
    for (pdb, reading), rows in sorted(data.items()):
        ca, pca, kf, th = zip(*rows)
        kw = dict(color=COLOURS[reading], ls=STYLES.get(pdb, ":"), lw=1.3,
                  marker="o", ms=3, label=f"{pdb} {reading}")
        axes[0].plot(ca, pca, **kw)
        axes[1].plot(ca, kf, **kw)
        if reading != "uncharged":
            axes[2].plot(ca, th, **kw)
    axes[0].axhline(15.2, color="#72cc80", lw=1.0, ls="--")
    axes[0].annotate("measured 15.2 at 10 mM", (10, 15.2), xytext=(0, 4),
                     textcoords="offset points", ha="center",
                     color="#72cc80", fontsize=7)
    axes[0].set_yscale("log")
    for ax, title, ylab in zip(
            axes, ("P_Ca:P_K at bi-ionic reversal",
                   "K+ current at the hold", "the site's occupancy"),
            ("P_Ca:P_K", "i_K / i_K with no Ca2+", "mean theta over the band")):
        ax.set_xscale("log")
        ax.set_xlabel("luminal CaCl2 (mM), 140 mM KCl both sides")
        ax.set_ylabel(ylab)
        dark(fig, ax)
        ax.set_title(title, fontsize=9)
        ax.axvline(10, color="#d7dbe3", lw=0.6, ls=":")
    axes[1].set_ylim(0, 1.1)
    axes[2].set_ylim(0, 1.05)
    handles, labels = axes[0].get_legend_handles_labels()
    leg = fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7,
                     frameon=False)
    for t in leg.get_texts():
        t.set_color("#d7dbe3")
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150, facecolor=fig.get_facecolor())
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(Path(sys.argv[1]) if len(sys.argv) > 1
                          else ROOT / "docs/img/molefrac.png"))
