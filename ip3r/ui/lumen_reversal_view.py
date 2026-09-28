"""The lumen box's plot at reversal (Round 7.24).

Three rows on S0's window, the reversal grid's lumen:

- the lumen's area, 3-D against the 1-D inscribed circle (as Round 7.10);
- each ion's plane-mean concentration on a log scale, its two baths marked
  at the ends, so a well (Ca²⁺ in the filter) or an exclusion (Cl⁻ in a
  charged wall) reads against the solution it came from;
- each ion's electrochemical drop beside the neutral pore's: where each
  rises steeply is where its resistance lies.

Round 7.26: with a candidate wall, the deposit's own reading of the same
experiment (``beside``) is drawn dashed on both ion rows, and the text sets
the two reversal potentials and each ion's peak side by side.

Round 7.28: with a Ca²⁺ site, its plane-mean occupancy θ is drawn on the
drop row (the same 0–1 scale), dotted.

Round 7.30: ``earlier`` is the same reading at another luminal CaCl₂; its
Ca²⁺ plane mean and (with a site) its θ are drawn thin beside the current
ones, and the text sets the two reversals and ions held side by side.
"""

from __future__ import annotations

import numpy as np

from .plot_canvas import PALETTE

__all__ = ["draw_reversal", "ION_COLOURS"]

#: One colour per ion, the same on every row.
ION_COLOURS = {"K+": PALETTE[0], "Cl-": PALETTE[2], "Ca2+": PALETTE[1]}


def draw_reversal(canvas, rev, s, beside=None, earlier=None) -> str:
    """Plot ``rev`` (a :class:`~ip3r.physics.lumen_reversal.ReversalLumen`)
    for the channel summary ``s`` on ``canvas``, with ``beside`` (another on
    the same grid, the deposit's own wall) dashed and ``earlier`` (the same
    reading at another luminal CaCl₂) thin; returns the panel's text."""
    from ..parameters import PARAMETERS as _P
    f = rev.lumen
    axes = canvas.reset(3, 1)
    area, conc, drop = axes[0, 0], axes[1, 0], axes[2, 0]
    area.plot(f.z, f.area_3d, color=PALETTE[0], lw=1.4, label="3-D: lumen region")
    area.plot(f.z, f.area_1d, color=PALETTE[2], lw=1.0,
              label="1-D: π (r_free − r_ion)²")
    area.set_ylabel("area (Å²)")
    at = f", {rev.ca * 1e3:.3g} mM CaCl2" if rev.ca > 0 else ""
    area.set_title(f"{f.name}: {rev.experiment} experiment{at} at reversal "
                   f"({rev.reading}), V = {rev.v * 1e3:+.1f} mV", fontsize=8)
    for ion in rev.species:
        col = ION_COLOURS[ion]
        c = rev.conc_3d(ion)
        conc.plot(f.z, np.where(c > 0, c, np.nan), color=col, lw=1.4,
                  label=f"{ion}: plane mean")
        lum, cyt = rev.baths[ion]
        for z, b in ((f.z[0], lum), (f.z[-1], cyt)):
            if b > 0:
                conc.plot([z], [b], marker="o", ms=4, color=col, ls="none")
        drop.plot(f.z, rev.drop_3d(ion), color=col, lw=1.4, label=ion)
        if beside is not None and ion in beside.conc:
            # one legend entry for every dashed line, the ions share colours
            first = ion == next(i for i in rev.species if i in beside.conc)
            b = beside.conc_3d(ion)
            conc.plot(beside.z, np.where(b > 0, b, np.nan), color=col, lw=0.9,
                      ls="--", label=(f"dashed: deposit ({beside.reading})"
                                      if first else "_nolegend_"))
            drop.plot(beside.z, beside.drop_3d(ion), color=col, lw=0.9,
                      ls="--", label="_nolegend_")
    then = f"{earlier.ca * 1e3:.3g} mM" if earlier is not None else ""
    if earlier is not None and "Ca2+" in earlier.conc:
        b = earlier.conc_3d("Ca2+")
        conc.plot(earlier.z, np.where(b > 0, b, np.nan),
                  color=ION_COLOURS["Ca2+"], lw=0.8, alpha=0.6,
                  label=f"Ca2+ at {then}")
    conc.set_yscale("log")
    conc.set_ylabel("c (M), ● bath")
    drop.plot(f.z, f.drop_3d, color="#8a8f99", lw=1.0, ls="--",
              label="neutral pore (Laplace)")
    occ = rev.occupancy_3d()
    if occ is not None:
        drop.plot(f.z, occ, color="#d7dbe3", lw=1.2, ls=":",
                  label="site occupancy θ (plane mean)")
    occ_then = None if earlier is None else earlier.occupancy_3d()
    if occ_then is not None:
        drop.plot(earlier.z, occ_then, color="#8a8f99", lw=1.0, ls=":",
                  label=f"θ at {then}")
    drop.set_ylabel("ion drop" if occ is None else "ion drop, θ")
    drop.set_ylim(-0.02, 1.02)
    drop.set_xlabel("z along the four-fold axis (Å; luminal ← → cytosolic)")
    w = _P.value("lumen.constriction_half_width")
    for c in s.constrictions.values():
        for ax in axes[:, 0]:
            ax.axvline(c.z, color="#8a8f99", lw=0.6, ls=":")
        drop.annotate(c.name, (c.z, 0.05), xytext=(3, 0),
                      textcoords="offset points", color="#d7dbe3", fontsize=7,
                      xycoords=("data", "axes fraction"))
    for ax in axes[:, 0]:
        canvas.legend(ax, loc="upper left")
    canvas.draw_now()
    lines = [rev.summary() + "."]
    if beside is not None:
        peaks = ", ".join(f"{ion} {rev.peak(ion)[0]:.3g} M against "
                          f"{beside.peak(ion)[0]:.3g} M" for ion in rev.species
                          if ion in beside.conc)
        lines.append(f"Beside the deposit's own wall ({beside.reading}, dashed):"
                     f" V_rev {rev.v * 1e3:+.2f} against {beside.v * 1e3:+.2f} "
                     f"mV; peaks {peaks}.")
    if earlier is not None:
        held = (f"; the site holds {rev.held:.2f} Ca2+ against "
                f"{earlier.held:.2f}" if np.isfinite(rev.held)
                and np.isfinite(earlier.held) else "")
        lines.append(f"Beside the same reading at {then} luminal CaCl2 (thin):"
                     f" V_rev {rev.v * 1e3:+.2f} against {earlier.v * 1e3:+.2f}"
                     f" mV{held}.")
    for c in s.constrictions.values():
        shares = ", ".join(f"{ion} {rev.drop_across(ion, c.z, w):.0%}"
                           for ion in rev.species)
        lines.append(f"{c.name} (z {c.z:+.1f} ± {w:.0f} Å): {shares} of each "
                     "ion's drop")
    note = (f"Round 7.23's steady state on its own {f.volume.spacing:g} Å grid "
            f"({f.species}'s volume, the electrostatic one), not the "
            "equilibrium box's. Concentrations are what Poisson counts; an ion "
            "too large for a voxel takes its nearest own voxel's n there. At "
            "reversal the ions' currents cancel, each still flows.")
    return "<br>".join(lines) + f"<br><i>{note}</i>"
