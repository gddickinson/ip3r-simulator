"""The in-application guide as data (Help → Guide, F1; Round 7.18).

Topics are HTML fragments; ``SHORTCUTS`` is the one table of keys, which
the menus' shortcuts are tested against, so the guide cannot list a key
the window does not have. ``DOCS`` names the shipped documents Help can
open.
"""

from __future__ import annotations

__all__ = ["TOPICS", "SHORTCUTS", "DOCS"]

#: (keys, what they do). Menu shortcuts as Qt writes them.
SHORTCUTS = (
    ("Ctrl+O", "open a saved session"),
    ("Ctrl+Shift+S", "save the session (the view, never results)"),
    ("Ctrl+S", "save a screenshot (with the HUD)"),
    ("Ctrl+1", "side view, cytosol up"),
    ("Ctrl+2", "top view, down the pore"),
    ("Ctrl+3", "the first IP3 site on a visible subunit"),
    ("Ctrl+0", "fit the visible subunits (resumes the automatic fit)"),
    ("Ctrl+Shift+Q", "the sequence window"),
    ("F11", "full screen: only the 3-D view and its HUD (also the "
            "platform's own full-screen key)"),
    ("Ctrl+M", "measure distances: click two atoms"),
    ("Esc", "leave full screen; otherwise clear the selection and stop measuring"),
    ("Ctrl+Shift+P", "the parameter editor"),
    ("F1", "this guide"),
    ("Space", "spin on / off (viewport focused)"),
    ("R", "frame everything (viewport focused)"),
    ("O", "orthographic / perspective (viewport focused)"),
    ("+ / −", "atom size (viewport focused)"),
)

_MOUSE = """<table>
<tr><td>left-drag</td><td>rotate</td></tr>
<tr><td>shift-drag, middle-drag</td><td>pan</td></tr>
<tr><td>wheel, right-drag</td><td>zoom</td></tr>
<tr><td>click</td><td>select the residue; the status line names it</td></tr>
<tr><td>shift-click</td><td>add or remove a residue</td></tr>
<tr><td>right-click</td><td>the context menu for what is under the cursor</td></tr>
</table>"""

TOPICS = {
    "Getting around": f"""
<p>The viewport shows one deposit, framed side-on with the cytosol up. The
camera refits itself on resize and when subunits are hidden, until you move
it; <b>Ctrl+0</b> resumes the fit.</p>{_MOUSE}
<p>The HUD (View → HUD) carries a scale bar (in perspective exact only in the
plane through the pivot, and labelled so; in orthographic, <b>O</b>, exact
at every depth), a gnomon with the four-fold axis's cytosolic end (<b>cyt</b>), and
the selection and last distance. An amber line appears whenever an
AlphaFold prediction is drawn: it is not switchable.</p>
<p>Panels are docks: drag them anywhere, float or close them; View →
Panels brings one back and View → Reset layout restores the shipped
arrangement. The layout is remembered between runs.</p>
<p><b>Full screen</b> (View → Full screen, F11, or the right-click menu)
leaves only the 3-D view and its HUD; the parameters banner stays if
parameters are modified. Every shortcut still works. Esc or F11 leaves,
and each panel comes back as it was, closed ones closed.</p>""",

    "Selection and measuring": """
<p>A click selects a residue, drawn as gold spheres in their own layer, so
it never disturbs the Structure panel's site highlights. The right-click
menu selects a residue on every subunit (the C4 copies), a whole chain,
centres on an atom, or starts a distance.</p>
<p><b>Ctrl+M</b> arms measuring: every two clicks close a distance, drawn as
a rod labelled in Å. Selections and distances follow a morph or mode
frame, and are dropped on a new deposit: a residue number means nothing on
another structure.</p>""",

    "Sequence window": """
<p>View → Sequence shows one chain's <i>construct</i>, every residue the
mmCIF lists, in the deposit's own author numbers. Unresolved residues are
kept, dimmed: a gap is information. Tracks paint the functional element,
the deep JSD on the fixed 0.50–0.95 scale, or which residues are resolved;
the IP3 contacts and the filter and gate lining are underlined.</p>
<p>Residue-keyed tracks paint only when the deposit is in that paralog's
numbering; rat 7LHF fits none, so it is grey rather than guessed. A drag
selects onto the model (tick "every subunit" for the C4 copies), and a
selection on the model shows here.</p>""",

    "Colourings and numbering": """
<p>Every quantitative colour scale is fixed, never auto-ranged per
structure, and a missing value is grey, never the low end of a scale.
Element, conservation, variants and shells are painted only on a deposit
in human numbering of the named paralog (Q14643 / Q14571 / Q14573).</p>
<p>Completeness adds AlphaFold's residues where the deposit has none,
coloured by pLDDT band; seams are drawn, red when broken.</p>
<p>Superpose draws another state of the same paralog on the one shown, in
orange, by the Transition tab's residue-matched fit (pore or global); the
panel gives the RMSD over the fitted sites and over all. Colour the shown
deposit Uniform to read the pair. It stays at the deposit's coordinates
through a morph, a fixed reference.</p>""",

    "Channel and the lumen": """
<p>The Channel tab measures the pore (S0's quantity), compares the states
of the paralog, and solves each state's K<sup>+</sup> conductance.
<b>Draw the lumen</b> voxelises the ion-accessible volume and solves where
the voltage falls, in 3-D and in the 1-D model on the same window.</p>
<p>Wall charge: none, Round 7.11's slice / local / pb placements, or the
dielectric closure (every charged group behind a protein of low ε).
<b>+ image</b> adds the image (Born) cost W of that wall. Colour by: the
voltage drop; the wall potential u; the <b>K<sup>+</sup> energy u + W</b>,
the well a cation actually feels (Round 7.18; u alone without the image);
or W itself. Scales are fixed (± display.lumen_potential_range kT,
0–display.lumen_image_range kT).</p>
<p><b>Steady state → at reversal</b> (Round 7.24) replaces the wall charge
with Round 7.23's bi-ionic experiment (Ca<sup>2+</sup>, or Cl<sup>−</sup>
for IP3R) solved to the voltage where no net current flows, under the
neutral, pb or pb + csc reading, on the 1 Å reversal grid (a minute or
two). Colour by each ion's <b>concentration</b> (fixed log scale,
display.lumen_conc_min–display.lumen_conc_max M) or its
<b>electrochemical drop</b> (0 lumen, 1 cytosol: where it rises steeply is
where that ion's resistance lies). The plot gives each ion's concentration
against its baths and its drop beside the neutral pore's.</p>
<p>The same selector offers Round 7.25's <b>candidate walls</b> (Round 7.26,
IP3R only), read as the search read them (point ions under Poisson): a
Ca<sup>2+</sup>-only well of wallsearch.gui_well_depth kT over the span,
that well with the deposit's charge, or the opposite-charge C4 ring pair
with the highest B (two minutes more the first time, for its search). Draw
the deposit's own reading first: the plot then sets it dashed beside the
candidate, with both reversal potentials and peaks in the text.</p>
<p>The last candidate is Round 7.27's <b>Ca<sup>2+</sup> site</b> (Round
7.28): four saturable sites over the span at casite.gui_depth kT,
compensated as they fill, whose occupancy blocks K<sup>+</sup> (the depth
at which it gives Vais's 15.2 on 8TKF). Colour it by the site's
<b>occupancy θ</b> (0–1) or by <b>K<sup>+</sup> block</b> −ln(1 − fθ)
(0–display.lumen_block_range kT); both are grey for any other reading.
The plot adds θ's plane mean, dotted, on the drop row.</p>""",

    "Analyses menu": """
<p>Analyses runs the command-line science in its own process, the same
command you would type (<code>python -m ip3r …</code>), so the window
prints exactly what the CLI prints. The command line is editable before
Run; <code>--help</code> lists a command's flags. Each result starts with
its stamp: the command, the deposit on screen, and the parameter set it
ran under. Stop ends it; Save keeps the text.</p>""",

    "Findings": """
<p>Each check re-derives one ip3r_genes result and says how:
<b>recomputed</b> (from coordinates, no shared code), <b>rederived</b>
(from the publication's input tables with this project's code) or
<b>read</b> (the table read and the prose tested against it). Every check
is calibrated by a planted change that flips its verdict. Checks refuse to
confirm under modified parameters.</p>""",

    "Parameters and sessions": """
<p>Every number a calculation uses is a registered parameter with a unit,
bounds and a citation (Help → Parameters). Editing one shows the amber
banner; results made under an edit are dropped when it is reset. A session
saves the view (deposit, style, panels, camera) and the parameter set it
was saved under, never a result: opening it re-measures.</p>""",
}

#: Shipped documents Help can open: (menu label, path from the project root).
DOCS = (
    ("README", "README.md"),
    ("Navigation map (INTERFACE.md)", "INTERFACE.md"),
    ("The science", "docs/SCIENCE.md"),
    ("The findings checks", "docs/SCIENCE_CHECKS.md"),
    ("Permeation", "docs/SCIENCE_PERM.md"),
    ("The image cost", "docs/SCIENCE_BORN.md"),
    ("Charge–space competition", "docs/SCIENCE_CSC.md"),
    ("Roadmap", "ROADMAP.md"),
)
