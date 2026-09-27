"""Round 7.23: the bi-ionic reversal in 3-D — the steady Poisson–Nernst–
Planck solve on the voxels (Gummel's loop) and the reversal root. Method
constants only: the baths, diffusivities and radii are the 1-D protocols'.
Imported into ``parameter_table.P``.
"""

from param_entry import entry as _p

REV3D = [
    _p("pnp3d.gummel_tolerance", "Gummel tolerance", 1e-4, "kT/e", "method",
       "pore3d", "method_choice", "Round 7.23: the loop between Nernst-Planck "
       "and Poisson stops when the potential moves less than this anywhere "
       "in the lumen.", "1e-4 kT/e is 2.6 uV, far below the 0.01 mV the "
       "reversal is read to", 1e-8, 1e-1),
    _p("pnp3d.gummel_damping", "Gummel damping", 1.0, "", "method",
       "pore3d", "method_choice", "Round 7.23: fraction of each Poisson "
       "update taken (1 = undamped).", "Undamped converges on every deposit "
       "and protocol measured; kept as a control for a harder wall", 0.05,
       1.0),
    _p("pnp3d.gummel_max_iterations", "Gummel iteration limit", 200, "",
       "method", "pore3d", "method_choice", "Round 7.23: the loop gives up "
       "(and the reading is flagged n.c.) after this many sweeps.", "A "
       "converged solve takes far fewer; the limit only bounds a failure",
       5, 5000),
    _p("reversal3d.spacing", "Reversal grid spacing", 1.0, "A", "method",
       "pore3d", "method_choice", "Round 7.23: voxel edge for the 3-D "
       "reversal solves (each reversal is a root over several full "
       "Poisson-Nernst-Planck solves).", "Round 7.19 found P_Ca:P_K at 1.0 A "
       "within 4 % of 0.5 A on 9HEO and 8TKF; the CLI's --spacing reads "
       "0.5 A", 0.25, 2.0),
    _p("reversal3d.voltage_tolerance", "Reversal voltage tolerance", 1e-5,
       "V", "method", "pore3d", "method_choice", "Round 7.23: the reversal "
       "root is found to this (0.01 mV).", "As the 1-D reversal's 1e-6 V "
       "order; 0.01 mV moves a ratio of 1 by < 0.1 %", 1e-8, 1e-3),
    _p("reversal3d.bracket", "Reversal search bracket", 0.1, "V", "method",
       "pore3d", "method_choice", "Round 7.23: the reversal is searched "
       "within +/- this.", "As the 1-D search (Round 7.4), which holds "
       "every reversal it has read; a root outside is refused by name", 0.01, 0.3),
    # Round 7.25: the wall-model search (physics/wall_search.py)
    _p("wallsearch.scale_min", "Search: smallest charge scale", 0.01, "",
       "method", "pore3d", "method_choice", "Round 7.25: the deposit's "
       "lining charge is scaled on a geometric grid from this to 1.", "7.23 "
       "found P_Cl:P_K halved already at 0.05; 0.01 reaches back to the "
       "uncharged pore's reading", 1e-4, 0.5),
    _p("wallsearch.scale_points", "Search: charge scales", 5, "", "method",
       "pore3d", "method_choice", "Round 7.25: points on that grid.",
       "Half-decade steps; each costs two reversal solves", 2, 20),
    _p("wallsearch.depth_step", "Search: Ca2+ well depth step", 2.0, "kT",
       "method", "pore3d", "method_choice", "Round 7.25: the Ca2+-only well "
       "is tried at this step up to the maximum.", "Resolves the rise and "
       "fall of P_Ca:P_K with depth (both within 2-8 kT on 8TKF)", 0.25,
       10.0),
    _p("wallsearch.depth_max", "Search: deepest Ca2+ well", 8.0, "kT",
       "method", "pore3d", "method_choice", "Round 7.25: the deepest "
       "Ca2+-only well tried.", "Past its peak on 8TKF: the well fills with "
       "Ca2+ and its own charge repels further entry", 1.0, 30.0),
    _p("wallsearch.ceiling_depth", "Search: well depth for the ceiling",
       15.0, "kT", "method", "pore3d", "method_choice", "Round 7.25: the "
       "Ca2+-only well's depth when reading how far any well can lift "
       "P_Ca:P_K in linear response.", "e^15 leaves the well's own "
       "resistance below 1e-6 of the rest; 8TKF's span ceiling moves 1.5 % "
       "from 10 to 20 kT, and deeper wells stall the weighted Laplace solve",
       5.0, 25.0),
    _p("wallsearch.ring_step", "Search: ring spacing", 5.0, "A", "method",
       "pore3d", "method_choice", "Round 7.25: opposite-charge ring pairs "
       "are tried at this axial step through the membrane span.", "About "
       "the 3-D charge's Gaussian width (pore_charge.smoothing) times two, "
       "so neighbouring heights overlap", 1.0, 20.0),
    _p("wallsearch.ring_charge_max", "Search: ring charge", 2.0, "e",
       "method", "pore3d", "method_choice", "Round 7.25: each ring is tried "
       "at 1 e per site and at this.", "Two charged groups per subunit at "
       "one height, as the densest lining rings of the deposits", 0.1,
       4.0),
    _p("wallsearch.gui_well_depth", "Drawn Ca2+ well depth", 4.0, "kT",
       "method", "pore3d", "method_choice", "Round 7.26: the depth of the "
       "Ca2+-only well over the span that the lumen box draws at reversal.",
       "Round 7.25: the best uncharged well on both 8TKF (score 1.30) and "
       "7T3T (1.10); 6 kT is within 2 % on 8TKF, 8 kT past the peak", 0.5,
       15.0),
]
