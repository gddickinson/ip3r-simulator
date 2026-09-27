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
]
