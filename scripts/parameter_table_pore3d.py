"""The voxelised pore and its 3-D ohmic solve (Round 7.6, the conductance
shortfall). Method constants only: the electrolyte (diffusivities, radii,
bath, temperature) is the 1-D model's, so the two readings differ in shape
alone. Imported into ``parameter_table.P``.
"""

from param_entry import entry as _p

PORE3D = [
    _p("pore3d.spacing", "Voxel spacing", 0.5, "A", "method", "pore3d",
       "method_choice", "Edge of the cubic voxels the ion-accessible volume "
       "is cut into.", "8TKF reads 90 / 95 / 98 pS at 1.0 / 0.75 / 0.5 A, "
       "converging upward roughly linearly in the spacing; the CLI reports "
       "the reading at this spacing and the extrapolation from twice it",
       0.1, 2.0),
    _p("pore3d.box_half_width", "Box half-width", 30.0, "A", "method",
       "pore3d", "method_choice", "Half-width of the grid about the pore "
       "axis; beyond the membrane its side faces are bath.", "40 A moves "
       "8TKF by < 0.1 %", 5.0, 100.0),
    _p("pore3d.bath_margin", "Bath beyond the span", 25.0, "A", "method",
       "pore3d", "method_choice", "Grid reach beyond each end of the "
       "pore-domain span, into the baths.", "40 A moves 8TKF by < 0.5 %; the "
       "1-D window (S0's 12 A) is shorter", 2.0, 100.0),
    _p("pore3d.seal_radius", "Membrane seal radius", 25.0, "A", "method",
       "pore3d", "method_choice", "Inside the span only voxels this close "
       "to the axis may conduct: beyond the pore-domain helices the atom "
       "model has empty space where a cell has bilayer.", "20 and 25 A give "
       "8TKF the same conductance; 35 A opens a path through the lipid "
       "space (1.3 nS), which is the failure this radius exists to prevent "
       "(scanned by `shortfall --scan`)", 5.0, 60.0),
    _p("pore3d.cg_tolerance", "Laplace solver tolerance", 1e-9, "",
       "method", "pore3d", "method_choice", "Relative residual at which the "
       "conjugate-gradient solve stops.", "Well below the grid error", 1e-14,
       1e-3),
    _p("pore3d.cg_max_iterations", "Laplace solver iteration cap", 50000,
       "", "method", "pore3d", "method_choice", "Iterations before the solve "
       "is reported unconverged.", "A 0.5 A deposit grid converges in a few "
       "thousand", 100, 1000000),
]
