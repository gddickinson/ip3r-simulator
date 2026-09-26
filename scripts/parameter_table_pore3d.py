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
    _p("lumen.constriction_half_width", "Drop read across a constriction",
       3.0, "A", "method", "pore3d", "method_choice", "Half-width of the "
       "z interval about a constriction over which the share of the "
       "window's voltage drop is read (Round 7.10).", "About one ion "
       "diameter either side; the 1-D and 3-D shares are read over the same "
       "interval, so the comparison does not depend on it", 0.5, 15.0),
    _p("display.lumen_radius", "Lumen drawn within", 15.0, "A",
       "convention", "display", "convention", "The lumen surface is drawn "
       "for conducting voxels this close to the axis, inside S0's window.",
       "A view choice: beyond it the cytosolic side opens into the space "
       "between domains, which would hide the pore; the solve itself uses "
       "the whole box", 3.0, 30.0),
    _p("display.lumen_smoothing", "Lumen surface smoothing", 0.7, "voxel",
       "convention", "display", "convention", "Gaussian width applied to "
       "the voxel mask before the surface is contoured at one half.",
       "Removes the voxel staircase. Drawing only: a one-voxel neck may "
       "look closed, but every number is read from the unsmoothed mask", 0.0, 2.0),
    _p("display.lumen_alpha", "Lumen surface opacity", 0.6, "",
       "convention", "display", "convention", "Opacity of the drawn lumen.",
       "Lets the lining side chains show through", 0.1, 1.0),
    _p("display.lumen_potential_range", "Lumen wall-potential scale", 5.0,
       "kT/e", "convention", "display", "convention", "The wall potential "
       "is drawn on a fixed diverging ramp from minus to plus this value "
       "(Round 7.12).", "Fixed, never auto-ranged: 8TKF's and 9HEO's "
       "readings under every closure lie within -4.9 to +3.9 kT/e, so no "
       "reading saturates and deposits compare by colour", 0.5, 20.0),
    _p("charge3d.gaussian_reach", "3-D charge Gaussian reach", 3.0,
       "width", "method", "pore3d", "method_choice", "A lining group's "
       "charge is spread over the lumen voxels within this many Gaussian "
       "widths (pore_charge.smoothing) of its centre, normalised over them "
       "(Round 7.11, the local and pb closures).", "Beyond 3 widths the "
       "Gaussian holds 1 % of its weight; the charge is conserved at any "
       "reach, only its shape is cut", 1.0, 6.0),
    _p("charge3d.newton_tolerance", "Poisson-Boltzmann Newton tolerance",
       1e-6, "kT/e", "method", "pore3d", "method_choice", "Largest "
       "potential update at which the 3-D Poisson-Boltzmann Newton "
       "iteration stops.", "A conductance moves as e^{-u}, so 1e-6 kT is "
       "a 1e-6 relative change", 1e-12, 1e-2),
    _p("charge3d.newton_max_iterations", "Poisson-Boltzmann Newton cap", 60,
       "", "method", "pore3d", "method_choice", "Newton iterations before "
       "the 3-D Poisson-Boltzmann solve is reported unconverged.",
       "With the step cap, 8TKF converges well inside it", 5, 1000),
    _p("charge3d.max_step", "Poisson-Boltzmann step cap", 2.0, "kT/e",
       "method", "pore3d", "method_choice", "Largest potential change one "
       "Newton step may make; a larger step is scaled down.", "The "
       "exponential nonlinearity overshoots from a far start; any cap of "
       "order kT gives the same root", 0.1, 20.0),
    _p("dielectric.eps_protein", "Protein relative permittivity", 4.0, "",
       "method", "pore3d", "schutz2001", "Permittivity of every voxel not "
       "open to the ion probe (protein, and the sealed membrane) in the "
       "dielectric closure (Round 7.13).", "The value for charge-charge "
       "interactions inside a protein depends on what the model treats "
       "explicitly; 2-4 for a model without relaxation, higher for one "
       "that folds it in. Swept 2-20, never tuned", 1.0, 40.0),
    _p("dielectric.charge_width", "Dielectric charge width", 1.0, "A",
       "method", "pore3d", "method_choice", "Gaussian width of each "
       "ionisable group's charge about its own centre in the dielectric "
       "closure, over every voxel it reaches (Round 7.13).", "About the "
       "spread of a carboxylate's or guanidinium's charge over its atoms; "
       "it regularises the point charge on the grid, and is swept", 0.25,
       6.0),
]
