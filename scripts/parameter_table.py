"""Every number a calculation in ``ip3r`` depends on, with its provenance.

Validated and written to ``ip3r/resources/parameters.json`` by
``build_parameters.py``. Consumed at call time as ``_P.value("key")``.

Kinds: ``physical`` (a property of the world), ``empirical`` (a fitted model
constant), ``method`` (an algorithmic choice) and ``convention``.

The gating constants are the De Young-Keizer set as reduced by Li & Rinzel.
Dissociation constants are ``d_i = b_i / a_i``; the reduced model needs
``d1, d2, d3, d5`` and ``a2``, and the stochastic subunit model also needs the
on-rates ``a1`` and ``a5`` of the fast IP3 and activating-Ca2+ sites.
"""

def _p(key, name, value, unit, kind, category, citation, description,
       source_note="", minimum=None, maximum=None):
    return {"key": key, "name": name, "value": value, "unit": unit,
            "kind": kind, "category": category, "citation": citation,
            "description": description, "source_note": source_note,
            "minimum": minimum, "maximum": maximum}


_DYK = ("De Young & Keizer 1992, Table 1; the same values are used "
        "unchanged by Li & Rinzel 1994")

P = [
    # ---------------------------------------------------------- gating (DYK)
    _p("gating.a1", "IP3 on-rate a1", 400.0, "1/(uM s)", "empirical",
       "gating", "deyoung1992", "Binding rate of IP3 to its site on one "
       "subunit.", _DYK, 1.0, 1e4),
    _p("gating.a2", "Inhibitory Ca2+ on-rate a2", 0.2, "1/(uM s)", "empirical",
       "gating", "deyoung1992", "Binding rate of Ca2+ to the slow inhibitory "
       "site; sets the time scale of inactivation (Li-Rinzel h gate).",
       _DYK, 0.01, 10.0),
    _p("gating.a5", "Activating Ca2+ on-rate a5", 20.0, "1/(uM s)",
       "empirical", "gating", "deyoung1992", "Binding rate of Ca2+ to the "
       "fast activating site.", _DYK, 0.1, 1e3),
    _p("gating.d1", "IP3 dissociation constant d1", 0.13, "uM", "empirical",
       "gating", "deyoung1992", "IP3 site affinity when the inhibitory site "
       "is empty.", _DYK, 0.001, 10.0),
    _p("gating.d2", "Inhibitory Ca2+ dissociation constant d2", 1.049, "uM",
       "empirical", "gating", "deyoung1992", "Inhibitory-site affinity with "
       "IP3 bound.", _DYK, 0.01, 100.0),
    _p("gating.d3", "IP3 dissociation constant d3", 0.9434, "uM", "empirical",
       "gating", "deyoung1992", "IP3 site affinity when the inhibitory site "
       "is occupied.", _DYK, 0.001, 100.0),
    _p("gating.d5", "Activating Ca2+ dissociation constant d5", 0.08234, "uM",
       "empirical", "gating", "deyoung1992", "Activating-site affinity.",
       _DYK, 0.001, 10.0),
    _p("gating.subunits_required", "Active subunits to open", 3.0, "",
       "empirical", "gating", "li1994", "A channel conducts when this many "
       "of its four subunits are in the active state. The deterministic "
       "Li-Rinzel flux uses the cube of the subunit activity, i.e. three.",
       "Li & Rinzel 1994 (the cubic exponent); Shuai & Jung 2002 use the same "
       "3-of-4 rule for single stochastic channels", 1.0, 4.0),

    # ------------------------------------------------------- cell (Li-Rinzel)
    _p("cell.v1", "Maximal channel flux rate v1", 6.0, "1/s", "empirical",
       "cell", "deyoung1992", "Rate constant of IP3R-mediated release.",
       _DYK, 0.0, 100.0),
    _p("cell.v2", "ER leak rate v2", 0.11, "1/s", "empirical", "cell",
       "deyoung1992", "Passive leak from the ER.", _DYK, 0.0, 10.0),
    _p("cell.v3", "SERCA maximal rate v3", 0.9, "uM/s", "empirical", "cell",
       "deyoung1992", "Maximal SERCA uptake.", _DYK, 0.0, 100.0),
    _p("cell.k3", "SERCA half-activation k3", 0.1, "uM", "empirical", "cell",
       "deyoung1992", "Ca2+ at half-maximal SERCA uptake (Hill 2).", _DYK,
       0.001, 10.0),
    _p("cell.c0", "Total cell Ca2+ c0", 2.0, "uM", "empirical", "cell",
       "deyoung1992", "Total free Ca2+ referred to the cytosolic volume; "
       "conserved in the closed-cell model.", _DYK, 0.1, 100.0),
    _p("cell.c1", "ER/cytosol volume ratio c1", 0.185, "", "empirical",
       "cell", "deyoung1992", "Ratio of ER to cytosolic volume.", _DYK,
       0.01, 10.0),

    # -------------------------------------------------------- puffs (method)
    _p("puff.n_channels", "Channels per cluster", 20.0, "", "method", "puff",
       "method_choice", "Number of IP3Rs in the simulated cluster.",
       "Order-of-magnitude choice: Smith & Parker 2009 resolve puffs made of "
       "a few to a few tens of channel openings; the cluster modelled here "
       "is a teaching default, not a measurement", 1.0, 200.0),
    _p("puff.ca_rest", "Resting cytosolic Ca2+", 0.1, "uM", "convention",
       "puff", "convention", "Background Ca2+ each channel sees when no "
       "channel in the cluster is open.",
       "The customary 100 nM resting level; the same value is the lower "
       "bound of the Li-Rinzel oscillation in its own figures", 0.0, 10.0),
    _p("puff.ca_per_open", "Cluster Ca2+ per open channel", 1.0, "uM",
       "method", "puff", "method_choice", "Mean-field increment of the Ca2+ "
       "every channel in the cluster sees for each open channel — the "
       "coupling that turns a blip into a puff.",
       "A mean-field simplification of Swillens et al. 1999, who show that "
       "puffs need inter-channel Ca2+ coupling; 1 uM per channel was chosen "
       "from a scan (0-2 uM at 0.05-0.2 uM IP3) as the value that most raises "
       "the Fano factor of simultaneous openings. The local nanodomain at an "
       "open pore is far higher and is not modelled", 0.0,
       20.0),
    _p("puff.dt", "Puff simulation time step", 0.2e-3, "s", "method", "puff",
       "method_choice", "Fixed step of the stochastic gate simulation.",
       "Chosen so the fastest transition probability per step (a1 x IP3 or "
       "a5 x cluster Ca2+) stays below about 0.1 at the default parameters",
       1e-6, 1e-2),

    # ------------------------------------------------------------ ANM
    _p("anm.cutoff", "ANM contact cutoff", 15.0, "A", "method", "anm",
       "atilgan2001", "C-alpha pairs closer than this are joined by a "
       "spring.", "Atilgan et al. 2001 use 13-15 A; the upper end keeps a "
       "coarse-grained tetramer connected", 6.0, 30.0),
    _p("anm.gamma", "ANM spring constant", 1.0, "arb.", "convention", "anm",
       "convention", "Uniform stiffness scale; only eigenvalue ratios are "
       "reported.", "Eigenvalues are arbitrary without a B-factor fit", 1e-3,
       1e3),
    _p("anm.d0", "ANM distance scale d0", 7.5, "A", "method", "anm",
       "yang2009", "Distance scale of the inverse-square spring weighting.",
       "Yang et al. 2009 parameter-free ANM weighting", 1.0, 20.0),
    _p("anm.n_modes", "Modes computed", 20.0, "", "method", "anm",
       "method_choice", "Non-trivial normal modes kept.",
       "Enough to cover the collective motions; more only costs time", 1.0,
       200.0),
    _p("anm.stride", "C-alpha stride", 3.0, "", "method", "anm",
       "method_choice", "Keep every n-th C-alpha of each subunit.",
       "A full tetramer is ~9,000 C-alphas; a stride of 3 keeps the "
       "shift-invert solve interactive while the low modes, which are "
       "collective by construction, are insensitive to it", 1.0, 10.0),
    _p("anm.symmetry_tolerance", "C4 character tolerance", 0.15, "",
       "method", "anm", "method_choice", "How close to an ideal character "
       "(+1, -1, 0) a mode must be to receive an irrep label.",
       "Coarse-grained deposits are not perfectly symmetric; beyond this "
       "the mode is reported as mixed rather than forced into a label",
       0.01, 0.5),

    _p("anm.min_collectivity", "Collective-mode threshold", 0.2, "",
       "method", "anm", "bruschweiler1995", "Modes whose collectivity "
       "(Bruschweiler's kappa, the fraction of sites effectively moving) is "
       "below this are reported as local network artefacts, not collective "
       "motions.",
       "The measure is Bruschweiler 1995; the threshold is ours. On the 8TKG "
       "network at strides 1-5 every mode is either <= 0.11 (weakly attached "
       "fragments, near-zero eigenvalues) or >= 0.27; 0.2 sits in that gap. "
       "At stride 4 the 'lowest A mode' was such a fragment (kappa 0.01, "
       "overlap 0.008) while the collective one below it had 0.49", 0.0,
       1.0),

    # ------------------------------------------------- transition and morph
    _p("transition.min_residues", "Minimum common residues", 50.0, "",
       "method", "transition", "method_choice", "Fewer residues common to "
       "both deposits than this and no transition is built.",
       "A sanity floor, not a tuning knob: two deposits of one paralog share "
       "thousands (8TKG/8TKF: 2,198 per subunit); fewer than 50 means the "
       "wrong pair or the wrong numbering", 4.0, 1000.0),
    _p("transition.tie_tolerance", "Subunit correspondence tie", 0.05, "A",
       "method", "transition", "method_choice", "Cyclic subunit "
       "correspondences whose fit RMSD is within this of the best are "
       "treated as equally good; the tie goes to the deposited chain labels.",
       "A C4-symmetric pair fits every cyclic relabelling to within "
       "~0.01 A (8TKG->8TKF: 14.18 A for all four); 0.05 A is above that "
       "and far below any real mismatch", 0.0, 1.0),
    _p("morph.n_frames", "Morph frames", 30.0, "", "method", "morph",
       "method_choice", "Solved frames between the two endpoints.",
       "As the PIEZO1 simulator: smooth playback; more only costs time", 2.0,
       200.0),
    _p("morph.iterations", "Morph restraint iterations", 60.0, "", "method",
       "morph", "method_choice", "SHAKE-style passes restoring peptide "
       "C-alpha-C-alpha distances per frame.",
       "As the PIEZO1 simulator; the reported bond error shows whether it "
       "was enough", 0.0, 1000.0),
    _p("morph.max_bond", "Peptide bond cut", 6.0, "A", "method", "morph",
       "method_choice", "Consecutive basis C-alphas farther apart than this "
       "in either endpoint are a chain break, not a bond, and are not "
       "restrained.",
       "A trans peptide C-alpha-C-alpha distance is 3.8 A; 6 A admits "
       "coordinate error and excludes any gap of one missing residue or "
       "more", 4.0, 10.0),
    _p("display.displacement_max", "Displacement colour scale top", 25.0,
       "A", "convention", "display", "convention", "Per-residue displacement "
       "at which the fixed displacement ramp saturates.",
       "Fixed, never auto-ranged, so two transitions painted with it can be "
       "compared; 8TKG->8TKF fitted on the pore domain has a 95th "
       "percentile of 24.1 A and a maximum of 30.6 A, so only the most "
       "mobile few per cent of the cytosolic assembly saturate", 1.0, 100.0),

    # ------------------------------------------------------ structure geometry
    _p("pore.step", "Pore profile step", 0.5, "A", "method", "pore",
       "ip3r_genes", "Axial sampling interval of the pore profile.",
       "The value ip3r_genes S0 used, so the two profiles share a grid", 0.1,
       5.0),
    _p("pore.slab", "Pore profile half-slab", 1.5, "A", "method", "pore",
       "ip3r_genes", "Half-thickness of the slab in which the minimum "
       "heavy-atom distance to the axis is taken.",
       "The value ip3r_genes S0 used for structure_pore.tsv", 0.25, 5.0),
    _p("pore.lining_slab", "Lining-residue half-slab", 2.0, "A", "method",
       "pore", "ip3r_genes", "Half-thickness of the slab searched for the "
       "residues forming a constriction.", "As ip3r_genes S0", 0.5, 6.0),
    _p("pore.lining_tol", "Lining-residue radial tolerance", 1.2, "A",
       "method", "pore", "ip3r_genes", "Atoms within this distance of the "
       "minimum radius count as lining the constriction.", "As ip3r_genes S0",
       0.1, 5.0),
    _p("ligand.contact_cutoff", "Ligand contact cutoff", 4.5, "A",
       "convention", "ligand", "ip3r_genes", "Heavy-atom distance at which a "
       "residue counts as an IP3 contact.",
       "The cutoff ip3r_genes S0 and S22 used for the ten measured contacts",
       2.5, 8.0),
    _p("ligand.module_min_coverage", "Module coverage floor", 0.5, "",
       "method", "ligand", "ip3r_genes", "Fraction of a module's reference "
       "columns an orthologue must resolve for its identity to enter the "
       "paired core-vs-pore test.", "S22's MIN_MODULE_COVERAGE; a tip "
       "truncated at one end would otherwise enter as an extreme divergence "
       "in that module", 0.0, 1.0),
    _p("ligand.shell_second_edge", "Second-shell outer edge", 8.0, "A",
       "convention", "ligand", "ip3r_genes", "All-atom distance to IP3 "
       "separating the second ligand shell from the third (the first shell "
       "is ligand.contact_cutoff).", "S22's SHELL_EDGES: equal 3.5 A steps "
       "out from the 4.5 A contact cutoff", 4.5, 15.0),
    _p("ligand.shell_third_edge", "Third-shell outer edge", 11.5, "A",
       "convention", "ligand", "ip3r_genes", "All-atom distance to IP3 "
       "separating the third ligand shell from the fourth.",
       "S22's SHELL_EDGES", 4.5, 15.0),
    _p("ligand.shell_radius", "Ligand-shell search radius", 15.0, "A",
       "convention", "ligand", "ip3r_genes", "Residues farther than this "
       "from IP3 belong to no shell (absent, not an open last bin).",
       "S22's SEARCH_RADIUS_A: past 15 A no side chain of the binding core "
       "reaches the ligand", 8.0, 30.0),
    _p("align.gap_open", "Gap-open cost", 10.0, "BLOSUM62 units",
       "method", "structure", "rice2000", "Cost of opening a gap in the "
       "pairwise alignment that carries residue numbers between paralogs.",
       "EMBOSS needle's default with BLOSUM62 [henikoff1992]; the recursion "
       "is Gotoh's [gotoh1982]", 1.0, 30.0),
    _p("align.gap_extend", "Gap-extend cost", 0.5, "BLOSUM62 units",
       "method", "structure", "rice2000", "Cost of each further residue of "
       "a gap in the paralog-transfer alignment.",
       "EMBOSS needle's default", 0.0, 5.0),
    _p("numbering.min_identity", "Numbering-check identity", 0.95, "",
       "method", "structure", "ip3r_genes", "Fraction of residues a "
       "structure must share, by number, with a reference sequence to be "
       "declared in that reference's numbering.",
       "The rule ip3r_genes S24 (D64) applied before painting variants on "
       "structures; strict, because an offset still aligns most residues "
       "with something", 0.5, 1.0),

    _p("constraint.min_occupancy", "Minimum column occupancy", 0.5, "",
       "method", "constraint", "ip3r_genes", "A conservation score is used "
       "only where at least this fraction of (weighted) sequences have a "
       "residue in the column.", "S17's MIN_OCCUPANCY, so re-derived AUCs are "
       "computed on the positions S17 scored", 0.0, 1.0),

    # ------------------------------------------------------ check tolerances
    _p("check.length_tol", "Geometry agreement tolerance", 0.05, "A",
       "method", "checks", "method_choice", "Largest difference in a "
       "radius, residual or distance that still counts as reproducing the "
       "published value.", "The published values are quoted to 0.01 A; "
       "0.05 A allows for rounding and axis-fit differences and nothing "
       "else", 0.0, 1.0),
    _p("check.alpha", "Significance level", 0.05, "",
       "convention", "checks", "convention", "Level at which a check reads "
       "a published or re-derived p (or q) as significant when testing the "
       "pattern the prose states.", "The level ip3r_genes' prose uses "
       "throughout", 0.0, 1.0),
    _p("check.log_p_tol", "P-value agreement tolerance", 0.02, "decades",
       "method", "checks", "method_choice", "Largest difference in log10 "
       "of a re-derived and a published p-value that still counts as "
       "agreement.", "Published p-values carry six significant figures, "
       "and a normal approximation written independently of scipy's agrees "
       "with it to well under 0.01 decades; 0.02 leaves room for rounding "
       "only", 0.0, 1.0),
    _p("check.stat_tol", "Statistic agreement tolerance", 0.002, "",
       "method", "checks", "method_choice", "Largest difference in a "
       "re-derived AUC, mean or fraction that still counts as agreement.",
       "Published tables round to four decimals; 0.002 allows for rounding "
       "and tie-handling differences", 0.0, 0.1),
    # ------------------------------------------------------ tree support
    _p("tree.alrt_min", "SH-aLRT support bar", 80.0, "%", "convention",
       "phylogeny", "guindon2010", "A tree node counts as supported only if "
       "its SH-like approximate likelihood-ratio support is at least this.",
       "The IQ-TREE recommendation for SH-aLRT, and half of the joint bar "
       "ip3r_genes Paper 2 holds every claim to", 0.0, 100.0),
    _p("tree.ufboot_min", "UFBoot support bar", 95.0, "%", "convention",
       "phylogeny", "hoang2018", "A tree node counts as supported only if "
       "its ultrafast bootstrap support is at least this.",
       "Hoang et al.'s threshold for a clade to be considered real, and the "
       "other half of Paper 2's joint bar", 0.0, 100.0),
    # ------------------------------------------------------ genome grid
    _p("genomes.contiguity_bar_bp", "Contiguity bar (contig N50)", 142212.0,
       "bp", "method", "genomes", "ip3r_genes", "An assembly counts as able "
       "to hold the gene when its contig N50 is at least this; the grid "
       "splits genomes above and below it.", "ip3r_genes D4's bar, fixed "
       "before any error rate was known: the median measured genomic span "
       "of an ITPR gene (Paper 3, results/methods/contiguity_summary.json "
       "d4_bar)", 1.0, 1e9),
    # ------------------------------------------------------ Paper 1 range
    _p("range.absence_min_proteomes", "Absence-clade proteome floor", 10.0,
       "proteomes", "method", "range", "ip3r_genes", "A phylum or class is "
       "an absence target when at least this many of its reference proteomes "
       "were swept and none carries a call.", "ip3r_genes S23 manifest rule "
       "G3 (results/s23_scope/report.md: 'one per clade S20 swept >= 10 "
       "proteomes of and found 0 ITPR in')", 1.0, 1e4),
    _p("range.substantial_evalue", "Substantial-match E-value", 1e-5, "",
       "method", "range", "ip3r_genes", "A relaxed-search hit is substantial "
       "only at or below this full-sequence E-value (with the coverage bar).",
       "ip3r_genes S20's primary threshold, the strict side of its "
       "one-search-two-sensitivities design (results/s20_sweep/report.md)",
       0.0, 10.0),
    _p("range.substantial_coverage", "Substantial-match model coverage", 0.5,
       "", "method", "range", "ip3r_genes", "A relaxed-search hit is "
       "substantial only if it spans at least this fraction of the profile.",
       "ip3r_genes S20: 'spanning at least half the model' "
       "(results/s20_sweep/report.md, negative-claims table)", 0.0, 1.0),
    _p("range.family_floor_aa", "Full-length family floor", 2000.0, "aa",
       "method", "range", "ip3r_genes", "A plant or fungal record shorter "
       "than this is a fragment of a gene, not a full-length receptor.",
       "ip3r_genes S20 chase threshold min_length_aa "
       "(results/s20_sweep/verdict_summary.json)", 1.0, 1e4),
    _p("range.contaminant_pident", "Contaminant identity", 95.0, "%",
       "method", "range", "ip3r_genes", "A record this identical (or more) "
       "to one relative outside its kingdom, over the coverage bar, is an "
       "assembly-contamination suspect.", "ip3r_genes S20 chase threshold "
       "contaminant_pident (results/s20_sweep/verdict_summary.json)",
       0.0, 100.0),
    _p("range.contaminant_qcov", "Contaminant coverage", 0.5, "",
       "method", "range", "ip3r_genes", "Query coverage the contaminant "
       "identity must hold over.", "ip3r_genes S20 chase threshold "
       "contaminant_qcov (results/s20_sweep/verdict_summary.json)", 0.0, 1.0),
]
