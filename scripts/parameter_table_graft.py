"""Filling unresolved stretches from an AlphaFold model: the anchor, seam
and clash rules and AlphaFold's own confidence bands. Imported into
``parameter_table.P``.
"""

from param_entry import entry as _p

_BANDS = ("AlphaFold DB's four pLDDT bands (very low < 50 <= low < 70 <= "
          "confident < 90 <= very high), as defined in Jumper 2021 and used "
          "on every AlphaFold DB entry page")

GRAFT = [
    _p("graft.anchor_window", "Anchor window", 6, "residues", "method",
       "graft", "method_choice", "Residues searched on each side of an "
       "unresolved stretch for anchors: resolved in the deposit, not a "
       "backbone+CB stub, and the same amino acid in the prediction.",
       "Short enough that the fit is local to the seam (PIEZO1's blade fit "
       "went from 19 A globally to 2.4 A locally), long enough to hold "
       "two helical turns", 2, 50),
    _p("graft.min_anchor", "Minimum anchors per side", 3, "residues",
       "method", "graft", "method_choice", "An internal stretch is filled "
       "only with at least this many anchors on each side; a terminus, "
       "which has one side, needs twice as many there.",
       "Three points fix a rigid-body placement; fewer on one side lets "
       "the fill pivot about the other", 1, 20),
    _p("graft.join_tolerance", "Seam join tolerance", 5.5, "A", "method",
       "graft", "measured_here", "A seam joins if the placed fill's end "
       "C-alpha lies within this distance of the deposit's flanking "
       "C-alpha; otherwise the seam is drawn and reported as broken.",
       "A trans peptide spans 3.80 A. Hiding 31 stretches that 8TKG/8TKF "
       "resolve and filling them (structure/graft_calibration.py) gave seams "
       "up to 5.42 A (90th percentile 4.74) for fills of median 1.5 A RMSD; "
       "the first value, 4.5 A, failed 9 of those 62 true seams", 3.0, 20.0),
    _p("graft.align_min_identity", "Alignment fill identity", 0.99, "",
       "method", "graft", "measured_here", "A model not in the deposit's "
       "numbering may fill it through an alignment of the deposit's whole "
       "construct only if the two agree over the aligned pairs at least "
       "this well: the same protein, another isoform.",
       "7LHF's construct aligned to each downloaded model "
       "(structure/graft_numbering.py): rat ITPR1 isoform 8 1.000 (the "
       "same protein; only the SI and SII splice segments unpaired), human "
       "ITPR1 isoform 4 0.988 (the ortholog, ~32 substitutions), ITPR3 "
       "0.658; 9YKK (ITPR2) reaches 0.713 at best. The bar admits only the "
       "same protein's own isoform", 0.5, 1.0),
    _p("graft.clash_distance", "Clash distance", 2.2, "A", "method",
       "graft", "method_choice", "A filled residue clashes if any of its "
       "heavy atoms is closer than this to a deposited heavy atom (the "
       "residues flanking its own seams excluded).",
       "Well inside any non-bonded heavy-atom contact (>= ~3 A); only "
       "interpenetration counts", 1.0, 4.0),
    _p("graft.plddt_confident", "pLDDT confident", 70.0, "", "convention",
       "graft", "jumper2021", "A predicted residue is 'confident' at or "
       "above this pLDDT; a fill's confident fraction is reported.",
       _BANDS, 0.0, 100.0),
    _p("display.plddt_low", "pLDDT low / very low edge", 50.0, "",
       "convention", "graft", "varadi2022", "pLDDT below this is drawn "
       "in the very-low colour.", _BANDS, 0.0, 100.0),
    _p("display.plddt_very_high", "pLDDT very-high edge", 90.0, "",
       "convention", "graft", "varadi2022", "pLDDT at or above this is "
       "drawn in the very-high colour.", _BANDS, 0.0, 100.0),
]
