"""Round 7.8's publication views: Paper 3's lesion strata (S15a/S15b, the
within-genome identity-matched sign test) and the uncertainty drawn on
Paper 5's VUS thresholds. Imported into ``parameter_table.P``.
"""

from param_entry import entry as _p

VIEWS = [
    _p("lesion.coverage_bar", "Lesion scoring coverage bar", 0.70, "",
       "method", "lesion", "ip3r_genes", "A located locus is scored for "
       "frameshifts and internal stops only at or above this coverage of "
       "the bait; the genome's best-covered scored locus stands for the cell.",
       "ip3r_genes S15a COV_FULL (scripts/s15_integrity.py), the sweep's own "
       "COV_FOUND", 0.0, 1.0),
    _p("lesion.identity_window", "Sibling identity window", 0.02, "",
       "method", "lesion", "ip3r_genes", "A sibling locus enters a cell's "
       "within-genome comparison only if its identity to its own bait is "
       "within this of the cell's (D47's identity matching).",
       "ip3r_genes S15a IDENTITY_WINDOW (scripts/s15_integrity.py)", 0.0, 1.0),
    _p("lesion.min_untied", "Minimum untied pairs for a p", 8, "pairs",
       "convention", "lesion", "ip3r_genes", "A class stratum with fewer "
       "untied pairs is reported with its counts and no p-value, and stays "
       "out of the Benjamini-Hochberg family.",
       "ip3r_genes S15b MIN_CLASS_N (scripts/s15b_fossils.py)", 1, 1000),
    _p("vus.median_level", "Threshold interval level", 0.95, "",
       "convention", "variants", "convention", "Coverage of the exact "
       "order-statistic interval drawn around each labelled class's median "
       "(Paper 5 §8's thresholds); a VUS inside an interval is 'near a "
       "median'. Too few positions for this level leave the interval "
       "unbounded.", "Conover, Practical Nonparametric Statistics §3.2: the "
       "distribution-free binomial interval for a median; 0.95 by convention",
       0.5, 0.999),
]
