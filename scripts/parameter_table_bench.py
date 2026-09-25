"""Paper 1's family-call benchmark (S1; Round 7.7): the discovery scorer's
component points, its promotion threshold and the labelled-bait margin.
They are ip3r_genes' own constants (``src/discovery/candidates.py``), used
here to rederive each control's score from the components that fired.
Imported into ``parameter_table.P``.
"""

from param_entry import entry as _p


def _points(key, name, value, what):
    return _p(f"bench.points_{key}", f"Score: {name}", value, "points",
              "convention", "benchmark", "ip3r_genes",
              f"Points the S1 discovery scorer adds when {what}.",
              "src/discovery/candidates.py, the additive score S1 benchmarked",
              0, 100)


BENCH = [
    _points("outlier", "twilight-zone outlier", 20,
            "identity to the nearest known paralog sits in the twilight zone"),
    _points("fold", "structural fold", 20, "a Foldseek hit links it to a known "
            "paralog at TM-score >= 0.5"),
    _points("pfam", "Pfam signature", 20, "it carries a family-defining Pfam"),
    _points("signature", "MSA-signature fallback", 15, "it lacks InterPro "
            "records but hits the MSA-derived family signature blocks"),
    _points("size", "size band", 15, "its length is in the family's band"),
    _points("breadth", "taxonomic breadth", 15, "a sibling in another species "
            "is corroborated"),
    _points("cluster", "cluster exclusion", 10, "it does not co-cluster with "
            "a known paralog"),
    _points("homology", "homology provenance", 10, "a family-linked homology "
            "search surfaced it"),
    _points("split", "split annotation", 10, "it sits beside another family "
            "locus (a fragmented gene model)"),
    _p("bench.score_cap", "Score ceiling", 100, "points", "convention",
       "benchmark", "ip3r_genes", "The additive score is capped here.",
       "src/discovery/candidates.py", 1, 1000),
    _p("bench.promotion", "Promotion threshold", 40, "points", "convention",
       "benchmark", "ip3r_genes", "A candidate is promoted (called a family "
       "member) at this score or above; a sister-family call or a missing "
       "family-specific component caps it one point below.",
       "S1's threshold, 'worth manual review'", 1, 100),
    _p("bench.sister_margin", "Labelled-bait margin", 0.10, "", "convention",
       "benchmark", "ip3r_genes", "A candidate is assigned to the family whose "
       "labelled baits it is closer to only when it wins by more than this "
       "identity; closer to the ryanodine receptors by more is a sister call.",
       "Roadmap decision D7 ('two baits within ~10 % is no call'), applied "
       "by S1's sister test", 0.0, 1.0),
]
