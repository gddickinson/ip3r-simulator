"""The literature every registered parameter may cite.

Declared here, written to ``ip3r/resources/references.json`` by
``build_parameters.py``. A parameter whose ``citation`` is not a key below (or
one of the sentinels) fails the build.

Identifiers are DOIs, which are checkable by anyone with a browser; nothing
bibliographic here was computed. ``ip3r_genes`` is this application's sister
project, cited when a number is *its* measurement or method so that the two
projects compute the same quantity the same way.
"""

REFERENCES = [
    {"key": "deyoung1992", "authors": "De Young GW, Keizer J", "year": 1992,
     "title": "A single-pool inositol 1,4,5-trisphosphate-receptor-based "
              "model for agonist-stimulated oscillations in Ca2+ concentration",
     "journal": "Proc Natl Acad Sci USA 89:9895-9899",
     "doi": "10.1073/pnas.89.20.9895"},
    {"key": "li1994", "authors": "Li YX, Rinzel J", "year": 1994,
     "title": "Equations for InsP3 receptor-mediated [Ca2+]i oscillations "
              "derived from a detailed kinetic model: a Hodgkin-Huxley like "
              "formalism",
     "journal": "J Theor Biol 166:461-473", "doi": "10.1006/jtbi.1994.1041"},
    {"key": "bezprozvanny1991",
     "authors": "Bezprozvanny I, Watras J, Ehrlich BE", "year": 1991,
     "title": "Bell-shaped calcium-response curves of Ins(1,4,5)P3- and "
              "calcium-gated channels from endoplasmic reticulum of "
              "cerebellum",
     "journal": "Nature 351:751-754", "doi": "10.1038/351751a0"},
    {"key": "mak1998", "authors": "Mak DO, McBride S, Foskett JK", "year": 1998,
     "title": "Inositol 1,4,5-trisphosphate activation of inositol "
              "trisphosphate receptor Ca2+ channel by ligand tuning of Ca2+ "
              "inhibition",
     "journal": "Proc Natl Acad Sci USA 95:15821-15825",
     "doi": "10.1073/pnas.95.26.15821"},
    {"key": "shuai2002", "authors": "Shuai JW, Jung P", "year": 2002,
     "title": "Stochastic properties of Ca2+ release of inositol "
              "1,4,5-trisphosphate receptor clusters",
     "journal": "Biophys J 83:87-97",
     "doi": "10.1016/S0006-3495(02)75151-5"},
    {"key": "swillens1999",
     "authors": "Swillens S, Dupont G, Combettes L, Champeil P", "year": 1999,
     "title": "From calcium blips to calcium puffs: theoretical analysis of "
              "the requirements for interchannel communication",
     "journal": "Proc Natl Acad Sci USA 96:13750-13755",
     "doi": "10.1073/pnas.96.24.13750"},
    {"key": "smith2009", "authors": "Smith IF, Parker I", "year": 2009,
     "title": "Imaging the quantal substructure of single IP3R channel "
              "activity during Ca2+ puffs in intact mammalian cells",
     "journal": "Proc Natl Acad Sci USA 106:6404-6409",
     "doi": "10.1073/pnas.0810799106"},
    {"key": "foskett2007",
     "authors": "Foskett JK, White C, Cheung KH, Mak DO", "year": 2007,
     "title": "Inositol trisphosphate receptor Ca2+ release channels",
     "journal": "Physiol Rev 87:593-658",
     "doi": "10.1152/physrev.00035.2006"},
    {"key": "paknejad2018", "authors": "Paknejad N, Hite RK", "year": 2018,
     "title": "Structural basis for the regulation of inositol "
              "trisphosphate receptors by Ca2+ and IP3",
     "journal": "Nat Struct Mol Biol 25:660-668",
     "doi": "10.1038/s41594-018-0089-6"},
    {"key": "atilgan2001",
     "authors": "Atilgan AR, Durell SR, Jernigan RL, Demirel MC, Keskin O, "
                "Bahar I", "year": 2001,
     "title": "Anisotropy of fluctuation dynamics of proteins with an "
              "elastic network model",
     "journal": "Biophys J 80:505-515",
     "doi": "10.1016/S0006-3495(01)76033-X"},
    {"key": "bruschweiler1995", "authors": "Brüschweiler R", "year": 1995,
     "title": "Collective protein dynamics and nuclear spin relaxation",
     "journal": "J Chem Phys 102:3396-3403", "doi": "10.1063/1.469213"},
    {"key": "tama2001", "authors": "Tama F, Sanejouand YH", "year": 2001,
     "title": "Conformational change of proteins arising from normal mode "
              "calculations",
     "journal": "Protein Eng 14:1-6", "doi": "10.1093/protein/14.1.1"},
    {"key": "yang2009", "authors": "Yang L, Song G, Jernigan RL", "year": 2009,
     "title": "Protein elastic network models and the ranges of "
              "cooperativity",
     "journal": "Proc Natl Acad Sci USA 106:12347-12352",
     "doi": "10.1073/pnas.0902159106"},
    {"key": "kabsch1976", "authors": "Kabsch W", "year": 1976,
     "title": "A solution for the best rotation to relate two sets of "
              "vectors",
     "journal": "Acta Cryst A32:922-923", "doi": "10.1107/S0567739476001873"},
    {"key": "hanley1982", "authors": "Hanley JA, McNeil BJ", "year": 1982,
     "title": "The meaning and use of the area under a receiver operating "
              "characteristic (ROC) curve",
     "journal": "Radiology 143:29-36",
     "doi": "10.1148/radiology.143.1.7063747"},
    {"key": "ip3r_genes", "authors": "ip3r_genes project", "year": 2026,
     "title": "Census, evolution and constraint of the IP3 receptor family "
              "(the publication project this application illustrates)",
     "journal": "github.com/gddickinson/ip3r_genes", "doi": ""},
]
