"""IP3R dynamic structural simulator.

Subpackages, in dependency order (``io -> core -> structure -> physics ->
analysis``; ``render`` and ``ui`` consume all of them and nothing below them
imports either):

- ``io``        mmCIF reading, the structure registry, downloads
- ``core``      the Structure container, paralog annotation, ip3r_genes tables
- ``structure`` four-fold symmetry, pore profile, ligand contacts, numbering
- ``physics``   elastic-network modes, IP3R gating, cell Ca2+, puffs
- ``analysis``  the findings checks that re-derive ip3r_genes results
- ``render``    moderngl renderer (impostor spheres/cylinders, cartoons)
- ``ui``        the PyQt6 application

See INTERFACE.md for the file-by-file map.
"""

__version__ = "0.1.0"
