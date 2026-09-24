"""The junctional cleft: steady Ca2+ from each open RyR1, at every other one.

Stern, Pizarro & Rios 1997 (J Gen Physiol 110:415) computed their spark
couplings this way, and the geometry here is theirs (Table I, Fig. 7 B,
Appendix):

* A rectangular cleft ``spark.cleft_width`` (60 nm) wide and
  ``spark.cleft_height`` (15 nm) high, as long as the couplon. Diffusion is
  two-dimensional: nothing varies across the height.
* Two rows of release sites, ``spark.channel_spacing`` (30 nm) apart in
  both directions, V (voltage-coupled) and C (Ca2+-gated) channels
  alternating like a chessboard. Only the C channels are simulated; the V
  sites are there because they fix where the C channels sit.
* An open channel releases ``spark.unitary_current`` uniformly over a disc
  ``spark.source_diameter`` (30 nm) across: "the foot process".
* Ca2+ leaves across every edge at a rate proportional to the edge
  concentration, ``D_inf = 2 pi D / (h ln(R_max / h))`` (their Eq. 13): the
  edge treated as a line source that reaches background at
  ``spark.r_max``.

Stern et al. used a Galerkin finite-element solve. This one is a finite-
volume grid of ``spark.cleft_grid`` cells, so conservation is exact: the
flux in equals the flux out across the edges (tested).

**Steady state.** Their Fig. 10 shows the cleft settles in microseconds, so
the Ca2+ at each channel is taken as the superposition of every open
channel's steady field, recomputed only when a channel opens or closes.
Fixed, fast buffers do not change a steady field. A mobile buffer (their
fura-2 runs) does, and needs time-dependent diffusion; it is not here.

**What a channel sees.** The Ca2+ at its own centre, like their Fig. 9
profiles "along the line of centers". That includes its own release while
it is open (the diagonal of the coupling matrix). The authors note that
where the sensor really sits is unknown, so the self-coupling is the least
certain number in the model.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from ..parameters import PARAMETERS as _P

__all__ = ["CleftGeometry", "Couplon", "couplon", "edge_transfer",
           "coupling_matrix", "v_coupling_matrix", "field", "nearest_coupling"]

F_FARADAY = 96485.33212          # C/mol


@dataclass(frozen=True)
class CleftGeometry:
    n_c: int                     # Ca2+-gated channels (half the sites)
    width: float                 # nm
    height: float                # nm
    spacing: float               # nm
    source_diameter: float       # nm
    r_max: float                 # nm
    grid: float                  # nm
    current: float               # pA per open channel
    d_ca: float                  # cm^2/s

    @classmethod
    def from_parameters(cls, n_c: int | None = None) -> "CleftGeometry":
        n = int(round(_P.value("spark.n_channels"))) if n_c is None else int(n_c)
        return cls(n, _P.value("spark.cleft_width"), _P.value("spark.cleft_height"),
                   _P.value("spark.channel_spacing"), _P.value("spark.source_diameter"),
                   _P.value("spark.r_max"), _P.value("spark.cleft_grid"),
                   _P.value("spark.unitary_current"), _P.value("spark.d_ca"))

    @property
    def per_row(self) -> int:
        """Sites per row: 2 n_c sites over two rows."""
        return self.n_c

    @property
    def length(self) -> float:
        return self.per_row * self.spacing


@dataclass(frozen=True)
class Couplon:
    sites: np.ndarray            # (2 n_c, 2) site centres, nm
    is_c: np.ndarray             # (2 n_c,) True for a Ca2+-gated channel

    @property
    def c_sites(self) -> np.ndarray:
        return self.sites[self.is_c]

    @property
    def v_sites(self) -> np.ndarray:
        return self.sites[~self.is_c]


def couplon(g: CleftGeometry) -> Couplon:
    """Two rows, one spacing apart, V and C alternating in a chessboard.

    Row centres sit half a spacing in from each long edge only when the
    width is two spacings (Stern's 60 nm); otherwise the rows are centred.
    """
    xs = (np.arange(g.per_row) + 0.5) * g.spacing
    y0 = 0.5 * (g.width - g.spacing)
    sites, is_c = [], []
    for row in (0, 1):
        for k, x in enumerate(xs):
            sites.append((x, y0 + row * g.spacing))
            is_c.append((k + row) % 2 == 1)
    return Couplon(np.array(sites), np.array(is_c))


def edge_transfer(g: CleftGeometry) -> float:
    """D_inf / D in 1/nm: the edge flux per unit area of edge face is
    ``D_inf c`` (Stern 1997 Eq. 13)."""
    return 2.0 * np.pi / (g.height * np.log(g.r_max / g.height))


def _grid(g: CleftGeometry):
    nx = int(round(g.length / g.grid))
    ny = int(round(g.width / g.grid))
    return nx, ny, g.length / nx, g.width / ny


def _operator(g: CleftGeometry):
    """Finite-volume operator for sum of face fluxes (units D h, nm)."""
    nx, ny, dx, dy = _grid(g)
    kappa = edge_transfer(g)
    idx = np.arange(nx * ny).reshape(nx, ny)
    rows, cols, vals = [], [], []
    diag = np.zeros(nx * ny)
    # Interior faces, x then y; conductance per unit D h = face length / gap.
    for a, b, cond in ((idx[:-1, :], idx[1:, :], dy / dx), (idx[:, :-1], idx[:, 1:], dx / dy)):
        a, b = a.ravel(), b.ravel()
        rows += [a, b]
        cols += [b, a]
        vals += [np.full(a.size, -cond)] * 2
        np.add.at(diag, a, cond)
        np.add.at(diag, b, cond)
    # Edge faces: half a cell to the face, then the Robin transfer, in series.
    for cells, face, gap in ((idx[0, :], dy, dx), (idx[-1, :], dy, dx),
                             (idx[:, 0], dx, dy), (idx[:, -1], dx, dy)):
        np.add.at(diag, cells.ravel(), face / (gap / 2.0 + 1.0 / kappa))
    rows.append(np.arange(nx * ny))
    cols.append(np.arange(nx * ny))
    vals.append(diag)
    a = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(nx * ny, nx * ny)).tocsc()
    return a


def _source(g: CleftGeometry, centre, sub: int = 4) -> np.ndarray:
    """Fraction of one channel's release landing in each cell (sums to 1)."""
    nx, ny, dx, dy = _grid(g)
    off = (np.arange(sub) + 0.5) / sub
    x = ((np.arange(nx)[:, None] + off[None, :]) * dx).ravel()
    y = ((np.arange(ny)[:, None] + off[None, :]) * dy).ravel()
    r2 = (x[:, None] - centre[0]) ** 2 + (y[None, :] - centre[1]) ** 2
    inside = (r2 <= (g.source_diameter / 2.0) ** 2).reshape(nx, sub, ny, sub)
    w = inside.sum(axis=(1, 3)).astype(float).ravel()
    return w / w.sum()


def _scale(g: CleftGeometry) -> float:
    """µM per unit of the solve: flux (mol/s) / (D h) in SI, to µM."""
    flux = g.current * 1e-12 / (2.0 * F_FARADAY)
    return flux / (g.d_ca * 1e-4 * g.height * 1e-9) * 1e3


def _sample(g: CleftGeometry, u: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Bilinear value of cell-centred fields ``u`` (cells, k) at points."""
    nx, ny, dx, dy = _grid(g)
    grid = u.reshape(nx, ny, -1)
    fx = np.clip(pts[:, 0] / dx - 0.5, 0, nx - 1.000001)
    fy = np.clip(pts[:, 1] / dy - 0.5, 0, ny - 1.000001)
    i, j = fx.astype(int), fy.astype(int)
    tx, ty = (fx - i)[:, None], (fy - j)[:, None]
    i1, j1 = np.minimum(i + 1, nx - 1), np.minimum(j + 1, ny - 1)
    return ((1 - tx) * (1 - ty) * grid[i, j] + tx * (1 - ty) * grid[i1, j]
            + (1 - tx) * ty * grid[i, j1] + tx * ty * grid[i1, j1])


@lru_cache(maxsize=8)
def _solved(g: CleftGeometry):
    lu = splu(_operator(g))
    src = np.stack([_source(g, c) for c in couplon(g).c_sites], axis=1)
    return lu.solve(src) * _scale(g), src


def coupling_matrix(g: CleftGeometry | None = None) -> np.ndarray:
    """G[i, j]: µM at C channel i's centre while C channel j is open
    (diagonal: a channel's own release). Memoised per geometry."""
    g = g or CleftGeometry.from_parameters()
    u, _ = _solved(g)
    return _sample(g, u, couplon(g).c_sites)


@lru_cache(maxsize=8)
def _solved_v(g: CleftGeometry):
    src = np.stack([_source(g, c) for c in couplon(g).v_sites], axis=1)
    return splu(_operator(g)).solve(src) * _scale(g)


def v_coupling_matrix(v_current: float, g: CleftGeometry | None = None
                      ) -> np.ndarray:
    """H[i, k]: µM at C channel i's centre while V channel k is open and
    passing ``v_current`` pA (Stern's V channels are not Ca2+-gated, so
    nothing is needed at the V sites)."""
    g = g or CleftGeometry.from_parameters()
    return _sample(g, _solved_v(g), couplon(g).c_sites) * (v_current / g.current)


def field(open_mask, g: CleftGeometry | None = None) -> np.ndarray:
    """The cleft Ca2+ above background (µM), shape (nx, ny), for the C
    channels ``open_mask`` (bool per C channel)."""
    g = g or CleftGeometry.from_parameters()
    u, _ = _solved(g)
    nx, ny, _, _ = _grid(g)
    return (u @ np.asarray(open_mask, float)).reshape(nx, ny)


def nearest_coupling(g_matrix: np.ndarray) -> float:
    """The largest coupling between two different channels, µM."""
    off = g_matrix - np.diag(np.diag(g_matrix))
    return float(off.max())
