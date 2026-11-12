"""Compact JSON encodings shared by the precompute (``export``) and the API.

Coordinates and mode vectors are stored as little-endian int16 in base64,
with a stated scale, so a tetramer's C-alpha trace is ~70 kB instead of
~400 kB of decimal text. ``decode_i16`` is the exact inverse used by the
tests and by the reduced-model fixture, so Python and the browser always
work on the same (quantised) numbers.
"""

from __future__ import annotations

import base64

import numpy as np

__all__ = ["COORD_SCALE", "encode_i16", "decode_i16", "encode_coords",
           "decode_coords"]

#: Å per integer step of a stored coordinate (0.01 Å; |x| < 327 Å fits int16).
COORD_SCALE = 0.01


def encode_i16(values: np.ndarray, scale: float) -> str:
    q = np.rint(np.asarray(values, np.float64) / scale)
    if q.size and np.abs(q).max() > 32767:
        raise ValueError(f"value out of int16 range at scale {scale}")
    return base64.b64encode(q.astype("<i2").tobytes()).decode("ascii")


def decode_i16(text: str, scale: float) -> np.ndarray:
    return np.frombuffer(base64.b64decode(text), dtype="<i2").astype(np.float64) * scale


def encode_coords(xyz: np.ndarray) -> str:
    """``(n, 3)`` Å -> base64 int16 at :data:`COORD_SCALE`."""
    return encode_i16(np.asarray(xyz).reshape(-1), COORD_SCALE)


def decode_coords(text: str) -> np.ndarray:
    return decode_i16(text, COORD_SCALE).reshape(-1, 3)
