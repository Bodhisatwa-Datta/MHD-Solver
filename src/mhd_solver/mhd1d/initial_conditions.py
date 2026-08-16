"""Initial conditions for one-dimensional ideal-MHD benchmarks."""

from __future__ import annotations

import numpy as np


def brio_wu_shock_tube(x: np.ndarray, discontinuity: float = 0.0) -> np.ndarray:
    """Return the standard Brio-Wu primitive states.

    The corresponding configuration uses ``gamma=2`` and ``B_x=0.75``.
    Primitive ordering is ``(rho, vx, vy, vz, p, By, Bz)``.
    """
    left = x < discontinuity
    zeros = np.zeros_like(x)
    return np.stack(
        (
            np.where(left, 1.0, 0.125),
            zeros,
            zeros,
            zeros,
            np.where(left, 1.0, 0.1),
            np.where(left, 1.0, -1.0),
            zeros,
        )
    )
