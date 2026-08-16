"""Initial conditions for standard hydrodynamic tests."""

from __future__ import annotations

import numpy as np


def sod_shock_tube(x: np.ndarray, discontinuity: float = 0.5) -> np.ndarray:
    """Return Sod's standard primitive states ``(rho, velocity, pressure)``."""
    left = x < discontinuity
    return np.stack(
        (np.where(left, 1.0, 0.125), np.zeros_like(x), np.where(left, 1.0, 0.1))
    )

