"""Initial conditions for validating the two-dimensional MHD framework."""

from __future__ import annotations

import numpy as np


def divergence_perturbation(
    x: np.ndarray, y: np.ndarray, amplitude: float = 1.0e-3
) -> np.ndarray:
    """Return a small sinusoidal ``div(B)`` perturbation for GLM testing."""
    xx, _ = np.meshgrid(x, y)
    shape = xx.shape
    zeros = np.zeros(shape)
    return np.stack(
        (
            np.ones(shape),
            zeros,
            zeros,
            zeros,
            np.ones(shape),
            amplitude * np.sin(2.0 * np.pi * xx),
            zeros,
            zeros,
            zeros,
        )
    )
