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


def orszag_tang_vortex(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return the standard periodic Orszag-Tang vortex in ``mu_0 = 1`` units.

    The normalization uses ``rho=25/(36*pi)``, ``p=5/(12*pi)``, and magnetic
    amplitudes divided by ``sqrt(4*pi)``, with ``gamma=5/3``.
    """
    xx, yy = np.meshgrid(x, y)
    shape = xx.shape
    zeros = np.zeros(shape)
    magnetic_scale = 1.0 / np.sqrt(4.0 * np.pi)
    return np.stack(
        (
            np.full(shape, 25.0 / (36.0 * np.pi)),
            -np.sin(2.0 * np.pi * yy),
            np.sin(2.0 * np.pi * xx),
            zeros,
            np.full(shape, 5.0 / (12.0 * np.pi)),
            -magnetic_scale * np.sin(2.0 * np.pi * yy),
            magnetic_scale * np.sin(4.0 * np.pi * xx),
            zeros,
            zeros,
        )
    )
