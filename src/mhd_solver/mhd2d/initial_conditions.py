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


def magnetic_rotor(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return the standard fast magnetic-rotor initial condition.

    A dense rotor of radius 0.1 spins with angular velocity 20. Density and
    velocity taper linearly to the ambient state between radii 0.1 and 0.115.
    The corresponding configuration uses ``gamma=1.4`` and outflow boundaries.
    """
    xx, yy = np.meshgrid(x, y)
    delta_x, delta_y = xx - 0.5, yy - 0.5
    radius = np.sqrt(delta_x**2 + delta_y**2)
    inner_radius, outer_radius = 0.1, 0.115
    taper = np.where(
        radius <= inner_radius,
        1.0,
        np.where(
            radius < outer_radius,
            (outer_radius - radius) / (outer_radius - inner_radius),
            0.0,
        ),
    )
    density = 1.0 + 9.0 * taper
    velocity_x = -20.0 * taper * delta_y
    velocity_y = 20.0 * taper * delta_x
    zeros = np.zeros_like(xx)
    return np.stack(
        (
            density,
            velocity_x,
            velocity_y,
            zeros,
            np.ones_like(xx),
            np.full_like(xx, 5.0 / np.sqrt(4.0 * np.pi)),
            zeros,
            zeros,
            zeros,
        )
    )
