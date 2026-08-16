"""Derived fields and integral diagnostics for two-dimensional MHD."""

from __future__ import annotations

import numpy as np


def current_density_z(
    primitive: np.ndarray, dx: float, dy: float, periodic: bool = True
) -> np.ndarray:
    """Return the out-of-plane current ``J_z = dB_y/dx - dB_x/dy``."""
    bx, by = primitive[5], primitive[6]
    if periodic:
        derivative_by_x = (np.roll(by, -1, axis=1) - np.roll(by, 1, axis=1)) / (2.0 * dx)
        derivative_bx_y = (np.roll(bx, -1, axis=0) - np.roll(bx, 1, axis=0)) / (2.0 * dy)
    else:
        derivative_by_x = np.gradient(by, dx, axis=1, edge_order=2)
        derivative_bx_y = np.gradient(bx, dy, axis=0, edge_order=2)
    return derivative_by_x - derivative_bx_y


def kinetic_energy(primitive: np.ndarray, dx: float, dy: float) -> float:
    density = primitive[0]
    velocity_squared = primitive[1] ** 2 + primitive[2] ** 2 + primitive[3] ** 2
    return float(0.5 * np.sum(density * velocity_squared) * dx * dy)


def magnetic_energy(primitive: np.ndarray, dx: float, dy: float) -> float:
    field_squared = primitive[5] ** 2 + primitive[6] ** 2 + primitive[7] ** 2
    return float(0.5 * np.sum(field_squared) * dx * dy)
