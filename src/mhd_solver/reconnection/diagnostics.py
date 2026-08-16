"""Diagnostics for two-dimensional Harris-sheet reconnection."""

from __future__ import annotations

import numpy as np


def current_density_z(primitive: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Return J_z for periodic x and outflow y boundaries."""
    bx, by = primitive[5], primitive[6]
    derivative_by_x = (
        np.roll(by, -1, axis=1) - np.roll(by, 1, axis=1)
    ) / (2.0 * dx)
    derivative_bx_y = np.gradient(bx, dy, axis=0, edge_order=2)
    return derivative_by_x - derivative_bx_y


def reconnected_flux(
    primitive: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    *,
    o_point: float = 0.0,
    x_point: float = 0.5,
) -> float:
    r"""Return magnetic flux between the prescribed O- and X-points.

    On the sheet midplane, ``A_z(X)-A_z(O) = -integral B_y dx``. The absolute
    midpoint-rule integral is used because finite-volume fields are cell
    centred. The requested interval must not cross the periodic boundary.
    """
    if not x[0] < x_point < x[-1] or o_point > x[0]:
        raise ValueError("diagnostic expects the O-point at the left periodic boundary")
    dx = float(x[1] - x[0])
    middle_rows = np.argsort(np.abs(y))[:2]
    midplane_by = np.mean(primitive[6, middle_rows, :], axis=0)
    interval = (x >= o_point) & (x < x_point)
    return float(abs(np.sum(midplane_by[interval]) * dx))


def reconnection_electric_field(
    primitive: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    resistivity: float,
    *,
    x_point: float = 0.5,
) -> float:
    r"""Return ``E_z = -(v cross B)_z + eta J_z`` at the seeded X-point."""
    if resistivity < 0.0:
        raise ValueError("resistivity must be non-negative")
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    current = current_density_z(primitive, dx, dy)
    electric_z = -(
        primitive[1] * primitive[6] - primitive[2] * primitive[5]
    ) + resistivity * current
    columns = np.argsort(np.abs(x - x_point))[:2]
    rows = np.argsort(np.abs(y))[:2]
    return float(np.mean(electric_z[np.ix_(rows, columns)]))
