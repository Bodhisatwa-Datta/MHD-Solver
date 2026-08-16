"""Uniform-resistivity terms for the two-dimensional MHD equations."""

from __future__ import annotations

import numpy as np


def _first_derivative(
    values: np.ndarray, spacing: float, axis: int, boundary: str
) -> np.ndarray:
    if boundary == "periodic":
        return (
            np.roll(values, -1, axis=axis) - np.roll(values, 1, axis=axis)
        ) / (2.0 * spacing)
    return np.gradient(values, spacing, axis=axis, edge_order=2)


def _laplacian(
    values: np.ndarray, dx: float, dy: float, boundary_x: str, boundary_y: str
) -> np.ndarray:
    if boundary_x == "periodic":
        second_x = (
            np.roll(values, -1, axis=1) - 2.0 * values + np.roll(values, 1, axis=1)
        ) / dx**2
    else:
        padded_x = np.pad(values, ((0, 0), (1, 1)), mode="edge")
        second_x = (padded_x[:, 2:] - 2.0 * values + padded_x[:, :-2]) / dx**2
    if boundary_y == "periodic":
        second_y = (
            np.roll(values, -1, axis=0) - 2.0 * values + np.roll(values, 1, axis=0)
        ) / dy**2
    else:
        padded_y = np.pad(values, ((1, 1), (0, 0)), mode="edge")
        second_y = (padded_y[2:] - 2.0 * values + padded_y[:-2]) / dy**2
    return second_x + second_y


def resistive_rhs(
    primitive: np.ndarray,
    resistivity: float,
    dx: float,
    dy: float,
    boundary: str | tuple[str, str],
) -> np.ndarray:
    """Return uniform-resistivity contributions to conserved-variable rates.

    Magnetic components receive ``eta*laplacian(B)``. Total energy receives the
    conservative flux divergence ``-div[eta*(J cross B)]`` so periodic-domain
    total energy is conserved while magnetic energy is converted to heat.
    """
    rhs = np.zeros_like(primitive)
    if resistivity == 0.0:
        return rhs
    boundary_x, boundary_y = (
        (boundary, boundary) if isinstance(boundary, str) else boundary
    )
    bx, by, bz = primitive[5], primitive[6], primitive[7]
    current_x = _first_derivative(bz, dy, axis=0, boundary=boundary_y)
    current_y = -_first_derivative(bz, dx, axis=1, boundary=boundary_x)
    current_z = _first_derivative(
        by, dx, axis=1, boundary=boundary_x
    ) - _first_derivative(
        bx,
        dy,
        axis=0,
        boundary=boundary_y,
    )
    poynting_x = resistivity * (current_y * bz - current_z * by)
    poynting_y = resistivity * (current_z * bx - current_x * bz)
    rhs[4] = -(
        _first_derivative(poynting_x, dx, axis=1, boundary=boundary_x)
        + _first_derivative(poynting_y, dy, axis=0, boundary=boundary_y)
    )
    rhs[5] = resistivity * _laplacian(bx, dx, dy, boundary_x, boundary_y)
    rhs[6] = resistivity * _laplacian(by, dx, dy, boundary_x, boundary_y)
    rhs[7] = resistivity * _laplacian(bz, dx, dy, boundary_x, boundary_y)
    return rhs
