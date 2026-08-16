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
    values: np.ndarray, dx: float, dy: float, boundary: str
) -> np.ndarray:
    if boundary == "periodic":
        return (
            (np.roll(values, -1, axis=1) - 2.0 * values + np.roll(values, 1, axis=1))
            / dx**2
            + (np.roll(values, -1, axis=0) - 2.0 * values + np.roll(values, 1, axis=0))
            / dy**2
        )
    padded = np.pad(values, ((1, 1), (1, 1)), mode="edge")
    return (
        (padded[1:-1, 2:] - 2.0 * values + padded[1:-1, :-2]) / dx**2
        + (padded[2:, 1:-1] - 2.0 * values + padded[:-2, 1:-1]) / dy**2
    )


def resistive_rhs(
    primitive: np.ndarray,
    resistivity: float,
    dx: float,
    dy: float,
    boundary: str,
) -> np.ndarray:
    """Return uniform-resistivity contributions to conserved-variable rates.

    Magnetic components receive ``eta*laplacian(B)``. Total energy receives the
    conservative flux divergence ``-div[eta*(J cross B)]`` so periodic-domain
    total energy is conserved while magnetic energy is converted to heat.
    """
    rhs = np.zeros_like(primitive)
    if resistivity == 0.0:
        return rhs
    bx, by, bz = primitive[5], primitive[6], primitive[7]
    current_x = _first_derivative(bz, dy, axis=0, boundary=boundary)
    current_y = -_first_derivative(bz, dx, axis=1, boundary=boundary)
    current_z = _first_derivative(by, dx, axis=1, boundary=boundary) - _first_derivative(
        bx, dy, axis=0, boundary=boundary
    )
    poynting_x = resistivity * (current_y * bz - current_z * by)
    poynting_y = resistivity * (current_z * bx - current_x * bz)
    rhs[4] = -(
        _first_derivative(poynting_x, dx, axis=1, boundary=boundary)
        + _first_derivative(poynting_y, dy, axis=0, boundary=boundary)
    )
    rhs[5] = resistivity * _laplacian(bx, dx, dy, boundary)
    rhs[6] = resistivity * _laplacian(by, dx, dy, boundary)
    rhs[7] = resistivity * _laplacian(bz, dx, dy, boundary)
    return rhs
