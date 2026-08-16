"""Piecewise-linear reconstruction for finite-volume states."""

from __future__ import annotations

import numpy as np


def minmod(left_difference: np.ndarray, right_difference: np.ndarray) -> np.ndarray:
    """Return the component-wise minmod-limited slope."""
    same_sign = left_difference * right_difference > 0.0
    return np.where(
        same_sign,
        np.sign(left_difference) * np.minimum(np.abs(left_difference), np.abs(right_difference)),
        0.0,
    )


def muscl_interface_states(primitive_with_ghosts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Reconstruct primitive states at interfaces using two ghosts per side.

    A side whose reconstructed density or pressure is non-positive locally
    reverts to its cell average. This suppresses the offending slope rather than
    clipping a physical variable.
    """
    values = primitive_with_ghosts
    slopes = np.zeros_like(values)
    slopes[:, 1:-1] = minmod(
        values[:, 1:-1] - values[:, :-2], values[:, 2:] - values[:, 1:-1]
    )
    left_cells = values[:, 1:-2]
    right_cells = values[:, 2:-1]
    left = left_cells + 0.5 * slopes[:, 1:-2]
    right = right_cells - 0.5 * slopes[:, 2:-1]
    bad_left = (left[0] <= 0.0) | (left[2] <= 0.0)
    bad_right = (right[0] <= 0.0) | (right[2] <= 0.0)
    left[:, bad_left] = left_cells[:, bad_left]
    right[:, bad_right] = right_cells[:, bad_right]
    return left, right

