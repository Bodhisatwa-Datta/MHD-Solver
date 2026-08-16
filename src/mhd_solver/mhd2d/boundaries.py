"""Boundary fills for two-dimensional cell-centred states."""

from __future__ import annotations

import numpy as np


def fill_boundaries(
    values: np.ndarray,
    boundary_x: str,
    boundary_y: str | None = None,
    ghosts: int = 2,
) -> np.ndarray:
    """Fill periodic or transmissive ghosts independently in x and y."""
    boundary_y = boundary_x if boundary_y is None else boundary_y
    if boundary_x not in ("periodic", "outflow") or boundary_y not in (
        "periodic",
        "outflow",
    ):
        raise ValueError("boundaries must be 'periodic' or 'outflow'")
    mode_y = "wrap" if boundary_y == "periodic" else "edge"
    mode_x = "wrap" if boundary_x == "periodic" else "edge"
    padded = np.pad(values, ((0, 0), (ghosts, ghosts), (0, 0)), mode=mode_y)
    return np.pad(padded, ((0, 0), (0, 0), (ghosts, ghosts)), mode=mode_x)
