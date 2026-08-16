"""Boundary fills for two-dimensional cell-centred states."""

from __future__ import annotations

import numpy as np


def fill_boundaries(values: np.ndarray, boundary: str, ghosts: int = 2) -> np.ndarray:
    """Fill periodic or transmissive ghost cells on both coordinate axes."""
    if boundary == "periodic":
        return np.pad(values, ((0, 0), (ghosts, ghosts), (ghosts, ghosts)), mode="wrap")
    if boundary == "outflow":
        return np.pad(values, ((0, 0), (ghosts, ghosts), (ghosts, ghosts)), mode="edge")
    raise ValueError("boundary must be 'periodic' or 'outflow'")
