"""Boundary fills shared by one-dimensional finite-volume solvers."""

from __future__ import annotations

import numpy as np


def apply_outflow_boundaries(conserved: np.ndarray, ghost_cells: int = 2) -> np.ndarray:
    """Pad cell averages with transmissive constant-extrapolation ghosts."""
    return np.pad(conserved, ((0, 0), (ghost_cells, ghost_cells)), mode="edge")


def apply_periodic_boundaries(conserved: np.ndarray, ghost_cells: int = 2) -> np.ndarray:
    """Pad cell averages by wrapping values from the opposite domain edge."""
    return np.pad(conserved, ((0, 0), (ghost_cells, ghost_cells)), mode="wrap")
