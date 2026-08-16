"""Ideal-gas equation-of-state conversions for 1D Euler flow."""

from __future__ import annotations

import numpy as np


class UnphysicalStateError(RuntimeError):
    """Raised when a conserved or primitive state is not physically admissible."""


def _require_finite(values: np.ndarray, name: str) -> None:
    if not np.all(np.isfinite(values)):
        bad = np.argwhere(~np.isfinite(values))[0]
        raise UnphysicalStateError(f"{name} contains a non-finite value at index {tuple(bad)}")


def primitive_to_conserved(primitive: np.ndarray, gamma: float) -> np.ndarray:
    """Convert ``(density, velocity, pressure)`` to conserved variables."""
    primitive = np.asarray(primitive, dtype=float)
    _require_finite(primitive, "primitive state")
    density, velocity, pressure = primitive
    if np.any(density <= 0.0):
        raise UnphysicalStateError(f"density must be positive; minimum={density.min():.6e}")
    if np.any(pressure <= 0.0):
        raise UnphysicalStateError(f"pressure must be positive; minimum={pressure.min():.6e}")
    energy = pressure / (gamma - 1.0) + 0.5 * density * velocity**2
    return np.stack((density, density * velocity, energy))


def conserved_to_primitive(conserved: np.ndarray, gamma: float) -> np.ndarray:
    """Convert conserved variables to ``(density, velocity, pressure)``."""
    conserved = np.asarray(conserved, dtype=float)
    _require_finite(conserved, "conserved state")
    density, momentum, energy = conserved
    if np.any(density <= 0.0):
        raise UnphysicalStateError(f"density must be positive; minimum={density.min():.6e}")
    velocity = momentum / density
    pressure = (gamma - 1.0) * (energy - 0.5 * momentum * velocity)
    if np.any(pressure <= 0.0):
        raise UnphysicalStateError(f"pressure must be positive; minimum={pressure.min():.6e}")
    primitive = np.stack((density, velocity, pressure))
    _require_finite(primitive, "primitive state")
    return primitive


def sound_speed(primitive: np.ndarray, gamma: float) -> np.ndarray:
    """Return the ideal-gas adiabatic sound speed."""
    density, _, pressure = primitive
    return np.sqrt(gamma * pressure / density)

