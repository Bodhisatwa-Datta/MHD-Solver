"""Initial conditions for standard hydrodynamic tests."""

from __future__ import annotations

import numpy as np


def sod_shock_tube(x: np.ndarray, discontinuity: float = 0.5) -> np.ndarray:
    """Return Sod's standard primitive states ``(rho, velocity, pressure)``."""
    left = x < discontinuity
    return np.stack(
        (np.where(left, 1.0, 0.125), np.zeros_like(x), np.where(left, 1.0, 0.1))
    )


def smooth_density_wave(
    x: np.ndarray,
    time: float = 0.0,
    velocity: float = 1.0,
    amplitude: float = 0.2,
    pressure: float = 1.0,
    wavelength: float = 1.0,
) -> np.ndarray:
    """Return an exact periodic entropy wave advected by uniform flow.

    Uniform velocity and pressure make the sinusoidal density profile an exact
    Euler solution: ``rho(x, t) = rho(x - velocity*t, 0)``.
    """
    phase = 2.0 * np.pi * (x - velocity * time) / wavelength
    density = 1.0 + amplitude * np.sin(phase)
    return np.stack(
        (density, np.full_like(x, velocity), np.full_like(x, pressure))
    )
