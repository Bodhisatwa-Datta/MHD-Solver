"""Physical and approximate Riemann fluxes for Euler flow."""

from __future__ import annotations

import numpy as np

from mhd_solver.common.eos import primitive_to_conserved, sound_speed


def euler_flux(primitive: np.ndarray, gamma: float) -> np.ndarray:
    """Return the physical Euler flux for one or more primitive states."""
    density, velocity, pressure = primitive
    energy = primitive_to_conserved(primitive, gamma)[2]
    return np.stack(
        (density * velocity, density * velocity**2 + pressure, velocity * (energy + pressure))
    )


def hll_flux(left: np.ndarray, right: np.ndarray, gamma: float) -> np.ndarray:
    """Compute HLL fluxes with extremal acoustic signal-speed estimates."""
    conserved_left = primitive_to_conserved(left, gamma)
    conserved_right = primitive_to_conserved(right, gamma)
    flux_left = euler_flux(left, gamma)
    flux_right = euler_flux(right, gamma)
    speed_left = np.minimum(
        left[1] - sound_speed(left, gamma), right[1] - sound_speed(right, gamma)
    )
    speed_right = np.maximum(
        left[1] + sound_speed(left, gamma), right[1] + sound_speed(right, gamma)
    )
    denominator = speed_right - speed_left
    middle = (
        speed_right * flux_left
        - speed_left * flux_right
        + speed_left * speed_right * (conserved_right - conserved_left)
    ) / np.where(denominator > 0.0, denominator, 1.0)
    return np.where(
        speed_left >= 0.0, flux_left, np.where(speed_right <= 0.0, flux_right, middle)
    )
