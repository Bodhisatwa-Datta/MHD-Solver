"""Initial conditions for one-dimensional ideal-MHD benchmarks."""

from __future__ import annotations

import numpy as np


def brio_wu_shock_tube(x: np.ndarray, discontinuity: float = 0.0) -> np.ndarray:
    """Return the standard Brio-Wu primitive states.

    The corresponding configuration uses ``gamma=2`` and ``B_x=0.75``.
    Primitive ordering is ``(rho, vx, vy, vz, p, By, Bz)``.
    """
    left = x < discontinuity
    zeros = np.zeros_like(x)
    return np.stack(
        (
            np.where(left, 1.0, 0.125),
            zeros,
            zeros,
            zeros,
            np.where(left, 1.0, 0.1),
            np.where(left, 1.0, -1.0),
            zeros,
        )
    )


def circularly_polarized_alfven_wave(
    x: np.ndarray,
    time: float = 0.0,
    density: float = 1.0,
    pressure: float = 0.1,
    longitudinal_field: float = 1.0,
    transverse_amplitude: float = 0.1,
    wavelength: float = 1.0,
) -> np.ndarray:
    """Return an exact right-travelling circularly polarized Alfvén wave.

    Constant density, gas pressure, and transverse-field magnitude make this a
    finite-amplitude exact ideal-MHD solution. For positive ``B_x`` it propagates
    in the positive x direction with speed ``B_x/sqrt(rho)``.
    """
    alfven_speed = longitudinal_field / np.sqrt(density)
    phase = 2.0 * np.pi * (x - alfven_speed * time) / wavelength
    field_y = transverse_amplitude * np.cos(phase)
    field_z = transverse_amplitude * np.sin(phase)
    velocity_y = -field_y / np.sqrt(density)
    velocity_z = -field_z / np.sqrt(density)
    return np.stack(
        (
            np.full_like(x, density),
            np.zeros_like(x),
            velocity_y,
            velocity_z,
            np.full_like(x, pressure),
            field_y,
            field_z,
        )
    )
