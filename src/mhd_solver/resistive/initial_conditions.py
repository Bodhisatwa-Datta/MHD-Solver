"""Initial and analytical states for resistive-MHD validation."""

from __future__ import annotations

import numpy as np


def magnetic_diffusion_wave(
    x: np.ndarray,
    y: np.ndarray,
    time: float = 0.0,
    resistivity: float = 0.1,
    amplitude: float = 1.0e-3,
) -> np.ndarray:
    """Return a weak sinusoidal field with analytical resistive decay.

    ``B_y = A exp(-eta*k^2*t) sin(k*x)`` with ``k=2*pi``. The weak amplitude
    isolates magnetic diffusion from Lorentz-force feedback in the full MHD run.
    """
    xx, _ = np.meshgrid(x, y)
    shape = xx.shape
    zeros = np.zeros(shape)
    decay = np.exp(-resistivity * (2.0 * np.pi) ** 2 * time)
    field_y = amplitude * decay * np.sin(2.0 * np.pi * xx)
    return np.stack(
        (
            np.ones(shape),
            zeros,
            zeros,
            zeros,
            np.ones(shape),
            zeros,
            field_y,
            zeros,
            zeros,
        )
    )
