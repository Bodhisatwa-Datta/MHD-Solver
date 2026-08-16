"""Initial conditions for magnetic-reconnection validation problems."""

from __future__ import annotations

import numpy as np


def harris_current_sheet(
    x: np.ndarray,
    y: np.ndarray,
    *,
    magnetic_field: float = 1.0,
    half_width: float = 0.05,
    background_pressure: float = 0.2,
    background_density: float = 1.0,
) -> np.ndarray:
    r"""Return an unperturbed, pressure-balanced Harris current sheet.

    The reversing field is ``B_x = B_0 tanh(y / L)``. Gas pressure supplies
    the missing magnetic pressure, and density follows pressure so that the
    temperature ``p/rho`` is uniform. The state is an ideal-MHD equilibrium;
    no reconnection perturbation or guide field is included.
    """
    if magnetic_field <= 0.0:
        raise ValueError("magnetic_field must be positive")
    if half_width <= 0.0:
        raise ValueError("half_width must be positive")
    if background_pressure <= 0.0 or background_density <= 0.0:
        raise ValueError("background pressure and density must be positive")

    xx, yy = np.meshgrid(x, y)
    normalized_y = yy / half_width
    hyperbolic_tangent = np.tanh(normalized_y)
    field_x = magnetic_field * hyperbolic_tangent
    pressure = background_pressure + 0.5 * magnetic_field**2 * (
        1.0 - hyperbolic_tangent**2
    )
    temperature = background_pressure / background_density
    density = pressure / temperature
    zeros = np.zeros_like(xx)
    return np.stack(
        (
            density,
            zeros,
            zeros,
            zeros,
            pressure,
            field_x,
            zeros,
            zeros,
            zeros,
        )
    )
