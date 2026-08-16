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


def perturbed_harris_sheet(
    x: np.ndarray,
    y: np.ndarray,
    *,
    magnetic_field: float = 1.0,
    half_width: float = 0.05,
    background_pressure: float = 0.2,
    background_density: float = 1.0,
    perturbation_amplitude: float = 0.01,
    perturbation_width: float = 0.1,
    x_point: float = 0.5,
) -> np.ndarray:
    r"""Return a Harris sheet with a divergence-free magnetic perturbation.

    A localized out-of-plane vector potential
    ``delta A_z = A_p cos[2 pi (x-x_X)/L_x] exp[-(y/w)^2]`` is sampled on
    cell centres. Its discrete curl is evaluated with the same directional
    derivatives as the divergence diagnostic, making the perturbation
    divergence-free to roundoff. The default places an X-point at ``x=0.5``
    and an O-point at the periodic x boundary.
    """
    if x.size < 4 or y.size < 4:
        raise ValueError("perturbed sheet requires at least four cells per direction")
    if not np.isfinite(perturbation_amplitude) or perturbation_amplitude == 0.0:
        raise ValueError("perturbation amplitude must be finite and non-zero")
    if not np.isfinite(x_point) or perturbation_width <= 0.0:
        raise ValueError("x_point must be finite and perturbation width positive")
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    if not np.allclose(np.diff(x), dx) or not np.allclose(np.diff(y), dy):
        raise ValueError("perturbed sheet requires a uniform Cartesian grid")

    primitive = harris_current_sheet(
        x,
        y,
        magnetic_field=magnetic_field,
        half_width=half_width,
        background_pressure=background_pressure,
        background_density=background_density,
    )
    xx, yy = np.meshgrid(x, y)
    domain_length_x = float(x[-1] - x[0] + dx)
    vector_potential = perturbation_amplitude * np.cos(
        2.0 * np.pi * (xx - x_point) / domain_length_x
    ) * np.exp(-(yy / perturbation_width) ** 2)
    perturbation_bx = np.gradient(vector_potential, dy, axis=0, edge_order=2)
    perturbation_by = -(
        np.roll(vector_potential, -1, axis=1)
        - np.roll(vector_potential, 1, axis=1)
    ) / (2.0 * dx)
    primitive[5] += perturbation_bx
    primitive[6] += perturbation_by
    return primitive
