"""State conversions, fluxes, and wave speeds for one-dimensional ideal MHD.

The longitudinal magnetic field ``B_x`` is spatially constant in one dimension
and is passed separately rather than included in the evolved state.
"""

from __future__ import annotations

import numpy as np

from mhd_solver.common.eos import UnphysicalStateError


def _require_finite(values: np.ndarray, name: str) -> None:
    if not np.all(np.isfinite(values)):
        bad = np.argwhere(~np.isfinite(values))[0]
        raise UnphysicalStateError(f"{name} contains a non-finite value at index {tuple(bad)}")


def primitive_to_conserved(
    primitive: np.ndarray, gamma: float, longitudinal_field: float
) -> np.ndarray:
    """Convert ``(rho,vx,vy,vz,p,By,Bz)`` to the seven evolved variables."""
    primitive = np.asarray(primitive, dtype=float)
    _require_finite(primitive, "MHD primitive state")
    density, velocity_x, velocity_y, velocity_z, pressure, field_y, field_z = primitive
    if np.any(density <= 0.0):
        raise UnphysicalStateError(f"density must be positive; minimum={density.min():.6e}")
    if np.any(pressure <= 0.0):
        raise UnphysicalStateError(f"pressure must be positive; minimum={pressure.min():.6e}")
    kinetic_energy = 0.5 * density * (
        velocity_x**2 + velocity_y**2 + velocity_z**2
    )
    magnetic_energy = 0.5 * (longitudinal_field**2 + field_y**2 + field_z**2)
    total_energy = pressure / (gamma - 1.0) + kinetic_energy + magnetic_energy
    return np.stack(
        (
            density,
            density * velocity_x,
            density * velocity_y,
            density * velocity_z,
            total_energy,
            field_y,
            field_z,
        )
    )


def conserved_to_primitive(
    conserved: np.ndarray, gamma: float, longitudinal_field: float
) -> np.ndarray:
    """Convert the seven evolved variables to ``(rho,vx,vy,vz,p,By,Bz)``."""
    conserved = np.asarray(conserved, dtype=float)
    _require_finite(conserved, "MHD conserved state")
    density, momentum_x, momentum_y, momentum_z, total_energy, field_y, field_z = conserved
    if np.any(density <= 0.0):
        raise UnphysicalStateError(f"density must be positive; minimum={density.min():.6e}")
    velocity_x = momentum_x / density
    velocity_y = momentum_y / density
    velocity_z = momentum_z / density
    kinetic_energy = 0.5 * (
        momentum_x * velocity_x + momentum_y * velocity_y + momentum_z * velocity_z
    )
    magnetic_energy = 0.5 * (longitudinal_field**2 + field_y**2 + field_z**2)
    pressure = (gamma - 1.0) * (total_energy - kinetic_energy - magnetic_energy)
    if np.any(pressure <= 0.0):
        raise UnphysicalStateError(f"pressure must be positive; minimum={pressure.min():.6e}")
    primitive = np.stack(
        (density, velocity_x, velocity_y, velocity_z, pressure, field_y, field_z)
    )
    _require_finite(primitive, "MHD primitive state")
    return primitive


def total_pressure(primitive: np.ndarray, longitudinal_field: float) -> np.ndarray:
    """Return gas plus magnetic pressure in units with magnetic permeability one."""
    return primitive[4] + 0.5 * (
        longitudinal_field**2 + primitive[5] ** 2 + primitive[6] ** 2
    )


def fast_magnetosonic_speed(
    primitive: np.ndarray, gamma: float, longitudinal_field: float
) -> np.ndarray:
    """Return the fast magnetosonic speed along the x direction."""
    density = primitive[0]
    sound_speed_squared = gamma * primitive[4] / density
    alfven_speed_squared = (
        longitudinal_field**2 + primitive[5] ** 2 + primitive[6] ** 2
    ) / density
    longitudinal_alfven_squared = longitudinal_field**2 / density
    discriminant = (sound_speed_squared + alfven_speed_squared) ** 2 - (
        4.0 * sound_speed_squared * longitudinal_alfven_squared
    )
    # The analytic discriminant is non-negative; guard only roundoff below zero.
    discriminant = np.maximum(discriminant, 0.0)
    return np.sqrt(
        0.5
        * (sound_speed_squared + alfven_speed_squared + np.sqrt(discriminant))
    )


def physical_flux(
    primitive: np.ndarray, gamma: float, longitudinal_field: float
) -> np.ndarray:
    """Return the x-directed ideal-MHD physical flux."""
    density, velocity_x, velocity_y, velocity_z, _, field_y, field_z = primitive
    conserved = primitive_to_conserved(primitive, gamma, longitudinal_field)
    energy = conserved[4]
    pressure_total = total_pressure(primitive, longitudinal_field)
    velocity_dot_field = (
        velocity_x * longitudinal_field + velocity_y * field_y + velocity_z * field_z
    )
    return np.stack(
        (
            density * velocity_x,
            density * velocity_x**2 + pressure_total - longitudinal_field**2,
            density * velocity_x * velocity_y - longitudinal_field * field_y,
            density * velocity_x * velocity_z - longitudinal_field * field_z,
            (energy + pressure_total) * velocity_x
            - longitudinal_field * velocity_dot_field,
            velocity_x * field_y - velocity_y * longitudinal_field,
            velocity_x * field_z - velocity_z * longitudinal_field,
        )
    )


def hll_flux(
    left: np.ndarray, right: np.ndarray, gamma: float, longitudinal_field: float
) -> np.ndarray:
    """Return the two-wave HLL ideal-MHD interface flux."""
    conserved_left = primitive_to_conserved(left, gamma, longitudinal_field)
    conserved_right = primitive_to_conserved(right, gamma, longitudinal_field)
    flux_left = physical_flux(left, gamma, longitudinal_field)
    flux_right = physical_flux(right, gamma, longitudinal_field)
    speed_left = np.minimum(
        left[1] - fast_magnetosonic_speed(left, gamma, longitudinal_field),
        right[1] - fast_magnetosonic_speed(right, gamma, longitudinal_field),
    )
    speed_right = np.maximum(
        left[1] + fast_magnetosonic_speed(left, gamma, longitudinal_field),
        right[1] + fast_magnetosonic_speed(right, gamma, longitudinal_field),
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
