"""Conservative equations for two-dimensional GLM ideal MHD.

Primitive ordering is ``(rho,vx,vy,vz,p,Bx,By,Bz,psi)``. Conserved ordering
replaces velocities by momenta and stores total energy at index four.
"""

from __future__ import annotations

import numpy as np

from mhd_solver.common.eos import UnphysicalStateError


def _require_finite(values: np.ndarray, name: str) -> None:
    if not np.all(np.isfinite(values)):
        bad = np.argwhere(~np.isfinite(values))[0]
        raise UnphysicalStateError(f"{name} contains a non-finite value at index {tuple(bad)}")


def primitive_to_conserved(
    primitive: np.ndarray, gamma: float, cleaning_speed: float
) -> np.ndarray:
    """Convert primitive GLM-MHD states to conserved states."""
    primitive = np.asarray(primitive, dtype=float)
    _require_finite(primitive, "2D MHD primitive state")
    density, vx, vy, vz, pressure, bx, by, bz, psi = primitive
    if np.any(density <= 0.0):
        raise UnphysicalStateError(f"density must be positive; minimum={density.min():.6e}")
    if np.any(pressure <= 0.0):
        raise UnphysicalStateError(f"pressure must be positive; minimum={pressure.min():.6e}")
    if cleaning_speed <= 0.0:
        raise ValueError("GLM cleaning speed must be positive")
    kinetic = 0.5 * density * (vx**2 + vy**2 + vz**2)
    magnetic = 0.5 * (bx**2 + by**2 + bz**2)
    cleaning = 0.5 * psi**2 / cleaning_speed**2
    energy = pressure / (gamma - 1.0) + kinetic + magnetic + cleaning
    return np.stack((density, density * vx, density * vy, density * vz, energy, bx, by, bz, psi))


def conserved_to_primitive(
    conserved: np.ndarray, gamma: float, cleaning_speed: float
) -> np.ndarray:
    """Convert conserved GLM-MHD states to primitive states with checks."""
    conserved = np.asarray(conserved, dtype=float)
    _require_finite(conserved, "2D MHD conserved state")
    density, mx, my, mz, energy, bx, by, bz, psi = conserved
    if np.any(density <= 0.0):
        raise UnphysicalStateError(f"density must be positive; minimum={density.min():.6e}")
    if cleaning_speed <= 0.0:
        raise ValueError("GLM cleaning speed must be positive")
    vx, vy, vz = mx / density, my / density, mz / density
    kinetic = 0.5 * (mx * vx + my * vy + mz * vz)
    magnetic = 0.5 * (bx**2 + by**2 + bz**2)
    cleaning = 0.5 * psi**2 / cleaning_speed**2
    pressure = (gamma - 1.0) * (energy - kinetic - magnetic - cleaning)
    if np.any(pressure <= 0.0):
        raise UnphysicalStateError(f"pressure must be positive; minimum={pressure.min():.6e}")
    primitive = np.stack((density, vx, vy, vz, pressure, bx, by, bz, psi))
    _require_finite(primitive, "2D MHD primitive state")
    return primitive


def fast_magnetosonic_speed(
    primitive: np.ndarray, gamma: float, direction: int
) -> np.ndarray:
    """Return the fast speed normal to x (0) or y (1)."""
    density = primitive[0]
    sound_squared = gamma * primitive[4] / density
    magnetic_squared = (primitive[5] ** 2 + primitive[6] ** 2 + primitive[7] ** 2) / density
    normal_field = primitive[5] if direction == 0 else primitive[6]
    normal_alfven_squared = normal_field**2 / density
    discriminant = np.maximum(
        (sound_squared + magnetic_squared) ** 2
        - 4.0 * sound_squared * normal_alfven_squared,
        0.0,
    )
    return np.sqrt(0.5 * (sound_squared + magnetic_squared + np.sqrt(discriminant)))


def physical_flux(
    primitive: np.ndarray, gamma: float, cleaning_speed: float, direction: int
) -> np.ndarray:
    """Return the x- or y-directed conservative GLM-MHD flux."""
    rho, vx, vy, vz, pressure, bx, by, bz, psi = primitive
    conserved = primitive_to_conserved(primitive, gamma, cleaning_speed)
    energy = conserved[4]
    physical_energy = energy - 0.5 * psi**2 / cleaning_speed**2
    total_pressure = pressure + 0.5 * (bx**2 + by**2 + bz**2)
    velocity_dot_field = vx * bx + vy * by + vz * bz
    if direction == 0:
        return np.stack(
            (
                rho * vx,
                rho * vx**2 + total_pressure - bx**2,
                rho * vx * vy - bx * by,
                rho * vx * vz - bx * bz,
                (physical_energy + total_pressure) * vx
                - bx * velocity_dot_field
                + psi * bx,
                psi,
                vx * by - vy * bx,
                vx * bz - vz * bx,
                cleaning_speed**2 * bx,
            )
        )
    if direction == 1:
        return np.stack(
            (
                rho * vy,
                rho * vy * vx - by * bx,
                rho * vy**2 + total_pressure - by**2,
                rho * vy * vz - by * bz,
                (physical_energy + total_pressure) * vy
                - by * velocity_dot_field
                + psi * by,
                vy * bx - vx * by,
                psi,
                vy * bz - vz * by,
                cleaning_speed**2 * by,
            )
        )
    raise ValueError("direction must be 0 (x) or 1 (y)")


def hll_flux(
    left: np.ndarray,
    right: np.ndarray,
    gamma: float,
    cleaning_speed: float,
    direction: int,
) -> np.ndarray:
    """Return a directional two-wave HLL GLM-MHD flux."""
    conserved_left = primitive_to_conserved(left, gamma, cleaning_speed)
    conserved_right = primitive_to_conserved(right, gamma, cleaning_speed)
    flux_left = physical_flux(left, gamma, cleaning_speed, direction)
    flux_right = physical_flux(right, gamma, cleaning_speed, direction)
    velocity_index = 1 if direction == 0 else 2
    fast_left = fast_magnetosonic_speed(left, gamma, direction)
    fast_right = fast_magnetosonic_speed(right, gamma, direction)
    outer_left = np.maximum(fast_left, cleaning_speed)
    outer_right = np.maximum(fast_right, cleaning_speed)
    speed_left = np.minimum(left[velocity_index] - outer_left, right[velocity_index] - outer_right)
    speed_right = np.maximum(left[velocity_index] + outer_left, right[velocity_index] + outer_right)
    denominator = speed_right - speed_left
    middle = (
        speed_right * flux_left
        - speed_left * flux_right
        + speed_left * speed_right * (conserved_right - conserved_left)
    ) / np.where(denominator > 0.0, denominator, 1.0)
    return np.where(speed_left >= 0.0, flux_left, np.where(speed_right <= 0.0, flux_right, middle))
