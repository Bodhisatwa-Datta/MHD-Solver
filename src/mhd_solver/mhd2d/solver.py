"""Unsplit Cartesian finite-volume solver for 2D ideal GLM-MHD."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np

from mhd_solver.common.reconstruction import minmod
from mhd_solver.mhd2d.boundaries import fill_boundaries
from mhd_solver.mhd2d.equations import (
    conserved_to_primitive,
    fast_magnetosonic_speed,
    hll_flux,
    primitive_to_conserved,
)
from mhd_solver.resistive.resistivity import resistive_rhs

Order = Literal[1, 2]
Boundary = Literal["outflow", "periodic"]


@dataclass(frozen=True)
class MHD2DConfig:
    nx: int = 64
    ny: int = 64
    x_min: float = 0.0
    x_max: float = 1.0
    y_min: float = 0.0
    y_max: float = 1.0
    final_time: float = 1.0
    gamma: float = 5.0 / 3.0
    cfl: float = 0.35
    order: Order = 2
    boundary: Boundary = "periodic"
    cleaning_speed: float = 2.0
    cleaning_damping_rate: float = 2.0
    resistivity: float = 0.0
    diffusion_cfl: float = 0.8
    max_steps: int = 1_000_000

    def __post_init__(self) -> None:
        if self.nx < 4 or self.ny < 4 or self.x_max <= self.x_min or self.y_max <= self.y_min:
            raise ValueError("invalid 2D grid dimensions or extent")
        if self.final_time < 0.0 or self.gamma <= 1.0 or not 0.0 < self.cfl <= 1.0:
            raise ValueError("require final_time >= 0, gamma > 1, and 0 < CFL <= 1")
        if self.order not in (1, 2) or self.boundary not in ("outflow", "periodic"):
            raise ValueError("invalid order or boundary type")
        if self.cleaning_speed <= 0.0 or self.cleaning_damping_rate < 0.0:
            raise ValueError("cleaning speed must be positive and damping non-negative")
        if self.resistivity < 0.0 or not 0.0 < self.diffusion_cfl <= 1.0:
            raise ValueError("resistivity must be non-negative and 0 < diffusion_cfl <= 1")


@dataclass(frozen=True)
class MHD2DResult:
    x: np.ndarray
    y: np.ndarray
    conserved: np.ndarray
    primitive: np.ndarray
    time: float
    steps: int
    initial_integrals: np.ndarray
    final_integrals: np.ndarray
    divergence_times: np.ndarray
    divergence_l2: np.ndarray


def cell_centres(config: MHD2DConfig) -> tuple[np.ndarray, np.ndarray]:
    dx = (config.x_max - config.x_min) / config.nx
    dy = (config.y_max - config.y_min) / config.ny
    x = config.x_min + (np.arange(config.nx) + 0.5) * dx
    y = config.y_min + (np.arange(config.ny) + 0.5) * dy
    return x, y


def conserved_integrals(conserved: np.ndarray, dx: float, dy: float) -> np.ndarray:
    return np.sum(conserved, axis=(1, 2)) * dx * dy


def magnetic_divergence(
    primitive: np.ndarray, dx: float, dy: float, boundary: Boundary
) -> np.ndarray:
    """Return a cell-centred second-order estimate of ``div(B)``."""
    bx, by = primitive[5], primitive[6]
    if boundary == "periodic":
        derivative_x = (np.roll(bx, -1, axis=1) - np.roll(bx, 1, axis=1)) / (2.0 * dx)
        derivative_y = (np.roll(by, -1, axis=0) - np.roll(by, 1, axis=0)) / (2.0 * dy)
    else:
        derivative_x = np.gradient(bx, dx, axis=1, edge_order=2)
        derivative_y = np.gradient(by, dy, axis=0, edge_order=2)
    return derivative_x + derivative_y


def _reconstruct(primitive: np.ndarray, direction: int) -> tuple[np.ndarray, np.ndarray]:
    slopes = np.zeros_like(primitive)
    if direction == 0:
        slopes[:, :, 1:-1] = minmod(
            primitive[:, :, 1:-1] - primitive[:, :, :-2],
            primitive[:, :, 2:] - primitive[:, :, 1:-1],
        )
        left_cells = primitive[:, 2:-2, 1:-2]
        right_cells = primitive[:, 2:-2, 2:-1]
        left = left_cells + 0.5 * slopes[:, 2:-2, 1:-2]
        right = right_cells - 0.5 * slopes[:, 2:-2, 2:-1]
    else:
        slopes[:, 1:-1, :] = minmod(
            primitive[:, 1:-1, :] - primitive[:, :-2, :],
            primitive[:, 2:, :] - primitive[:, 1:-1, :],
        )
        left_cells = primitive[:, 1:-2, 2:-2]
        right_cells = primitive[:, 2:-1, 2:-2]
        left = left_cells + 0.5 * slopes[:, 1:-2, 2:-2]
        right = right_cells - 0.5 * slopes[:, 2:-1, 2:-2]
    bad_left = (left[0] <= 0.0) | (left[4] <= 0.0)
    bad_right = (right[0] <= 0.0) | (right[4] <= 0.0)
    left[:, bad_left] = left_cells[:, bad_left]
    right[:, bad_right] = right_cells[:, bad_right]
    return left, right


def _spatial_operator(conserved: np.ndarray, dx: float, dy: float, config: MHD2DConfig) -> np.ndarray:
    extended = fill_boundaries(conserved, config.boundary)
    primitive = conserved_to_primitive(extended, config.gamma, config.cleaning_speed)
    if config.order == 1:
        x_left, x_right = primitive[:, 2:-2, 1:-2], primitive[:, 2:-2, 2:-1]
        y_left, y_right = primitive[:, 1:-2, 2:-2], primitive[:, 2:-1, 2:-2]
    else:
        x_left, x_right = _reconstruct(primitive, 0)
        y_left, y_right = _reconstruct(primitive, 1)
    flux_x = hll_flux(x_left, x_right, config.gamma, config.cleaning_speed, 0)
    flux_y = hll_flux(y_left, y_right, config.gamma, config.cleaning_speed, 1)
    ideal_rhs = -(
        (flux_x[:, :, 1:] - flux_x[:, :, :-1]) / dx
        + (flux_y[:, 1:, :] - flux_y[:, :-1, :]) / dy
    )
    return ideal_rhs + resistive_rhs(
        conserved_to_primitive(conserved, config.gamma, config.cleaning_speed),
        config.resistivity,
        dx,
        dy,
        config.boundary,
    )


def _stable_timestep(conserved: np.ndarray, dx: float, dy: float, config: MHD2DConfig) -> float:
    primitive = conserved_to_primitive(conserved, config.gamma, config.cleaning_speed)
    speed_x = np.max(
        np.abs(primitive[1])
        + np.maximum(fast_magnetosonic_speed(primitive, config.gamma, 0), config.cleaning_speed)
    )
    speed_y = np.max(
        np.abs(primitive[2])
        + np.maximum(fast_magnetosonic_speed(primitive, config.gamma, 1), config.cleaning_speed)
    )
    hyperbolic_dt = config.cfl / (speed_x / dx + speed_y / dy)
    if config.resistivity == 0.0:
        return hyperbolic_dt
    diffusion_limit = 1.0 / (
        2.0 * config.resistivity * (1.0 / dx**2 + 1.0 / dy**2)
    )
    return min(hyperbolic_dt, config.diffusion_cfl * diffusion_limit)


def _damp_cleaning_field(conserved: np.ndarray, dt: float, config: MHD2DConfig) -> np.ndarray:
    if config.cleaning_damping_rate == 0.0 or dt == 0.0:
        return conserved
    damped = conserved.copy()
    old_psi = conserved[8]
    new_psi = old_psi * np.exp(-config.cleaning_damping_rate * dt)
    damped[4] -= 0.5 * (old_psi**2 - new_psi**2) / config.cleaning_speed**2
    damped[8] = new_psi
    return damped


def solve(
    config: MHD2DConfig,
    initial_condition: Callable[[np.ndarray, np.ndarray], np.ndarray],
) -> MHD2DResult:
    """Advance a two-dimensional GLM-MHD initial-value problem."""
    x, y = cell_centres(config)
    dx = (config.x_max - config.x_min) / config.nx
    dy = (config.y_max - config.y_min) / config.ny
    conserved = primitive_to_conserved(
        initial_condition(x, y), config.gamma, config.cleaning_speed
    )
    initial_integrals = conserved_integrals(conserved, dx, dy)
    time = 0.0
    primitive = conserved_to_primitive(conserved, config.gamma, config.cleaning_speed)
    divergence_times = [time]
    divergence_l2 = [float(np.sqrt(np.mean(magnetic_divergence(primitive, dx, dy, config.boundary) ** 2)))]

    for step in range(1, config.max_steps + 1):
        if time >= config.final_time:
            return MHD2DResult(
                x, y, conserved, primitive, time, step - 1,
                initial_integrals, conserved_integrals(conserved, dx, dy),
                np.asarray(divergence_times), np.asarray(divergence_l2),
            )
        dt = min(_stable_timestep(conserved, dx, dy, config), config.final_time - time)
        conserved = _damp_cleaning_field(conserved, 0.5 * dt, config)
        first_stage = conserved + dt * _spatial_operator(conserved, dx, dy, config)
        conserved_to_primitive(first_stage, config.gamma, config.cleaning_speed)
        if config.order == 1:
            conserved = first_stage
        else:
            conserved = 0.5 * (
                conserved + first_stage + dt * _spatial_operator(first_stage, dx, dy, config)
            )
        conserved = _damp_cleaning_field(conserved, 0.5 * dt, config)
        primitive = conserved_to_primitive(conserved, config.gamma, config.cleaning_speed)
        time += dt
        divergence_times.append(time)
        divergence_l2.append(
            float(np.sqrt(np.mean(magnetic_divergence(primitive, dx, dy, config.boundary) ** 2)))
        )
    raise RuntimeError(f"maximum step count ({config.max_steps}) reached at t={time:.6e}")
