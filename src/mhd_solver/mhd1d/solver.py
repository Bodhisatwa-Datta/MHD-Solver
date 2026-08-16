"""Uniform-grid finite-volume solver for one-dimensional ideal MHD."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np

from mhd_solver.common.boundaries import apply_outflow_boundaries, apply_periodic_boundaries
from mhd_solver.common.reconstruction import muscl_interface_states
from mhd_solver.mhd1d.equations import (
    conserved_to_primitive,
    fast_magnetosonic_speed,
    hll_flux,
    primitive_to_conserved,
)

Order = Literal[1, 2]
Boundary = Literal["outflow", "periodic"]


@dataclass(frozen=True)
class MHD1DConfig:
    """Physical and numerical settings for a one-dimensional ideal-MHD run."""

    cells: int = 800
    x_min: float = -0.5
    x_max: float = 0.5
    final_time: float = 0.1
    gamma: float = 2.0
    cfl: float = 0.4
    order: Order = 2
    boundary: Boundary = "outflow"
    longitudinal_field: float = 0.75
    max_steps: int = 1_000_000

    def __post_init__(self) -> None:
        if self.cells < 4 or self.x_max <= self.x_min or self.final_time < 0.0:
            raise ValueError("invalid grid extent, cell count, or final time")
        if self.gamma <= 1.0 or not 0.0 < self.cfl <= 1.0 or self.order not in (1, 2):
            raise ValueError("require gamma > 1, 0 < CFL <= 1, and order 1 or 2")
        if self.boundary not in ("outflow", "periodic"):
            raise ValueError("boundary must be 'outflow' or 'periodic'")
        if not np.isfinite(self.longitudinal_field):
            raise ValueError("longitudinal magnetic field must be finite")


@dataclass(frozen=True)
class MHD1DResult:
    x: np.ndarray
    conserved: np.ndarray
    primitive: np.ndarray
    time: float
    steps: int
    initial_integrals: np.ndarray
    final_integrals: np.ndarray


def cell_centres(config: MHD1DConfig) -> np.ndarray:
    dx = (config.x_max - config.x_min) / config.cells
    return config.x_min + (np.arange(config.cells) + 0.5) * dx


def conserved_integrals(conserved: np.ndarray, dx: float) -> np.ndarray:
    """Integrate the seven evolved conserved variables over physical cells."""
    return np.sum(conserved, axis=1) * dx


def _spatial_operator(conserved: np.ndarray, dx: float, config: MHD1DConfig) -> np.ndarray:
    if config.boundary == "periodic":
        extended = apply_periodic_boundaries(conserved)
    else:
        extended = apply_outflow_boundaries(conserved)
    primitive = conserved_to_primitive(
        extended, config.gamma, config.longitudinal_field
    )
    if config.order == 1:
        left, right = primitive[:, 1:-2], primitive[:, 2:-1]
    else:
        left, right = muscl_interface_states(primitive, pressure_index=4)
    interface_flux = hll_flux(
        left, right, config.gamma, config.longitudinal_field
    )
    return -(interface_flux[:, 1:] - interface_flux[:, :-1]) / dx


def _stable_timestep(conserved: np.ndarray, dx: float, config: MHD1DConfig) -> float:
    primitive = conserved_to_primitive(
        conserved, config.gamma, config.longitudinal_field
    )
    maximum_speed = np.max(
        np.abs(primitive[1])
        + fast_magnetosonic_speed(
            primitive, config.gamma, config.longitudinal_field
        )
    )
    if not np.isfinite(maximum_speed) or maximum_speed <= 0.0:
        raise RuntimeError(f"invalid maximum MHD wave speed: {maximum_speed}")
    return config.cfl * dx / maximum_speed


def solve(
    config: MHD1DConfig, initial_condition: Callable[[np.ndarray], np.ndarray]
) -> MHD1DResult:
    """Advance a one-dimensional ideal-MHD initial-value problem."""
    x = cell_centres(config)
    dx = (config.x_max - config.x_min) / config.cells
    conserved = primitive_to_conserved(
        initial_condition(x), config.gamma, config.longitudinal_field
    )
    initial_integrals = conserved_integrals(conserved, dx)
    time = 0.0

    for step in range(1, config.max_steps + 1):
        if time >= config.final_time:
            return MHD1DResult(
                x=x,
                conserved=conserved,
                primitive=conserved_to_primitive(
                    conserved, config.gamma, config.longitudinal_field
                ),
                time=time,
                steps=step - 1,
                initial_integrals=initial_integrals,
                final_integrals=conserved_integrals(conserved, dx),
            )
        dt = min(_stable_timestep(conserved, dx, config), config.final_time - time)
        first_stage = conserved + dt * _spatial_operator(conserved, dx, config)
        conserved_to_primitive(first_stage, config.gamma, config.longitudinal_field)
        if config.order == 1:
            conserved = first_stage
        else:
            conserved = 0.5 * (
                conserved
                + first_stage
                + dt * _spatial_operator(first_stage, dx, config)
            )
            conserved_to_primitive(
                conserved, config.gamma, config.longitudinal_field
            )
        time += dt

    raise RuntimeError(f"maximum step count ({config.max_steps}) reached at t={time:.6e}")
