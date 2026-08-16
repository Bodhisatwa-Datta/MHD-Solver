"""Uniform-grid finite-volume solver for the 1D Euler equations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Literal

import numpy as np

from mhd_solver.common.eos import conserved_to_primitive, primitive_to_conserved, sound_speed
from mhd_solver.common.reconstruction import muscl_interface_states
from mhd_solver.common.riemann import hll_flux

Order = Literal[1, 2]
Boundary = Literal["outflow", "periodic"]


@dataclass(frozen=True)
class HydroConfig:
    """Physical and numerical settings for a 1D Euler calculation."""

    cells: int = 400
    x_min: float = 0.0
    x_max: float = 1.0
    final_time: float = 0.2
    gamma: float = 1.4
    cfl: float = 0.8
    order: Order = 2
    boundary: Boundary = "outflow"
    max_steps: int = 1_000_000

    def __post_init__(self) -> None:
        if self.cells < 4 or self.x_max <= self.x_min or self.final_time < 0.0:
            raise ValueError("invalid grid extent, cell count, or final time")
        if self.gamma <= 1.0 or not 0.0 < self.cfl <= 1.0 or self.order not in (1, 2):
            raise ValueError("require gamma > 1, 0 < CFL <= 1, and order 1 or 2")
        if self.boundary not in ("outflow", "periodic"):
            raise ValueError("boundary must be 'outflow' or 'periodic'")


@dataclass(frozen=True)
class HydroResult:
    x: np.ndarray
    conserved: np.ndarray
    primitive: np.ndarray
    time: float
    steps: int
    initial_integrals: np.ndarray
    final_integrals: np.ndarray


def cell_centres(config: HydroConfig) -> np.ndarray:
    dx = (config.x_max - config.x_min) / config.cells
    return config.x_min + (np.arange(config.cells) + 0.5) * dx


def apply_outflow_boundaries(conserved: np.ndarray, ghost_cells: int = 2) -> np.ndarray:
    """Pad cell averages with transmissive constant-extrapolation ghosts."""
    return np.pad(conserved, ((0, 0), (ghost_cells, ghost_cells)), mode="edge")


def apply_periodic_boundaries(conserved: np.ndarray, ghost_cells: int = 2) -> np.ndarray:
    """Pad cell averages by wrapping values from the opposite domain edge."""
    return np.pad(conserved, ((0, 0), (ghost_cells, ghost_cells)), mode="wrap")


def conserved_integrals(conserved: np.ndarray, dx: float) -> np.ndarray:
    """Integrate mass, momentum, and total energy over physical cells."""
    return np.sum(conserved, axis=1) * dx


def _spatial_operator(
    conserved: np.ndarray, dx: float, gamma: float, order: Order, boundary: Boundary
) -> np.ndarray:
    if boundary == "periodic":
        extended = apply_periodic_boundaries(conserved)
    else:
        extended = apply_outflow_boundaries(conserved)
    primitive = conserved_to_primitive(extended, gamma)
    if order == 1:
        left, right = primitive[:, 1:-2], primitive[:, 2:-1]
    else:
        left, right = muscl_interface_states(primitive)
    interface_flux = hll_flux(left, right, gamma)
    return -(interface_flux[:, 1:] - interface_flux[:, :-1]) / dx


def _stable_timestep(conserved: np.ndarray, dx: float, config: HydroConfig) -> float:
    primitive = conserved_to_primitive(conserved, config.gamma)
    maximum_speed = np.max(np.abs(primitive[1]) + sound_speed(primitive, config.gamma))
    if not np.isfinite(maximum_speed) or maximum_speed <= 0.0:
        raise RuntimeError(f"invalid maximum wave speed: {maximum_speed}")
    return config.cfl * dx / maximum_speed


def solve(config: HydroConfig, initial_condition: Callable[[np.ndarray], np.ndarray]) -> HydroResult:
    """Advance a 1D Euler initial-value problem to ``config.final_time``."""
    x = cell_centres(config)
    dx = (config.x_max - config.x_min) / config.cells
    conserved = primitive_to_conserved(initial_condition(x), config.gamma)
    initial_integrals = conserved_integrals(conserved, dx)
    time = 0.0

    for step in range(1, config.max_steps + 1):
        if time >= config.final_time:
            return HydroResult(
                x=x,
                conserved=conserved,
                primitive=conserved_to_primitive(conserved, config.gamma),
                time=time,
                steps=step - 1,
                initial_integrals=initial_integrals,
                final_integrals=conserved_integrals(conserved, dx),
            )
        dt = min(_stable_timestep(conserved, dx, config), config.final_time - time)
        first_stage = conserved + dt * _spatial_operator(
            conserved, dx, config.gamma, config.order, config.boundary
        )
        conserved_to_primitive(first_stage, config.gamma)
        if config.order == 1:
            conserved = first_stage
        else:
            conserved = 0.5 * (
                conserved
                + first_stage
                + dt
                * _spatial_operator(
                    first_stage, dx, config.gamma, config.order, config.boundary
                )
            )
            conserved_to_primitive(conserved, config.gamma)
        time += dt

    raise RuntimeError(f"maximum step count ({config.max_steps}) reached at t={time:.6e}")
