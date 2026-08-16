"""Run the 2D Orszag-Tang vortex with GLM divergence cleaning."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd2d.diagnostics import current_density_z, kinetic_energy, magnetic_energy
from mhd_solver.mhd2d.equations import conserved_to_primitive
from mhd_solver.mhd2d.initial_conditions import orszag_tang_vortex
from mhd_solver.mhd2d.solver import MHD2DConfig, magnetic_divergence, solve


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument("--comparison-resolution", type=int, default=128)
    parser.add_argument("--final-time", type=float, default=0.5)
    parser.add_argument("--cfl", type=float, default=0.3)
    parser.add_argument("--figure", type=Path, default=Path("figures/orszag_tang.png"))
    return parser.parse_args()


def run(resolution: int, args: argparse.Namespace):
    config = MHD2DConfig(
        nx=resolution,
        ny=resolution,
        final_time=args.final_time,
        cfl=args.cfl,
        order=2,
        boundary="periodic",
        cleaning_speed=2.0,
        cleaning_damping_rate=2.0,
    )
    return config, solve(config, orszag_tang_vortex)


def main() -> None:
    args = parse_arguments()
    if args.comparison_resolution <= args.resolution:
        raise ValueError("comparison resolution must exceed the primary resolution")
    if args.comparison_resolution % args.resolution != 0:
        raise ValueError("comparison resolution must be an integer multiple of resolution")
    coarse_config, coarse = run(args.resolution, args)
    fine_config, fine = run(args.comparison_resolution, args)
    primitive = fine.primitive
    dx = (fine_config.x_max - fine_config.x_min) / fine_config.nx
    dy = (fine_config.y_max - fine_config.y_min) / fine_config.ny
    density, vx, vy, _, pressure, bx, by, bz, _ = primitive
    magnetic_pressure = 0.5 * (bx**2 + by**2 + bz**2)
    current = current_density_z(primitive, dx, dy)
    divergence = magnetic_divergence(primitive, dx, dy, fine_config.boundary)
    velocity_magnitude = np.sqrt(vx**2 + vy**2)

    quantities = (
        (density, r"Density $\rho$", "viridis"),
        (pressure, r"Gas pressure $p$", "magma"),
        (magnetic_pressure, r"Magnetic pressure $B^2/2$", "cividis"),
        (np.abs(current), r"Current magnitude $|J_z|$", "inferno"),
        (velocity_magnitude, r"Velocity magnitude $|\mathbf{v}|$", "viridis"),
        (divergence, r"Magnetic divergence $\nabla\cdot\mathbf{B}$", "RdBu_r"),
    )
    figure, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for axis, (values, title, colour_map) in zip(axes.flat, quantities):
        image = axis.imshow(
            values,
            origin="lower",
            extent=(0.0, 1.0, 0.0, 1.0),
            cmap=colour_map,
            aspect="equal",
        )
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        figure.colorbar(image, ax=axis, shrink=0.82)
    stride = max(1, fine_config.nx // 20)
    axes[0, 0].quiver(
        fine.x[::stride],
        fine.y[::stride],
        vx[::stride, ::stride],
        vy[::stride, ::stride],
        color="white",
        alpha=0.65,
        scale=25,
        width=0.0025,
    )
    figure.suptitle(
        f"Orszag-Tang vortex at t = {fine.time:.2f} "
        f"(N = {fine_config.nx}, HLL + GLM)"
    )
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)

    ratio = fine_config.nx // coarse_config.nx
    restricted = fine.conserved.reshape(
        9, coarse_config.ny, ratio, coarse_config.nx, ratio
    ).mean(axis=(2, 4))
    restricted_primitive = conserved_to_primitive(
        restricted, coarse_config.gamma, coarse_config.cleaning_speed
    )
    field_rms = float(np.sqrt(np.mean(bx**2 + by**2 + bz**2)))
    record = {
        "problem": "Orszag-Tang vortex",
        "primary_configuration": asdict(coarse_config),
        "comparison_configuration": asdict(fine_config),
        "primary_steps": coarse.steps,
        "comparison_steps": fine.steps,
        "primitive_extrema_at_comparison_resolution": {
            "density": [float(density.min()), float(density.max())],
            "pressure": [float(pressure.min()), float(pressure.max())],
            "velocity_magnitude": [float(velocity_magnitude.min()), float(velocity_magnitude.max())],
            "magnetic_pressure": [float(magnetic_pressure.min()), float(magnetic_pressure.max())],
            "current_magnitude": [float(np.abs(current).min()), float(np.abs(current).max())],
        },
        "final_kinetic_energy": kinetic_energy(primitive, dx, dy),
        "final_magnetic_energy": magnetic_energy(primitive, dx, dy),
        "integral_changes": (fine.final_integrals - fine.initial_integrals).tolist(),
        "divergence": {
            "final_L2": float(fine.divergence_l2[-1]),
            "maximum_L2": float(fine.divergence_l2.max()),
            "dimensionless_dx_L2_over_B_rms": float(dx * fine.divergence_l2[-1] / field_rms),
        },
        "coarse_vs_restricted_fine_L1": {
            "density": float(np.mean(np.abs(coarse.primitive[0] - restricted_primitive[0]))),
            "pressure": float(np.mean(np.abs(coarse.primitive[4] - restricted_primitive[4]))),
            "field_x": float(np.mean(np.abs(coarse.primitive[5] - restricted_primitive[5]))),
            "field_y": float(np.mean(np.abs(coarse.primitive[6] - restricted_primitive[6]))),
        },
    }
    metadata_path = args.figure.with_suffix(".json")
    metadata_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Parameters and diagnostics: {metadata_path}")


if __name__ == "__main__":
    main()
