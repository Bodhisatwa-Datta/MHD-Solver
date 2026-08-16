"""Run and plot the one-dimensional Brio-Wu ideal-MHD shock tube."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd1d.equations import conserved_to_primitive
from mhd_solver.mhd1d.initial_conditions import brio_wu_shock_tube
from mhd_solver.mhd1d.solver import MHD1DConfig, solve


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=800)
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument("--order", type=int, choices=(1, 2), default=2)
    parser.add_argument(
        "--comparison-cells",
        type=int,
        default=1600,
        help="resolution for the overlaid sensitivity run (0 disables it)",
    )
    parser.add_argument("--output", type=Path, default=Path("figures/brio_wu.png"))
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    config = MHD1DConfig(
        cells=args.cells,
        final_time=args.final_time,
        cfl=args.cfl,
        order=args.order,
    )
    result = solve(config, brio_wu_shock_tube)
    comparison = None
    if args.comparison_cells:
        if args.comparison_cells <= args.cells:
            raise ValueError("comparison-cells must exceed cells, or be zero")
        comparison_config = MHD1DConfig(
            cells=args.comparison_cells,
            final_time=args.final_time,
            cfl=args.cfl,
            order=args.order,
        )
        comparison = solve(comparison_config, brio_wu_shock_tube)
    density, velocity_x, velocity_y, velocity_z, pressure, field_y, field_z = (
        result.primitive
    )
    quantities = (
        (0, r"Density $\rho$", density),
        (4, r"Gas pressure $p$", pressure),
        (1, r"Longitudinal velocity $v_x$", velocity_x),
        (2, r"Transverse velocity $v_y$", velocity_y),
        (5, r"Transverse field $B_y$", field_y),
        (6, r"Out-of-plane field $B_z$", field_z),
    )
    figure, axes = plt.subplots(3, 2, figsize=(10, 9), sharex=True)
    for axis, (primitive_index, label, values) in zip(axes.flat, quantities):
        if comparison is not None:
            axis.plot(
                comparison.x,
                comparison.primitive[primitive_index],
                "k--",
                linewidth=0.9,
                label=f"N = {args.comparison_cells}",
            )
        axis.plot(
            result.x,
            values,
            color="tab:blue",
            linewidth=1.1,
            label=f"N = {args.cells}",
        )
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False)
    axes[-1, 0].set_xlabel("Position x")
    axes[-1, 1].set_xlabel("Position x")
    figure.suptitle(
        f"Brio-Wu shock tube at t = {result.time:.3f} "
        f"(N = {config.cells}, HLL, order {config.order})"
    )
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=200, bbox_inches="tight")
    plt.close(figure)

    dx = (config.x_max - config.x_min) / config.cells
    magnetic_energy = 0.5 * np.sum(
        config.longitudinal_field**2 + field_y**2 + field_z**2
    ) * dx
    kinetic_energy = 0.5 * np.sum(
        density * (velocity_x**2 + velocity_y**2 + velocity_z**2)
    ) * dx
    record = {
        "problem": "Brio-Wu shock tube",
        "configuration": asdict(config),
        "boundary_conditions": "transmissive/outflow",
        "riemann_solver": "HLL",
        "steps": result.steps,
        "integral_changes": (
            result.final_integrals - result.initial_integrals
        ).tolist(),
        "final_kinetic_energy": float(kinetic_energy),
        "final_magnetic_energy": float(magnetic_energy),
        "primitive_extrema": {
            "density": [float(density.min()), float(density.max())],
            "velocity_x": [float(velocity_x.min()), float(velocity_x.max())],
            "velocity_y": [float(velocity_y.min()), float(velocity_y.max())],
            "velocity_z": [float(velocity_z.min()), float(velocity_z.max())],
            "pressure": [float(pressure.min()), float(pressure.max())],
            "field_y": [float(field_y.min()), float(field_y.max())],
            "field_z": [float(field_z.min()), float(field_z.max())],
        },
        "longitudinal_divergence_error": 0.0,
    }
    if comparison is not None:
        if comparison.x.size % result.x.size != 0:
            raise ValueError("comparison resolution must be an integer multiple of cells")
        ratio = comparison.x.size // result.x.size
        restricted = comparison.conserved.reshape(7, result.x.size, ratio).mean(axis=2)
        restricted_primitive = conserved_to_primitive(
            restricted, config.gamma, config.longitudinal_field
        )
        record["resolution_sensitivity"] = {
            "comparison_cells": comparison.x.size,
            "comparison_steps": comparison.steps,
            "coarse_vs_restricted_fine_L1": {
                "density": float(np.mean(np.abs(density - restricted_primitive[0]))),
                "velocity_x": float(
                    np.mean(np.abs(velocity_x - restricted_primitive[1]))
                ),
                "velocity_y": float(
                    np.mean(np.abs(velocity_y - restricted_primitive[2]))
                ),
                "pressure": float(np.mean(np.abs(pressure - restricted_primitive[4]))),
                "field_y": float(np.mean(np.abs(field_y - restricted_primitive[5]))),
            },
        }
    metadata_path = args.output.with_suffix(".json")
    metadata_path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.output}")
    print(f"Parameters and diagnostics: {metadata_path}")


if __name__ == "__main__":
    main()
