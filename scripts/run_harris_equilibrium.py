"""Validate the stationary, pressure-balanced Harris current sheet."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd2d.diagnostics import current_density_z
from mhd_solver.mhd2d.solver import MHD2DConfig, magnetic_divergence, solve
from mhd_solver.reconnection.initial_conditions import harris_current_sheet


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[64, 128])
    parser.add_argument("--final-time", type=float, default=0.1)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/harris_equilibrium.png")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("figures/harris_equilibrium.json")
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 16:
        raise ValueError("provide at least two resolutions of 16 or more cells")

    runs: list[dict[str, object]] = []
    finest = None
    finest_initial = None
    finest_config = None
    for cells in resolutions:
        config = MHD2DConfig(
            nx=cells,
            ny=cells,
            y_min=-0.5,
            y_max=0.5,
            final_time=args.final_time,
            gamma=5.0 / 3.0,
            cfl=0.3,
            order=2,
            boundary_x="periodic",
            boundary_y="outflow",
            cleaning_speed=2.0,
            cleaning_damping_rate=2.0,
            resistivity=0.0,
        )
        result = solve(config, harris_current_sheet)
        initial = harris_current_sheet(result.x, result.y)
        dx = (config.x_max - config.x_min) / config.nx
        dy = (config.y_max - config.y_min) / config.ny
        speed = np.sqrt(np.sum(result.primitive[1:4] ** 2, axis=0))
        total_pressure = result.primitive[4] + 0.5 * np.sum(
            result.primitive[5:8] ** 2, axis=0
        )
        divergence = magnetic_divergence(
            result.primitive,
            dx,
            dy,
            config.resolved_boundary_x,
            config.resolved_boundary_y,
        )
        field_rms = float(np.sqrt(np.mean(np.sum(result.primitive[5:8] ** 2, axis=0))))
        current = current_density_z(result.primitive, dx, dy, periodic=False)
        runs.append(
            {
                "resolution": cells,
                "steps": result.steps,
                "maximum_speed": float(speed.max()),
                "density_L1_drift": float(np.mean(np.abs(result.primitive[0] - initial[0]))),
                "pressure_L1_drift": float(np.mean(np.abs(result.primitive[4] - initial[4]))),
                "field_x_L1_drift": float(np.mean(np.abs(result.primitive[5] - initial[5]))),
                "maximum_total_pressure_error": float(np.max(np.abs(total_pressure - 0.7))),
                "maximum_absolute_current_density": float(np.max(np.abs(current))),
                "final_divergence_L2": float(np.sqrt(np.mean(divergence**2))),
                "dimensionless_divergence": float(
                    min(dx, dy) * np.sqrt(np.mean(divergence**2)) / field_rms
                ),
                "density_range": [
                    float(result.primitive[0].min()),
                    float(result.primitive[0].max()),
                ],
                "pressure_range": [
                    float(result.primitive[4].min()),
                    float(result.primitive[4].max()),
                ],
                "integral_changes": (
                    result.final_integrals - result.initial_integrals
                ).tolist(),
            }
        )
        finest, finest_initial, finest_config = result, initial, config

    assert finest is not None and finest_initial is not None and finest_config is not None
    record = {
        "problem": "unperturbed ideal Harris current-sheet equilibrium",
        "configuration_at_finest_resolution": asdict(finest_config),
        "parameters": {
            "magnetic_field": 1.0,
            "half_width": 0.05,
            "background_pressure": 0.2,
            "background_density": 1.0,
            "constant_temperature": 0.2,
            "equilibrium_total_pressure": 0.7,
        },
        "runs": runs,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    dx = (finest_config.x_max - finest_config.x_min) / finest_config.nx
    dy = (finest_config.y_max - finest_config.y_min) / finest_config.ny
    middle = finest_config.nx // 2
    speed = np.sqrt(np.sum(finest.primitive[1:4] ** 2, axis=0))
    total_pressure = finest.primitive[4] + 0.5 * np.sum(
        finest.primitive[5:8] ** 2, axis=0
    )
    current = current_density_z(finest.primitive, dx, dy, periodic=False)
    divergence = magnetic_divergence(
        finest.primitive, dx, dy, "periodic", "outflow"
    )

    figure, axes = plt.subplots(2, 3, figsize=(13, 7.6))
    axes[0, 0].plot(finest.y, finest_initial[5, :, middle], "k--", label="Initial")
    axes[0, 0].plot(finest.y, finest.primitive[5, :, middle], label="Final")
    axes[0, 0].set_title(r"Reversing field $B_x(y)$")
    axes[0, 0].set_xlabel("y")
    axes[0, 0].legend(frameon=False)
    axes[0, 1].plot(finest.y, total_pressure[:, middle], color="tab:purple")
    axes[0, 1].axhline(0.7, color="k", linestyle="--", linewidth=1)
    axes[0, 1].set_title("Final total pressure")
    axes[0, 1].set_xlabel("y")
    axes[0, 2].plot(finest.y, current[:, middle], color="tab:red")
    axes[0, 2].set_title(r"Current density $J_z$")
    axes[0, 2].set_xlabel("y")

    extent = [finest_config.x_min, finest_config.x_max, finest_config.y_min, finest_config.y_max]
    images = (
        axes[1, 0].imshow(finest.primitive[0], origin="lower", extent=extent, aspect="auto"),
        axes[1, 1].imshow(speed, origin="lower", extent=extent, aspect="auto"),
        axes[1, 2].imshow(divergence, origin="lower", extent=extent, aspect="auto", cmap="RdBu_r"),
    )
    for axis, image, title in zip(
        axes[1], images, ("Density", "Velocity magnitude", r"$\nabla\cdot\mathbf{B}$")
    ):
        axis.set_title(title)
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        figure.colorbar(image, ax=axis, shrink=0.82)
    for axis in axes[0]:
        axis.grid(alpha=0.25)
        axis.set_ylabel("Value")
    figure.suptitle(
        f"Unperturbed Harris-sheet equilibrium at t = {args.final_time:g}, "
        f"N = {resolutions[-1]}"
    )
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Diagnostics: {args.output}")


if __name__ == "__main__":
    main()
