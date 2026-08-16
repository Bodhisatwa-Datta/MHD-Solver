"""Run a seeded Harris-sheet reconnection experiment with an ideal control."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd2d.diagnostics import kinetic_energy, magnetic_energy
from mhd_solver.mhd2d.solver import MHD2DConfig, magnetic_divergence, solve
from mhd_solver.reconnection.diagnostics import (
    current_density_z,
    reconnected_flux,
    reconnection_electric_field,
)
from mhd_solver.reconnection.initial_conditions import perturbed_harris_sheet


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[64, 128])
    parser.add_argument("--final-time", type=float, default=1.0)
    parser.add_argument("--resistivity", type=float, default=0.001)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/harris_reconnection.png")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("figures/harris_reconnection.json")
    )
    return parser.parse_args()


def run(cells: int, final_time: float, resistivity: float):
    config = MHD2DConfig(
        nx=cells,
        ny=cells,
        y_min=-0.5,
        y_max=0.5,
        final_time=final_time,
        gamma=5.0 / 3.0,
        cfl=0.3,
        order=2,
        boundary_x="periodic",
        boundary_y="outflow",
        cleaning_speed=2.0,
        cleaning_damping_rate=2.0,
        resistivity=resistivity,
        diffusion_cfl=0.8,
    )
    return config, solve(config, perturbed_harris_sheet)


def summarize(config: MHD2DConfig, result) -> dict[str, object]:
    dx = (config.x_max - config.x_min) / config.nx
    dy = (config.y_max - config.y_min) / config.ny
    initial = perturbed_harris_sheet(result.x, result.y)
    initial_flux = reconnected_flux(initial, result.x, result.y)
    final_flux = reconnected_flux(result.primitive, result.x, result.y)
    speed = np.sqrt(np.sum(result.primitive[1:4] ** 2, axis=0))
    divergence = magnetic_divergence(
        result.primitive, dx, dy, "periodic", "outflow"
    )
    field_rms = float(
        np.sqrt(np.mean(np.sum(result.primitive[5:8] ** 2, axis=0)))
    )
    return {
        "resistivity": config.resistivity,
        "steps": result.steps,
        "initial_reconnected_flux": initial_flux,
        "final_reconnected_flux": final_flux,
        "flux_change": final_flux - initial_flux,
        "average_flux_change_rate": (final_flux - initial_flux) / result.time,
        "final_x_point_electric_field": reconnection_electric_field(
            result.primitive, result.x, result.y, config.resistivity
        ),
        "initial_kinetic_energy": kinetic_energy(initial, dx, dy),
        "final_kinetic_energy": kinetic_energy(result.primitive, dx, dy),
        "initial_magnetic_energy": magnetic_energy(initial, dx, dy),
        "final_magnetic_energy": magnetic_energy(result.primitive, dx, dy),
        "maximum_speed": float(speed.max()),
        "maximum_absolute_current_density": float(
            np.max(np.abs(current_density_z(result.primitive, dx, dy)))
        ),
        "density_range": [
            float(result.primitive[0].min()),
            float(result.primitive[0].max()),
        ],
        "pressure_range": [
            float(result.primitive[4].min()),
            float(result.primitive[4].max()),
        ],
        "final_divergence_L2": float(np.sqrt(np.mean(divergence**2))),
        "dimensionless_divergence": float(
            min(dx, dy) * np.sqrt(np.mean(divergence**2)) / field_rms
        ),
        "integral_changes": (
            result.final_integrals - result.initial_integrals
        ).tolist(),
    }


def main() -> None:
    args = parse_arguments()
    resolutions = sorted(set(args.resolutions))
    if not resolutions or resolutions[0] < 16:
        raise ValueError("resolutions must contain grids of at least 16 cells")
    if args.final_time <= 0.0 or args.resistivity <= 0.0:
        raise ValueError("final time and resistivity must be positive")

    records: list[dict[str, object]] = []
    results: dict[tuple[int, float], tuple[MHD2DConfig, object]] = {}
    for cells in resolutions:
        cases: dict[str, object] = {"resolution": cells}
        for label, resistivity in (("ideal_control", 0.0), ("resistive", args.resistivity)):
            config, result = run(cells, args.final_time, resistivity)
            results[cells, resistivity] = (config, result)
            cases[label] = summarize(config, result)
        ideal_flux = cases["ideal_control"]["final_reconnected_flux"]
        resistive_flux = cases["resistive"]["final_reconnected_flux"]
        cases["resistive_excess_final_flux"] = resistive_flux - ideal_flux
        records.append(cases)

    fine_cells = resolutions[-1]
    fine_config, fine = results[fine_cells, args.resistivity]
    _, fine_ideal = results[fine_cells, 0.0]
    dx = (fine_config.x_max - fine_config.x_min) / fine_config.nx
    dy = (fine_config.y_max - fine_config.y_min) / fine_config.ny
    current = current_density_z(fine.primitive, dx, dy)
    speed = np.sqrt(np.sum(fine.primitive[1:4] ** 2, axis=0))
    divergence = magnetic_divergence(fine.primitive, dx, dy, "periodic", "outflow")
    resolution_sensitivity = None
    if len(resolutions) >= 2:
        coarse_cells = resolutions[-2]
        if fine_cells % coarse_cells == 0:
            _, coarse = results[coarse_cells, args.resistivity]
            ratio = fine_cells // coarse_cells
            restricted_fine = fine.primitive.reshape(
                9, coarse_cells, ratio, coarse_cells, ratio
            ).mean(axis=(2, 4))
            resolution_sensitivity = {
                "coarse_resolution": coarse_cells,
                "fine_resolution": fine_cells,
                "density_L1_difference": float(
                    np.mean(np.abs(coarse.primitive[0] - restricted_fine[0]))
                ),
                "pressure_L1_difference": float(
                    np.mean(np.abs(coarse.primitive[4] - restricted_fine[4]))
                ),
                "field_x_L1_difference": float(
                    np.mean(np.abs(coarse.primitive[5] - restricted_fine[5]))
                ),
                "field_y_L1_difference": float(
                    np.mean(np.abs(coarse.primitive[6] - restricted_fine[6]))
                ),
            }
    record = {
        "problem": "seeded Harris-sheet reconnection with ideal control",
        "configuration_at_finest_resolution": asdict(fine_config),
        "equilibrium": {
            "magnetic_field": 1.0,
            "sheet_half_width": 0.05,
            "background_pressure": 0.2,
            "background_density": 1.0,
            "upstream_alfven_speed": 1.0,
        },
        "perturbation": {
            "vector_potential_amplitude": 0.01,
            "vertical_width": 0.1,
            "x_point": [0.5, 0.0],
            "o_point": [0.0, 0.0],
        },
        "nominal_lundquist_numbers": {
            "definition": "S = length * upstream_alfven_speed / resistivity",
            "based_on_sheet_half_width": 0.05 / args.resistivity,
            "based_on_sheet_half_length": 0.5 / args.resistivity,
        },
        "runs": records,
        "resistive_resolution_sensitivity": resolution_sensitivity,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    extent = [0.0, 1.0, -0.5, 0.5]
    figure, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    fields = (
        (fine.primitive[0], "Density", "viridis"),
        (current, r"Current density $J_z$", "RdBu_r"),
        (speed, r"Velocity magnitude $|\mathbf{v}|$", "magma"),
        (fine.primitive[6], r"Reconnected field $B_y$", "RdBu_r"),
        (fine.primitive[6] - fine_ideal.primitive[6], r"$B_y(\eta)-B_y(0)$", "RdBu_r"),
        (divergence, r"Magnetic divergence $\nabla\cdot\mathbf{B}$", "RdBu_r"),
    )
    for axis, (values, title, colour_map) in zip(axes.flat, fields):
        image = axis.imshow(
            values, origin="lower", extent=extent, aspect="equal", cmap=colour_map
        )
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        axis.set_title(title)
        figure.colorbar(image, ax=axis, shrink=0.82)
    stride = max(1, fine_cells // 32)
    axes[0, 0].streamplot(
        fine.x,
        fine.y,
        fine.primitive[5],
        fine.primitive[6],
        color="white",
        density=1.1,
        linewidth=0.55,
        arrowsize=0.55,
    )
    axes[0, 2].quiver(
        fine.x[::stride],
        fine.y[::stride],
        fine.primitive[1, ::stride, ::stride],
        fine.primitive[2, ::stride, ::stride],
        color="white",
        alpha=0.65,
    )
    figure.suptitle(
        f"Seeded Harris-sheet reconnection: eta = {args.resistivity:g}, "
        f"t = {args.final_time:g}, N = {fine_cells}"
    )
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Diagnostics: {args.output}")


if __name__ == "__main__":
    main()
