"""Measure time-resolved Harris reconnection across uniform resistivities."""

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
    parser.add_argument(
        "--resistivities", type=float, nargs="+", default=[0.0, 0.0005, 0.001, 0.002]
    )
    parser.add_argument("--final-time", type=float, default=1.0)
    parser.add_argument("--samples", type=int, default=11)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/reconnection_resistivity_sweep.png")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/convergence/reconnection_resistivity_sweep.json"),
    )
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="regenerate the figure from the existing output JSON without rerunning",
    )
    return parser.parse_args()


def run_case(cells: int, resistivity: float, final_time: float, samples: int):
    output_times = tuple(float(value) for value in np.linspace(0.0, final_time, samples))
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
        output_times=output_times,
    )
    return config, solve(config, perturbed_harris_sheet)


def summarize(config: MHD2DConfig, result) -> dict[str, object]:
    dx = (config.x_max - config.x_min) / config.nx
    dy = (config.y_max - config.y_min) / config.ny
    fluxes: list[float] = []
    electric_fields: list[float] = []
    kinetic_energies: list[float] = []
    magnetic_energies: list[float] = []
    divergence_norms: list[float] = []
    dimensionless_divergences: list[float] = []
    maximum_speeds: list[float] = []
    maximum_currents: list[float] = []
    density_ranges: list[list[float]] = []
    pressure_ranges: list[list[float]] = []
    for primitive in result.primitive_snapshots:
        divergence = magnetic_divergence(primitive, dx, dy, "periodic", "outflow")
        divergence_l2 = float(np.sqrt(np.mean(divergence**2)))
        field_rms = float(np.sqrt(np.mean(np.sum(primitive[5:8] ** 2, axis=0))))
        speed = np.sqrt(np.sum(primitive[1:4] ** 2, axis=0))
        fluxes.append(reconnected_flux(primitive, result.x, result.y))
        electric_fields.append(
            reconnection_electric_field(
                primitive, result.x, result.y, config.resistivity
            )
        )
        kinetic_energies.append(kinetic_energy(primitive, dx, dy))
        magnetic_energies.append(magnetic_energy(primitive, dx, dy))
        divergence_norms.append(divergence_l2)
        dimensionless_divergences.append(min(dx, dy) * divergence_l2 / field_rms)
        maximum_speeds.append(float(speed.max()))
        maximum_currents.append(
            float(np.max(np.abs(current_density_z(primitive, dx, dy))))
        )
        density_ranges.append([float(primitive[0].min()), float(primitive[0].max())])
        pressure_ranges.append([float(primitive[4].min()), float(primitive[4].max())])

    flux_changes = np.asarray(fluxes) - fluxes[0]
    late_time = result.snapshot_times >= 0.5 * config.final_time
    late_slope = float(
        np.polyfit(result.snapshot_times[late_time], flux_changes[late_time], 1)[0]
    )
    return {
        "configuration": asdict(config),
        "steps": result.steps,
        "times": result.snapshot_times.tolist(),
        "reconnected_flux": fluxes,
        "reconnected_flux_change": flux_changes.tolist(),
        "late_time_flux_slope": late_slope,
        "x_point_electric_field": electric_fields,
        "kinetic_energy": kinetic_energies,
        "magnetic_energy": magnetic_energies,
        "divergence_L2": divergence_norms,
        "dimensionless_divergence": dimensionless_divergences,
        "maximum_speed": maximum_speeds,
        "maximum_absolute_current_density": maximum_currents,
        "density_range": density_ranges,
        "pressure_range": pressure_ranges,
        "integral_changes": (
            result.final_integrals - result.initial_integrals
        ).tolist(),
    }


def plot_record(record: dict[str, object], figure_path: Path) -> None:
    resolutions = record["resolutions"]
    resistivities = record["resistivities"]
    coarse_cells, fine_cells = resolutions[-2], resolutions[-1]
    summaries = {
        (run["resolution"], run["resistivity"]): run for run in record["runs"]
    }
    figure, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    colours = plt.cm.viridis(np.linspace(0.1, 0.9, len(resistivities)))
    for column, cells in enumerate((coarse_cells, fine_cells)):
        for resistivity, colour in zip(resistivities, colours):
            summary = summaries[cells, resistivity]
            label = fr"$\eta={resistivity:g}$"
            axes[0, column].plot(
                summary["times"],
                summary["reconnected_flux_change"],
                color=colour,
                marker="o",
                markersize=3,
                label=label,
            )
            axes[1, column].plot(
                summary["times"],
                np.abs(summary["x_point_electric_field"]),
                color=colour,
                marker="o",
                markersize=3,
                label=label,
            )
        axes[0, column].set_title(f"Reconnected flux change, N = {cells}")
        axes[0, column].set_ylabel(r"$\Delta\Psi_{rec}$")
        axes[1, column].set_title(f"X-point electric field, N = {cells}")
        axes[1, column].set_ylabel(r"$|E_z(X)|$")
        for axis in axes[:, column]:
            axis.set_xlabel("Time")
            axis.grid(alpha=0.25)
            axis.legend(frameon=False, fontsize=8)

    for cells, marker in ((coarse_cells, "o"), (fine_cells, "s")):
        final_changes = [
            summaries[cells, eta]["reconnected_flux_change"][-1]
            for eta in resistivities
        ]
        late_slopes = [
            summaries[cells, eta]["late_time_flux_slope"] for eta in resistivities
        ]
        axes[0, 2].plot(
            resistivities, final_changes, marker=marker, label=f"N = {cells}"
        )
        axes[1, 2].plot(
            resistivities, late_slopes, marker=marker, label=f"N = {cells}"
        )
    axes[0, 2].set_title("Final flux change versus resistivity")
    axes[0, 2].set_xlabel(r"Resistivity $\eta$")
    axes[0, 2].set_ylabel(r"$\Delta\Psi_{rec}(t_f)$")
    axes[1, 2].set_title("Late-time flux slope")
    axes[1, 2].set_xlabel(r"Resistivity $\eta$")
    axes[1, 2].set_ylabel(r"$d\Psi_{rec}/dt$")
    tick_labels = [
        "0" if value == 0.0 else f"{value:.0e}" for value in resistivities
    ]
    for axis in axes[:, 2]:
        axis.set_xticks(resistivities, tick_labels)
        axis.grid(alpha=0.25)
        axis.legend(frameon=False)
    figure.suptitle("Time-resolved Harris-sheet resistivity sweep")
    figure_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_arguments()
    if args.plot_only:
        record = json.loads(args.output.read_text(encoding="utf-8"))
        plot_record(record, args.figure)
        print(f"Figure: {args.figure}")
        print(f"Diagnostics: {args.output}")
        return
    resolutions = sorted(set(args.resolutions))
    resistivities = sorted(set(args.resistivities))
    if len(resolutions) < 2 or resolutions[0] < 16:
        raise ValueError("provide at least two resolutions of 16 or more cells")
    if len(resistivities) < 3 or resistivities[0] != 0.0:
        raise ValueError("include an ideal control and at least two positive resistivities")
    if args.final_time <= 0.0 or args.samples < 5:
        raise ValueError("final time must be positive and at least five samples are required")

    run_records: list[dict[str, object]] = []
    results: dict[tuple[int, float], object] = {}
    summaries: dict[tuple[int, float], dict[str, object]] = {}
    for cells in resolutions:
        for resistivity in resistivities:
            config, result = run_case(
                cells, resistivity, args.final_time, args.samples
            )
            summary = summarize(config, result)
            results[cells, resistivity] = result
            summaries[cells, resistivity] = summary
            run_records.append(
                {"resolution": cells, "resistivity": resistivity, **summary}
            )

    resolution_sensitivity: list[dict[str, float]] = []
    coarse_cells, fine_cells = resolutions[-2], resolutions[-1]
    if fine_cells % coarse_cells == 0:
        ratio = fine_cells // coarse_cells
        for resistivity in resistivities:
            coarse = results[coarse_cells, resistivity]
            fine = results[fine_cells, resistivity]
            restricted = fine.primitive.reshape(
                9, coarse_cells, ratio, coarse_cells, ratio
            ).mean(axis=(2, 4))
            resolution_sensitivity.append(
                {
                    "resistivity": resistivity,
                    "final_flux_change_difference": float(
                        summaries[fine_cells, resistivity]["reconnected_flux_change"][-1]
                        - summaries[coarse_cells, resistivity]["reconnected_flux_change"][-1]
                    ),
                    "late_time_slope_difference": float(
                        summaries[fine_cells, resistivity]["late_time_flux_slope"]
                        - summaries[coarse_cells, resistivity]["late_time_flux_slope"]
                    ),
                    "density_L1_difference": float(
                        np.mean(np.abs(coarse.primitive[0] - restricted[0]))
                    ),
                    "pressure_L1_difference": float(
                        np.mean(np.abs(coarse.primitive[4] - restricted[4]))
                    ),
                    "field_x_L1_difference": float(
                        np.mean(np.abs(coarse.primitive[5] - restricted[5]))
                    ),
                    "field_y_L1_difference": float(
                        np.mean(np.abs(coarse.primitive[6] - restricted[6]))
                    ),
                }
            )

    scaling_fits: dict[str, object] = {}
    positive_resistivities = np.asarray(resistivities[1:])
    for cells in resolutions:
        slopes = np.asarray(
            [summaries[cells, eta]["late_time_flux_slope"] for eta in resistivities[1:]]
        )
        if np.all(slopes > 0.0):
            exponent, log_prefactor = np.polyfit(
                np.log(positive_resistivities), np.log(slopes), 1
            )
            scaling_fits[str(cells)] = {
                "available": True,
                "late_time_rate_exponent": float(exponent),
                "prefactor": float(np.exp(log_prefactor)),
            }
        else:
            scaling_fits[str(cells)] = {
                "available": False,
                "reason": "not all positive-resistivity late-time slopes are positive",
            }

    record = {
        "problem": "time-resolved Harris-sheet resistivity sweep",
        "resolutions": resolutions,
        "resistivities": resistivities,
        "final_time": args.final_time,
        "samples": args.samples,
        "equilibrium": {
            "magnetic_field": 1.0,
            "sheet_half_width": 0.05,
            "background_pressure": 0.2,
            "background_density": 1.0,
        },
        "perturbation": {
            "vector_potential_amplitude": 0.01,
            "vertical_width": 0.1,
            "x_point": [0.5, 0.0],
            "o_point": [0.0, 0.0],
        },
        "runs": run_records,
        "resolution_sensitivity": resolution_sensitivity,
        "positive_resistivity_scaling_fits": scaling_fits,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    plot_record(record, args.figure)
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Diagnostics: {args.output}")


if __name__ == "__main__":
    main()
