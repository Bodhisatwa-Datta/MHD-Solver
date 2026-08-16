"""Measure effective numerical resistivity from sinusoidal-mode decay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd2d.solver import MHD2DConfig, solve
from mhd_solver.resistive.diagnostics import (
    effective_resistivity,
    sinusoidal_mode_amplitude,
)
from mhd_solver.resistive.initial_conditions import magnetic_diffusion_wave


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[16, 32, 64, 128])
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--reference-resistivity", type=float, default=0.001)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/numerical_resistivity.png")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/convergence/numerical_resistivity.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 3 or resolutions[0] < 8:
        raise ValueError("provide at least three resolutions of eight or more cells")
    if args.final_time <= 0.0 or args.reference_resistivity <= 0.0:
        raise ValueError("final time and reference resistivity must be positive")

    schemes: dict[str, dict[str, object]] = {}
    for order in (1, 2):
        amplitudes: list[float] = []
        resistivities: list[float] = []
        steps: list[int] = []
        maximum_speeds: list[float] = []
        total_energy_changes: list[float] = []
        for cells in resolutions:
            config = MHD2DConfig(
                nx=cells,
                ny=cells,
                final_time=args.final_time,
                gamma=5.0 / 3.0,
                cfl=0.3,
                order=order,
                boundary="periodic",
                cleaning_speed=2.0,
                cleaning_damping_rate=0.0,
                resistivity=0.0,
            )
            result = solve(config, magnetic_diffusion_wave)
            initial = magnetic_diffusion_wave(result.x, result.y)
            initial_amplitude = sinusoidal_mode_amplitude(initial[6], result.x)
            final_amplitude = sinusoidal_mode_amplitude(result.primitive[6], result.x)
            amplitudes.append(final_amplitude)
            resistivities.append(
                effective_resistivity(
                    initial_amplitude, final_amplitude, result.time
                )
            )
            steps.append(result.steps)
            maximum_speeds.append(
                float(np.max(np.sqrt(np.sum(result.primitive[1:4] ** 2, axis=0))))
            )
            total_energy_changes.append(
                float(result.final_integrals[4] - result.initial_integrals[4])
            )

        spacings = 1.0 / np.asarray(resolutions, dtype=float)
        fit_order, log_coefficient = np.polyfit(
            np.log(spacings), np.log(resistivities), 1
        )
        local_orders: list[float | None] = [None]
        for coarse_eta, fine_eta, coarse_n, fine_n in zip(
            resistivities[:-1], resistivities[1:], resolutions[:-1], resolutions[1:]
        ):
            local_orders.append(
                float(np.log(coarse_eta / fine_eta) / np.log(fine_n / coarse_n))
            )
        schemes[f"order_{order}"] = {
            "final_mode_amplitudes": amplitudes,
            "effective_resistivities": resistivities,
            "effective_resistivity_over_reference": [
                value / args.reference_resistivity for value in resistivities
            ],
            "local_scaling_orders": local_orders,
            "power_law_fit_order": float(fit_order),
            "power_law_fit_coefficient": float(np.exp(log_coefficient)),
            "steps": steps,
            "maximum_speeds": maximum_speeds,
            "total_energy_changes": total_energy_changes,
        }

    reference_total_resistivities: list[float] = []
    reference_excess_resistivities: list[float] = []
    reference_amplitudes: list[float] = []
    for cells, numerical_resistivity in zip(
        resolutions, schemes["order_2"]["effective_resistivities"]
    ):
        config = MHD2DConfig(
            nx=cells,
            ny=cells,
            final_time=args.final_time,
            gamma=5.0 / 3.0,
            cfl=0.3,
            order=2,
            boundary="periodic",
            cleaning_speed=2.0,
            cleaning_damping_rate=0.0,
            resistivity=args.reference_resistivity,
        )
        result = solve(config, magnetic_diffusion_wave)
        initial = magnetic_diffusion_wave(result.x, result.y)
        initial_amplitude = sinusoidal_mode_amplitude(initial[6], result.x)
        final_amplitude = sinusoidal_mode_amplitude(result.primitive[6], result.x)
        total_resistivity = effective_resistivity(
            initial_amplitude, final_amplitude, result.time
        )
        reference_amplitudes.append(final_amplitude)
        reference_total_resistivities.append(total_resistivity)
        reference_excess_resistivities.append(
            total_resistivity - args.reference_resistivity
        )

    record = {
        "problem": "effective numerical resistivity from weak Fourier-mode decay",
        "explicit_resistivity": 0.0,
        "reference_reconnection_resistivity": args.reference_resistivity,
        "initial_mode_amplitude": 1.0e-3,
        "wavenumber": float(2.0 * np.pi),
        "final_time": args.final_time,
        "cfl": 0.3,
        "cleaning_speed": 2.0,
        "resolutions": resolutions,
        "schemes": schemes,
        "second_order_explicit_reference_case": {
            "final_mode_amplitudes": reference_amplitudes,
            "inferred_total_resistivities": reference_total_resistivities,
            "inferred_excess_over_explicit": reference_excess_resistivities,
            "excess_over_zero_eta_numerical_resistivity": [
                excess / numerical
                for excess, numerical in zip(
                    reference_excess_resistivities,
                    schemes["order_2"]["effective_resistivities"],
                )
            ],
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for order, colour in ((1, "tab:orange"), (2, "tab:blue")):
        values = schemes[f"order_{order}"]
        label = "First order" if order == 1 else "MUSCL + RK2"
        axes[0].loglog(
            resolutions,
            values["effective_resistivities"],
            "o-",
            color=colour,
            label=label,
        )
        axes[1].plot(
            resolutions,
            np.asarray(values["final_mode_amplitudes"]) / 1.0e-3,
            "o-",
            color=colour,
            label=label,
        )
    axes[0].axhline(
        args.reference_resistivity,
        color="k",
        linestyle="--",
        label=fr"Reconnection $\eta={args.reference_resistivity:g}$",
    )
    axes[0].loglog(
        resolutions,
        reference_total_resistivities,
        "s--",
        color="tab:purple",
        label="MUSCL + explicit reference",
    )
    axes[0].set_xlabel("Cells per direction N")
    axes[0].set_ylabel(r"Effective numerical resistivity $\eta_{num}$")
    axes[0].set_title("Numerical-resistivity scaling")
    axes[0].grid(which="both", alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].set_xlabel("Cells per direction N")
    axes[1].set_ylabel(r"Retained Fourier amplitude $A(t)/A_0$")
    axes[1].set_title(f"Mode retention at t = {args.final_time:g}")
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Diagnostics: {args.output}")


if __name__ == "__main__":
    main()
