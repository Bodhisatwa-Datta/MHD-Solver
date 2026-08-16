"""Measure convergence of a circularly polarized 1D Alfvén wave."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd1d.initial_conditions import circularly_polarized_alfven_wave
from mhd_solver.mhd1d.solver import MHD1DConfig, MHD1DResult, solve


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[50, 100, 200, 400])
    parser.add_argument("--final-time", type=float, default=1.0)
    parser.add_argument("--cfl", type=float, default=0.4)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/alfven_wave_convergence.png")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/convergence/alfven_wave_convergence.json"),
    )
    return parser.parse_args()


def experimental_orders(errors: list[float], resolutions: list[int]) -> list[float | None]:
    orders: list[float | None] = [None]
    for coarse_error, fine_error, coarse_n, fine_n in zip(
        errors[:-1], errors[1:], resolutions[:-1], resolutions[1:]
    ):
        orders.append(float(np.log(coarse_error / fine_error) / np.log(fine_n / coarse_n)))
    return orders


def transverse_errors(result: MHD1DResult) -> tuple[float, float]:
    exact = circularly_polarized_alfven_wave(result.x, time=result.time)
    magnetic_error = np.sqrt(
        (result.primitive[5] - exact[5]) ** 2
        + (result.primitive[6] - exact[6]) ** 2
    )
    velocity_error = np.sqrt(
        (result.primitive[2] - exact[2]) ** 2
        + (result.primitive[3] - exact[3]) ** 2
    )
    return float(np.mean(magnetic_error)), float(np.mean(velocity_error))


def main() -> None:
    args = parse_arguments()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 8:
        raise ValueError("provide at least two distinct resolutions of eight or more cells")

    magnetic_errors = {"first_order": [], "second_order": []}
    velocity_errors = {"first_order": [], "second_order": []}
    steps = {"first_order": [], "second_order": []}
    integral_changes = {"first_order": [], "second_order": []}
    finest_second_order: MHD1DResult | None = None
    for cells in resolutions:
        for order, label in ((1, "first_order"), (2, "second_order")):
            config = MHD1DConfig(
                cells=cells,
                x_min=0.0,
                x_max=1.0,
                final_time=args.final_time,
                gamma=5.0 / 3.0,
                cfl=args.cfl,
                order=order,
                boundary="periodic",
                longitudinal_field=1.0,
            )
            result = solve(config, circularly_polarized_alfven_wave)
            magnetic_error, velocity_error = transverse_errors(result)
            magnetic_errors[label].append(magnetic_error)
            velocity_errors[label].append(velocity_error)
            steps[label].append(result.steps)
            integral_changes[label].append(
                np.abs(result.final_integrals - result.initial_integrals).tolist()
            )
            if cells == resolutions[-1] and order == 2:
                finest_second_order = result

    orders = {
        label: experimental_orders(errors, resolutions)
        for label, errors in magnetic_errors.items()
    }
    reference_config = MHD1DConfig(
        cells=resolutions[-1],
        x_min=0.0,
        x_max=1.0,
        final_time=args.final_time,
        gamma=5.0 / 3.0,
        cfl=args.cfl,
        boundary="periodic",
        longitudinal_field=1.0,
    )
    record = {
        "problem": "circularly polarized Alfven wave",
        "configuration": asdict(reference_config),
        "wave": {
            "density": 1.0,
            "pressure": 0.1,
            "longitudinal_field": 1.0,
            "transverse_amplitude": 0.1,
            "wavelength": 1.0,
            "alfven_speed": 1.0,
        },
        "schemes": {
            "first_order": "piecewise constant + forward Euler + HLL",
            "second_order": "MUSCL minmod + SSP-RK2 + HLL",
        },
        "resolutions": resolutions,
        "transverse_magnetic_L1_errors": magnetic_errors,
        "transverse_velocity_L1_errors": velocity_errors,
        "magnetic_experimental_orders": orders,
        "steps": steps,
        "absolute_integral_changes": integral_changes,
        "longitudinal_divergence_error": 0.0,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    assert finest_second_order is not None
    exact = circularly_polarized_alfven_wave(
        finest_second_order.x, time=finest_second_order.time
    )
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    axes[0].plot(finest_second_order.x, exact[5], "k-", linewidth=1.4, label="Exact")
    axes[0].plot(
        finest_second_order.x,
        finest_second_order.primitive[5],
        color="tab:blue",
        linewidth=1.1,
        label=f"MUSCL + RK2, N = {resolutions[-1]}",
    )
    axes[0].set_xlabel("Position x")
    axes[0].set_ylabel(r"Transverse field $B_y$")
    axes[0].set_title(f"Alfvén wave at t = {args.final_time:g}")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].loglog(resolutions, magnetic_errors["first_order"], "o-", label="First order")
    axes[1].loglog(resolutions, magnetic_errors["second_order"], "s-", label="MUSCL + RK2")
    reference = magnetic_errors["second_order"][0] * (
        resolutions[0] / np.asarray(resolutions)
    ) ** 2
    axes[1].loglog(resolutions, reference, "k--", linewidth=1.0, label=r"$N^{-2}$")
    axes[1].set_xlabel("Number of cells N")
    axes[1].set_ylabel(r"Transverse magnetic $L_1$ error")
    axes[1].set_title("Convergence after one crossing")
    axes[1].grid(which="both", alpha=0.25)
    axes[1].legend(frameon=False)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)

    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Convergence data: {args.output}")


if __name__ == "__main__":
    main()
