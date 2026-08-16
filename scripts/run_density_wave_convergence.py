"""Measure convergence for an exactly advected periodic Euler density wave."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.hydro1d.initial_conditions import smooth_density_wave
from mhd_solver.hydro1d.solver import HydroConfig, solve


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[50, 100, 200, 400])
    parser.add_argument("--final-time", type=float, default=1.0)
    parser.add_argument("--cfl", type=float, default=0.6)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/density_wave_convergence.png")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/convergence/density_wave_convergence.json"),
    )
    return parser.parse_args()


def experimental_orders(errors: list[float], resolutions: list[int]) -> list[float | None]:
    """Return pairwise convergence orders, with ``None`` at the coarsest grid."""
    orders: list[float | None] = [None]
    for coarse_error, fine_error, coarse_n, fine_n in zip(
        errors[:-1], errors[1:], resolutions[:-1], resolutions[1:]
    ):
        orders.append(float(np.log(coarse_error / fine_error) / np.log(fine_n / coarse_n)))
    return orders


def main() -> None:
    args = parse_arguments()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 2 or resolutions[0] < 4:
        raise ValueError("provide at least two distinct resolutions of four or more cells")

    errors: dict[str, list[float]] = {"first_order": [], "second_order": []}
    steps: dict[str, list[int]] = {"first_order": [], "second_order": []}
    integral_errors: dict[str, list[list[float]]] = {"first_order": [], "second_order": []}
    for cells in resolutions:
        for order, label in ((1, "first_order"), (2, "second_order")):
            config = HydroConfig(
                cells=cells,
                final_time=args.final_time,
                cfl=args.cfl,
                order=order,
                boundary="periodic",
            )
            result = solve(config, smooth_density_wave)
            exact_density = smooth_density_wave(result.x, time=result.time)[0]
            errors[label].append(float(np.mean(np.abs(result.primitive[0] - exact_density))))
            steps[label].append(result.steps)
            integral_errors[label].append(
                np.abs(result.final_integrals - result.initial_integrals).tolist()
            )

    orders = {
        label: experimental_orders(values, resolutions) for label, values in errors.items()
    }
    reference_config = HydroConfig(
        cells=resolutions[-1],
        final_time=args.final_time,
        cfl=args.cfl,
        boundary="periodic",
    )
    record = {
        "problem": "periodic smooth density (entropy) wave",
        "configuration": asdict(reference_config),
        "schemes": {
            "first_order": "piecewise constant + forward Euler",
            "second_order": "MUSCL minmod + SSP-RK2",
        },
        "resolutions": resolutions,
        "density_L1_errors": errors,
        "experimental_orders": orders,
        "steps": steps,
        "absolute_integral_changes": integral_errors,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    figure, axis = plt.subplots(figsize=(6.5, 4.5))
    axis.loglog(resolutions, errors["first_order"], "o-", label="First order")
    axis.loglog(resolutions, errors["second_order"], "s-", label="MUSCL + RK2")
    scale = errors["second_order"][0] * (resolutions[0] / np.asarray(resolutions)) ** 2
    axis.loglog(resolutions, scale, "k--", linewidth=1.0, label=r"$N^{-2}$ reference")
    axis.set_xlabel("Number of cells N")
    axis.set_ylabel(r"Density $L_1$ error")
    axis.set_title(f"Periodic density-wave convergence at t = {args.final_time:g}")
    axis.grid(which="both", alpha=0.25)
    axis.legend(frameon=False)
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)

    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Convergence data: {args.output}")


if __name__ == "__main__":
    main()
