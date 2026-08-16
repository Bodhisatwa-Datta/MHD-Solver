"""Measure convergence of explicit resistive magnetic diffusion."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd2d.diagnostics import magnetic_energy
from mhd_solver.mhd2d.solver import MHD2DConfig, solve
from mhd_solver.resistive.initial_conditions import magnetic_diffusion_wave


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolutions", type=int, nargs="+", default=[16, 32, 64, 128])
    parser.add_argument("--final-time", type=float, default=0.05)
    parser.add_argument("--resistivity", type=float, default=0.1)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/magnetic_diffusion_convergence.png")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/convergence/magnetic_diffusion_convergence.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    resolutions = sorted(set(args.resolutions))
    if len(resolutions) < 3 or resolutions[0] < 8:
        raise ValueError("provide at least three resolutions of eight or more cells")
    errors: list[float] = []
    steps: list[int] = []
    total_energy_changes: list[float] = []
    magnetic_energies: list[float] = []
    finest = None
    finest_config = None
    for cells in resolutions:
        config = MHD2DConfig(
            nx=cells,
            ny=cells,
            final_time=args.final_time,
            gamma=5.0 / 3.0,
            cfl=0.35,
            order=2,
            boundary="periodic",
            cleaning_speed=2.0,
            cleaning_damping_rate=0.0,
            resistivity=args.resistivity,
            diffusion_cfl=0.8,
        )
        result = solve(config, magnetic_diffusion_wave)
        exact = magnetic_diffusion_wave(
            result.x,
            result.y,
            time=result.time,
            resistivity=args.resistivity,
        )
        errors.append(float(np.mean(np.abs(result.primitive[6] - exact[6]))))
        steps.append(result.steps)
        total_energy_changes.append(
            float(result.final_integrals[4] - result.initial_integrals[4])
        )
        spacing = 1.0 / cells
        magnetic_energies.append(magnetic_energy(result.primitive, spacing, spacing))
        finest, finest_config = result, config

    orders: list[float | None] = [None]
    for coarse_error, fine_error, coarse_n, fine_n in zip(
        errors[:-1], errors[1:], resolutions[:-1], resolutions[1:]
    ):
        orders.append(
            float(np.log(coarse_error / fine_error) / np.log(fine_n / coarse_n))
        )
    assert finest is not None and finest_config is not None
    exact = magnetic_diffusion_wave(
        finest.x,
        finest.y,
        time=finest.time,
        resistivity=args.resistivity,
    )
    initial_magnetic_energy = 0.25 * (1.0e-3) ** 2
    exact_magnetic_energy = initial_magnetic_energy * np.exp(
        -2.0 * args.resistivity * (2.0 * np.pi) ** 2 * args.final_time
    )
    record = {
        "problem": "sinusoidal magnetic diffusion",
        "configuration_at_finest_resolution": asdict(finest_config),
        "wave_amplitude": 1.0e-3,
        "resolutions": resolutions,
        "field_y_L1_errors": errors,
        "experimental_orders": orders,
        "steps": steps,
        "total_energy_changes": total_energy_changes,
        "final_magnetic_energies": magnetic_energies,
        "exact_final_magnetic_energy": float(exact_magnetic_energy),
        "finest_pressure_range": [
            float(finest.primitive[4].min()),
            float(finest.primitive[4].max()),
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    mid_y = finest_config.ny // 2
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    axes[0].plot(finest.x, exact[6, mid_y], "k-", linewidth=1.5, label="Analytical")
    axes[0].plot(
        finest.x,
        finest.primitive[6, mid_y],
        color="tab:blue",
        linewidth=1.1,
        label=f"Numerical, N = {resolutions[-1]}",
    )
    axes[0].set_xlabel("Position x")
    axes[0].set_ylabel(r"Magnetic field $B_y$")
    axes[0].set_title(f"Resistive decay at t = {args.final_time:g}")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].loglog(resolutions, errors, "o-", color="tab:blue", label="Measured")
    reference = errors[0] * (resolutions[0] / np.asarray(resolutions)) ** 2
    axes[1].loglog(resolutions, reference, "k--", label=r"$N^{-2}$ reference")
    axes[1].set_xlabel("Cells per direction N")
    axes[1].set_ylabel(r"$B_y$ mean absolute error")
    axes[1].set_title("Magnetic-diffusion convergence")
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
