"""Validate damped GLM cleaning using a controlled 2D divergence perturbation."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.mhd2d.initial_conditions import divergence_perturbation
from mhd_solver.mhd2d.solver import MHD2DConfig, magnetic_divergence, solve


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument("--final-time", type=float, default=1.0)
    parser.add_argument("--cleaning-speed", type=float, default=2.0)
    parser.add_argument("--damping-rate", type=float, default=2.0)
    parser.add_argument(
        "--figure", type=Path, default=Path("figures/divergence_cleaning.png")
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benchmarks/convergence/divergence_cleaning.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    common = dict(
        nx=args.resolution,
        ny=args.resolution,
        final_time=args.final_time,
        cleaning_speed=args.cleaning_speed,
        order=2,
        boundary="periodic",
    )
    undamped = solve(
        MHD2DConfig(**common, cleaning_damping_rate=0.0),
        divergence_perturbation,
    )
    damped_config = MHD2DConfig(**common, cleaning_damping_rate=args.damping_rate)
    damped = solve(damped_config, divergence_perturbation)
    dx = (damped_config.x_max - damped_config.x_min) / damped_config.nx
    dy = (damped_config.y_max - damped_config.y_min) / damped_config.ny
    final_divergence = magnetic_divergence(
        damped.primitive, dx, dy, damped_config.boundary
    )

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    axes[0].semilogy(
        undamped.divergence_times,
        undamped.divergence_l2,
        color="tab:orange",
        label="Hyperbolic only",
    )
    axes[0].semilogy(
        damped.divergence_times,
        damped.divergence_l2,
        color="tab:blue",
        label=f"Damped, rate = {args.damping_rate:g}",
    )
    axes[0].set_xlabel("Time")
    axes[0].set_ylabel(r"$\|\nabla\cdot\mathbf{B}\|_2$")
    axes[0].set_title("GLM divergence cleaning")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    image = axes[1].imshow(
        final_divergence,
        origin="lower",
        extent=(0.0, 1.0, 0.0, 1.0),
        cmap="RdBu_r",
        aspect="equal",
    )
    axes[1].set_xlabel("Position x")
    axes[1].set_ylabel("Position y")
    axes[1].set_title(r"Final $\nabla\cdot\mathbf{B}$")
    figure.colorbar(image, ax=axes[1], label=r"$\nabla\cdot\mathbf{B}$")
    figure.tight_layout()
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.figure, dpi=200, bbox_inches="tight")
    plt.close(figure)

    energy_change = float(damped.final_integrals[4] - damped.initial_integrals[4])
    record = {
        "problem": "sinusoidal magnetic-divergence perturbation",
        "configuration": asdict(damped_config),
        "initial_divergence_L2": float(damped.divergence_l2[0]),
        "final_damped_divergence_L2": float(damped.divergence_l2[-1]),
        "final_undamped_divergence_L2": float(undamped.divergence_l2[-1]),
        "damped_reduction_factor": float(
            damped.divergence_l2[0] / damped.divergence_l2[-1]
        ),
        "damped_total_energy_change": energy_change,
        "minimum_density": float(damped.primitive[0].min()),
        "minimum_pressure": float(damped.primitive[4].min()),
        "damped_steps": damped.steps,
        "undamped_steps": undamped.steps,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(record, indent=2))
    print(f"Figure: {args.figure}")
    print(f"Diagnostics: {args.output}")


if __name__ == "__main__":
    main()
