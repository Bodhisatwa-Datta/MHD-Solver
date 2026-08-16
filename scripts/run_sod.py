"""Run and plot first- and second-order Sod shock-tube calculations."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from mhd_solver.common.eos import primitive_to_conserved
from mhd_solver.hydro1d.exact_sod import exact_sod
from mhd_solver.hydro1d.initial_conditions import sod_shock_tube
from mhd_solver.hydro1d.solver import HydroConfig, solve


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=400)
    parser.add_argument("--final-time", type=float, default=0.2)
    parser.add_argument("--cfl", type=float, default=0.8)
    parser.add_argument("--output", type=Path, default=Path("figures/sod_comparison.png"))
    return parser.parse_args()


def main() -> None:
    args = parse_arguments()
    base = HydroConfig(cells=args.cells, final_time=args.final_time, cfl=args.cfl)
    first = solve(replace(base, order=1), sod_shock_tube)
    second = solve(replace(base, order=2), sod_shock_tube)
    exact_primitive = exact_sod(second.x, second.time, base.gamma)
    exact_energy = primitive_to_conserved(exact_primitive, base.gamma)[2]

    quantities = (
        ("Density", first.primitive[0], second.primitive[0], exact_primitive[0]),
        ("Velocity", first.primitive[1], second.primitive[1], exact_primitive[1]),
        ("Pressure", first.primitive[2], second.primitive[2], exact_primitive[2]),
        ("Total energy", first.conserved[2], second.conserved[2], exact_energy),
    )
    figure, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for axis, (label, first_values, second_values, exact_values) in zip(axes.flat, quantities):
        axis.plot(second.x, exact_values, "k-", linewidth=1.5, label="Exact")
        axis.plot(first.x, first_values, color="tab:orange", linewidth=1.0, label="First order")
        axis.plot(second.x, second_values, color="tab:blue", linewidth=1.0, label="MUSCL + RK2")
        axis.set_ylabel(label)
        axis.grid(alpha=0.2)
    axes[1, 0].set_xlabel("Position x")
    axes[1, 1].set_xlabel("Position x")
    axes[0, 0].legend(frameon=False)
    figure.suptitle(f"Sod shock tube at t = {second.time:.3f} (N = {base.cells})")
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=200, bbox_inches="tight")
    plt.close(figure)

    errors = {}
    for name, index in (("density", 0), ("velocity", 1), ("pressure", 2)):
        errors[name] = {
            "first_order_L1": float(np.mean(np.abs(first.primitive[index] - exact_primitive[index]))),
            "second_order_L1": float(np.mean(np.abs(second.primitive[index] - exact_primitive[index]))),
        }
    metadata = {
        "problem": "Sod shock tube",
        "configuration": asdict(base),
        "boundary_conditions": "transmissive/outflow",
        "riemann_solver": "HLL",
        "first_order_steps": first.steps,
        "second_order_steps": second.steps,
        "L1_errors": errors,
        "integral_changes": {
            "first_order": (first.final_integrals - first.initial_integrals).tolist(),
            "second_order": (second.final_integrals - second.initial_integrals).tolist(),
        },
    }
    metadata_path = args.output.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print(f"Figure: {args.output}")
    print(f"Parameters and diagnostics: {metadata_path}")


if __name__ == "__main__":
    main()
