# MHD Solver

A progressively developed finite-volume solver for compressible hydrodynamics and
magnetohydrodynamics. This is a readable, reproducible computational-physics
codebase rather than a production simulation package.

## Current status: Phase 1 hydrodynamics

The implemented model is the one-dimensional compressible Euler system

```text
dU/dt + dF(U)/dx = 0,
U = (rho, rho v, E)^T,
F = (rho v, rho v^2 + p, v(E + p))^T,
p = (gamma - 1) [E - rho v^2/2].
```

Implemented numerical features:

- cell-centred finite volumes on a uniform one-dimensional mesh;
- transmissive (zero-gradient) boundary conditions;
- HLL approximate Riemann flux and CFL-controlled explicit timesteps;
- first-order piecewise-constant reconstruction with forward Euler;
- second-order MUSCL reconstruction in primitive variables, a minmod limiter,
  and two-stage SSP-RK2 time integration;
- checks for non-finite values, non-positive density, and non-positive pressure;
- mass, momentum, and total-energy diagnostics;
- a reproducible Sod shock-tube benchmark with an exact Riemann solution.

No MHD equations are implemented yet. They remain deliberately out of scope until
the hydrodynamics baseline has been validated.

## Installation

Python 3.10 or newer is recommended. From the repository root:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
python -m pip install -e .
```

## Run the Sod benchmark

```bash
python scripts/run_sod.py
```

This runs 400-cell first- and second-order calculations to `t = 0.2`, prints
conservation and error diagnostics, and writes `figures/sod_comparison.png`.
Run `python scripts/run_sod.py --help` for configurable parameters.

## Tests

```bash
python -m pytest
```

The suite is also compatible with the Python standard library runner:

```bash
python -m unittest discover -s tests
```

The tests cover equation-of-state conversions, fluxes, reconstruction, outflow
boundaries, uniform-state preservation, Sod evolution, and conservation.

## Repository layout

```text
src/mhd_solver/common/     equation of state, reconstruction, HLL flux
src/mhd_solver/hydro1d/    finite-volume solver, initial data, exact Sod solution
scripts/                   reproducible benchmark entry points
tests/                     automated unit and integration tests
docs/                      equations, method, and validation record
figures/                   generated benchmark figures
benchmarks/                future convergence and performance studies
```

## Reproducibility and limitations

Settings use an immutable `HydroConfig` object and are exposed by the benchmark
script. The uniform-grid solver has no AMR, source terms, non-ideal equation of
state, or multidimensional support. HLL is robust but smears contact
discontinuities. MUSCL is formally second order only in smooth regions.

See [the equations](docs/equations.md), [the numerical method](docs/numerical_method.md),
and [the validation record](docs/validation.md) for scientific details.

## Roadmap

1. Hydrodynamics baseline (current)
2. 1D ideal MHD: HLL, Brio-Wu, and Alfvén-wave convergence
3. 2D ideal MHD with divergence cleaning
4. Explicit resistive MHD validated against magnetic diffusion
5. Harris-sheet reconnection experiments
