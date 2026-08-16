# MHD Solver

A progressively developed finite-volume solver for compressible hydrodynamics and
magnetohydrodynamics. This is a readable, reproducible computational-physics
codebase rather than a production simulation package.

## Current status: Phase 2, validated 1D ideal MHD

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

The validated hydrodynamics baseline is retained. The first 1D ideal-MHD
increment adds seven-variable state conversions, ideal-MHD fluxes, the fast
magnetosonic speed, an HLL MHD solver, the Brio-Wu shock tube, and exact
Alfvén-wave convergence. The
longitudinal field is constant in 1D, satisfying the divergence constraint in
this restricted geometry.

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

Run the smooth periodic-wave convergence study with:

```bash
python scripts/run_density_wave_convergence.py
```

It advects an exact sinusoidal entropy wave for one domain crossing at 50, 100,
200, and 400 cells. The measured first-order rates approach one and the
MUSCL–RK2 rates approach two. Results are written to
`benchmarks/convergence/density_wave_convergence.json` and
`figures/density_wave_convergence.png`.

Run the Brio-Wu MHD shock tube and doubled-resolution comparison with:

```bash
python scripts/run_brio_wu.py
```

The default calculation uses 800 cells plus a 1600-cell sensitivity run and
writes `figures/brio_wu.png` and `figures/brio_wu.json`.

Run the exact circularly polarized Alfvén-wave convergence study with:

```bash
python scripts/run_alfven_convergence.py
```

It evolves one periodic domain crossing at 50, 100, 200, and 400 cells, compares
against the exact nonlinear solution, and writes
`figures/alfven_wave_convergence.png` plus
`benchmarks/convergence/alfven_wave_convergence.json`.

## Tests

```bash
python -m pytest
```

The suite is also compatible with the Python standard library runner:

```bash
python -m unittest discover -s tests
```

The tests cover hydrodynamic and MHD state conversions, fluxes, wave speeds,
reconstruction, boundaries, uniform-state preservation, Sod and Brio-Wu
evolution, smooth-wave convergence, and conservation.

## Repository layout

```text
src/mhd_solver/common/     equation of state, reconstruction, HLL flux
src/mhd_solver/hydro1d/    Euler solver, initial data, exact Sod solution
src/mhd_solver/mhd1d/      ideal-MHD equations, solver, Brio-Wu and Alfvén data
scripts/                   reproducible benchmark entry points
tests/                     automated unit and integration tests
docs/                      equations, method, and validation record
figures/                   generated benchmark figures
benchmarks/                future convergence and performance studies
```

## Reproducibility and limitations

Settings use immutable `HydroConfig` and `MHD1DConfig` objects and are exposed by
the benchmark scripts. The uniform-grid solvers have no AMR, source terms,
non-ideal equation of state, or multidimensional support. Both outflow and
periodic boundaries are available. HLL is robust but smears contact and
Alfvénic waves. MUSCL is formally second order only in smooth regions. The 1D
formulation keeps \(B_x\) constant by construction; multidimensional divergence
control remains unimplemented.

See [the equations](docs/equations.md), [the numerical method](docs/numerical_method.md),
and [the validation record](docs/validation.md) for scientific details.

## Roadmap

1. Hydrodynamics baseline (validated)
2. 1D ideal MHD: HLL, Brio-Wu, and Alfvén convergence (validated)
3. 2D ideal MHD with divergence cleaning
4. Explicit resistive MHD validated against magnetic diffusion
5. Harris-sheet reconnection experiments
