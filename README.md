# MHD Solver

A progressively developed finite-volume solver for compressible hydrodynamics and
magnetohydrodynamics. This is a readable, reproducible computational-physics
codebase rather than a production simulation package.

## Current status: Phase 5, seeded Harris-sheet reconnection

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

The first two-dimensional increment adds an unsplit Cartesian finite-volume
solver, x/y HLL fluxes, multidimensional CFL control, two-dimensional boundary
fills, and damped generalized-Lagrange-multiplier (GLM) divergence cleaning. A
controlled divergence perturbation, nonlinear Orszag–Tang vortex, and fast
magnetic rotor have been validated. The 2D benchmarks are resolution-sensitive
and use GLM cleaning rather than constrained transport.

Uniform explicit resistivity is included in the induction and total-energy
equations and has passed a second-order sinusoidal magnetic-diffusion convergence
test. The first reconnection-phase increment adds independent x/y boundaries
and a constant-temperature, total-pressure-balanced Harris current sheet. Its
unperturbed ideal evolution has been checked at 64² and 128². A seeded,
resistive experiment now produces an X-point, magnetic island, current layer,
and reconnection outflows. An ideal control exposes substantial numerical
reconnection, so resistivity scaling and Sweet–Parker validation remain pending.

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

Run the two-dimensional GLM divergence-cleaning validation with:

```bash
python scripts/run_divergence_cleaning.py
```

This compares damped and undamped evolution of a controlled sinusoidal
\(\nabla\cdot\mathbf B\) perturbation and writes
`figures/divergence_cleaning.png` and
`benchmarks/convergence/divergence_cleaning.json`.

Run the Orszag–Tang vortex and its 64/128-cell resolution comparison with:

```bash
python scripts/run_orszag_tang.py
```

The generated figure contains density, gas and magnetic pressure, current
density, velocity, and magnetic-divergence fields. Parameters and diagnostics
are saved beside it in `figures/orszag_tang.json`.

Run the fast magnetic-rotor benchmark with:

```bash
python scripts/run_magnetic_rotor.py
```

The default command compares 64² and 128² runs at \(t=0.15\), using outflow
boundaries, and saves `figures/magnetic_rotor.png` and its JSON diagnostics.

Run the analytical magnetic-diffusion convergence test with:

```bash
python scripts/run_magnetic_diffusion.py
```

The benchmark uses 16² through 128² grids, verifies the explicit parabolic
timestep restriction and total-energy conservation, and saves its measured
orders in `benchmarks/convergence/magnetic_diffusion_convergence.json`.

Validate the unperturbed Harris current-sheet equilibrium with:

```bash
python scripts/run_harris_equilibrium.py
```

The benchmark uses periodic boundaries along the sheet and outflow boundaries
across it, compares 64² and 128² grids through `t=0.1`, and records numerical
drift before any reconnection perturbation is introduced. It writes
`figures/harris_equilibrium.png` and `figures/harris_equilibrium.json`.

Run the seeded resistive experiment and its ideal control with:

```bash
python scripts/run_reconnection.py
```

The default calculation compares 64² and 128² grids through `t=1` at
`eta=0.001`, measures reconnected flux and X-point electric field, and saves
`figures/harris_reconnection.png` plus its JSON record. The ideal control is
mandatory because HLL contributes significant numerical resistivity at these
resolutions.

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
src/mhd_solver/mhd2d/      Cartesian GLM-MHD equations, solver, boundaries
src/mhd_solver/resistive/  resistive operator and analytical diffusion state
src/mhd_solver/reconnection/  Harris-sheet and future reconnection states
scripts/                   reproducible benchmark entry points
tests/                     automated unit and integration tests
docs/                      equations, method, and validation record
figures/                   generated benchmark figures
benchmarks/                future convergence and performance studies
```

## Reproducibility and limitations

Settings use immutable configuration objects and are exposed by the benchmark
scripts. The uniform-grid solvers have no AMR, viscosity, thermal conduction,
or non-ideal equation of state. Outflow and periodic boundaries may be selected
independently in x and y. HLL is robust but smears contacts, Alfvénic waves, and
thin current sheets. MUSCL is formally second order only in smooth regions.
Multidimensional divergence is controlled with GLM cleaning rather than
constrained transport. The Harris state is not discretely well-balanced, so its
small equilibrium drift decreases with refinement but is not identically zero.

See [the equations](docs/equations.md), [the numerical method](docs/numerical_method.md),
and [the validation record](docs/validation.md) for scientific details.

## Roadmap

1. Hydrodynamics baseline (validated)
2. 1D ideal MHD: HLL, Brio-Wu, and Alfvén convergence (validated)
3. 2D ideal MHD: GLM, Orszag–Tang, and magnetic rotor (validated)
4. Uniform resistive MHD: magnetic diffusion validated
5. Harris current sheet: equilibrium and first seeded experiment validated;
   numerical-resistivity study and Sweet–Parker scaling pending
