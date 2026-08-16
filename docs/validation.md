# Validation record

## Sod shock tube

The initial discontinuity is at \(x=0.5\) on \([0,1]\), with

| State | Density | Velocity | Pressure |
|---|---:|---:|---:|
| Left | 1.0 | 0.0 | 1.0 |
| Right | 0.125 | 0.0 | 0.1 |

The gas has \(\gamma=1.4\), boundaries are transmissive, and the target time is
\(t=0.2\). `scripts/run_sod.py` runs both first-order HLL and MUSCL-minmod HLL with
SSP-RK2, compares cell-centred values with the exact self-similar Riemann
solution, and records mean absolute (discrete \(L_1\)) errors. It also writes the
configuration, step counts, and changes in the domain integrals of mass,
momentum, and total energy alongside the figure as JSON.

Numerical results are added below only after running the committed code.

<!-- BENCHMARK_RESULTS -->

## Scope of validation

The Sod test exercises a rarefaction, contact discontinuity, and shock. Agreement
with it is necessary but not sufficient to validate a general Euler solver.
Future hydrodynamic work should include smooth-wave convergence and stronger
shock tests before extending the numerical framework to MHD.

