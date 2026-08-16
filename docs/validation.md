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

### Measured result (2026-08-16)

The committed benchmark was run at 400 cells with CFL 0.8. Both schemes took 218
steps. The measured cell-centred errors were:

| Quantity | First-order \(L_1\) | MUSCL + RK2 \(L_1\) |
|---|---:|---:|
| Density | 6.75094e-3 | 2.61565e-3 |
| Velocity | 8.29800e-3 | 4.62889e-3 |
| Pressure | 4.78999e-3 | 1.78994e-3 |

The second-order configuration reduced all three measured errors. Density and
pressure remained positive, and the numerical profiles reproduce the expected
left rarefaction, contact discontinuity, and right-moving shock. Mass and total
energy changed by zero to printed precision. Domain momentum increased by 0.18;
this is the correct boundary-flux balance \((p_L-p_R)t=(1-0.1)0.2=0.18\), not a
failure of conservation in the open, pressure-loaded domain.

The generated figure and machine-readable run record are
`figures/sod_comparison.png` and `figures/sod_comparison.json`.

## Scope of validation

The Sod test exercises a rarefaction, contact discontinuity, and shock. Agreement
with it is necessary but not sufficient to validate a general Euler solver.
Future hydrodynamic work should include smooth-wave convergence and stronger
shock tests before extending the numerical framework to MHD.
