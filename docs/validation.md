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

## Smooth periodic density wave

A sinusoidal density profile with amplitude 0.2, uniform velocity 1, and uniform
pressure 1 is advected across the unit periodic domain. At \(t=1\) it has
completed one crossing and must equal its initial state exactly. Calculations
used \(\gamma=1.4\), CFL 0.6, HLL fluxes, and resolutions 50, 100, 200, and 400.

| Cells | First-order \(L_1\) | Order | MUSCL + RK2 \(L_1\) | Order |
|---:|---:|---:|---:|---:|
| 50  | 3.91613e-2 | —    | 7.46898e-3 | —    |
| 100 | 2.14014e-2 | 0.872 | 2.12666e-3 | 1.812 |
| 200 | 1.11960e-2 | 0.935 | 5.90722e-4 | 1.848 |
| 400 | 5.72710e-3 | 0.967 | 1.60072e-4 | 1.884 |

The first-order method approaches its expected order of one. MUSCL with minmod
and SSP-RK2 approaches second order; the remaining reduction is consistent with
the limiter becoming first order at smooth extrema. Across every run, absolute
changes in integrated mass, momentum, and energy were at most
\(4.44\times10^{-16}\), demonstrating conservative periodic evolution to
roundoff. The figure and full machine-readable record are
`figures/density_wave_convergence.png` and
`benchmarks/convergence/density_wave_convergence.json`.

## Scope of validation

Together, the Sod problem and smooth-wave study exercise shocks, a rarefaction,
contact transport, formal convergence, and periodic conservation. A stronger
shock test would further broaden the hydrodynamic validation envelope, but the
present baseline now covers the principal requirements needed before the 1D MHD
extension.

## Brio-Wu ideal-MHD shock tube

The standard left and right primitive states are

| State | \(\rho\) | \(p\) | \(B_y\) | \(v_x,v_y,v_z\) | \(B_z\) |
|---|---:|---:|---:|---:|---:|
| Left  | 1.0   | 1.0 | 1.0  | 0 | 0 |
| Right | 0.125 | 0.1 | -1.0 | 0 | 0 |

with \(B_x=0.75\), \(\gamma=2\), a discontinuity at \(x=0\), and transmissive
boundaries on \([-0.5,0.5]\). The benchmark originates with Brio and Wu,
*Journal of Computational Physics* 75 (1988), 400–422,
[doi:10.1016/0021-9991(88)90120-9](https://doi.org/10.1016/0021-9991(88)90120-9).

The committed run reached \(t=0.1\) using HLL, MUSCL-minmod, SSP-RK2, CFL 0.4,
and 800 cells. A separate 1600-cell run was overlaid and restricted to the
800-cell mesh for resolution sensitivity:

| Quantity | 800 vs restricted-1600 \(L_1\) difference |
|---|---:|
| Density | 1.78909e-3 |
| Longitudinal velocity | 3.17995e-3 |
| Transverse velocity | 4.53666e-3 |
| Gas pressure | 1.69215e-3 |
| Transverse field \(B_y\) | 2.30918e-3 |

The state remained finite, with minimum density 0.11707 and minimum gas pressure
0.08773. The profiles show the multi-wave MHD structure, and the resolutions
agree closely away from discontinuities. This is resolution sensitivity, not an
error measured against an exact solution.

Mass and total energy changed by at most roundoff. Integrated x-momentum changed
by 0.09, matching the boundary total-pressure imbalance over \(t=0.1\), while
y-momentum changed by -0.15, matching the boundary magnetic-stress flux. Since
\(B_x\) is constant by construction, the one-dimensional divergence error is
zero. Complete parameters and extrema are stored in `figures/brio_wu.json`.
