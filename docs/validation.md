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

## Circularly polarized Alfvén wave

The convergence test uses \(\rho=1\), \(p=0.1\), \(B_x=1\), transverse
amplitude 0.1, wavelength 1, and \(\gamma=5/3\) on the periodic unit domain.
The Alfvén speed is one, so the solution returns to its initial position at
\(t=1\). Runs used HLL fluxes and CFL 0.4.

The error norm is the mean magnitude of the transverse magnetic-field error,

\[
L_1=\frac{1}{N}\sum_i\sqrt{(B_{y,i}-B_{y,i}^{exact})^2+
(B_{z,i}-B_{z,i}^{exact})^2}.
\]

| Cells | First-order \(L_1\) | Order | MUSCL + RK2 \(L_1\) | Order |
|---:|---:|---:|---:|---:|
| 50  | 2.13331e-2 | —     | 5.13395e-3 | —     |
| 100 | 1.13071e-2 | 0.916 | 1.54999e-3 | 1.728 |
| 200 | 5.82433e-3 | 0.957 | 4.31538e-4 | 1.845 |
| 400 | 2.95661e-3 | 0.978 | 1.17674e-4 | 1.875 |

The first-order method approaches order one and MUSCL with SSP-RK2 approaches
order two. The limiter reduces the observed rate near smooth extrema, as in the
hydrodynamic density-wave test. At 400 cells, the second-order transverse
velocity error was \(1.16983\times10^{-4}\). Across all resolutions, changes in
the seven evolved integrals were no larger than \(2.22\times10^{-16}\), and
\(\partial_x B_x=0\) exactly by construction.

The measured data and figure are stored in
`benchmarks/convergence/alfven_wave_convergence.json` and
`figures/alfven_wave_convergence.png`.

## Two-dimensional GLM divergence cleaning

A small controlled perturbation
\(B_x=10^{-3}\sin(2\pi x)\), with \(B_y=B_z=0\), was initialized on a periodic
64×64 grid with \(\rho=p=1\). Both runs used MUSCL, SSP-RK2, HLL, CFL 0.35, and
cleaning speed \(c_h=2\); one used no parabolic damping and the other used
\(\kappa=2\).

At \(t=1\), which compares the same cleaning-wave phase, the measured results
were:

| Diagnostic | Value |
|---|---:|
| Initial \(\|\nabla\cdot\mathbf B\|_2\) | 4.43575e-3 |
| Final undamped norm | 4.33942e-3 |
| Final damped norm | 1.59183e-3 |
| Damped reduction factor | 2.7866 |
| Damped total-energy change | -2.13037e-7 |

The norm oscillates because the hyperbolic subsystem exchanges divergence error
with \(\psi\); its envelope decays only in the damped run. The small energy loss
is intentional removal of cleaning/divergence energy and is comparable to the
initial perturbation magnetic energy. Minimum density and pressure remained
0.999999995 and 0.999999994. Uniform 2D states are preserved to roundoff for
both numerical orders.

This test validates the cleaning subsystem in isolation. It does not yet
validate nonlinear multidimensional MHD dynamics. Orszag–Tang is the next
required benchmark. Data and the diagnostic figure are stored in
`benchmarks/convergence/divergence_cleaning.json` and
`figures/divergence_cleaning.png`.

## Orszag–Tang vortex

The periodic unit-square initial condition uses \(\gamma=5/3\),

\[
\rho=\frac{25}{36\pi},\qquad p=\frac{5}{12\pi},
\]

\[
\mathbf v=(-\sin 2\pi y,\ \sin 2\pi x,\ 0),
\]

\[
\mathbf B=\frac{1}{\sqrt{4\pi}}
(-\sin 2\pi y,\ \sin 4\pi x,\ 0).
\]

This is the standard compressible setup used by established MHD codes such as
the [Athena test suite](https://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html).
The discrete initial magnetic divergence is zero to roundoff. Runs used HLL,
MUSCL-minmod, SSP-RK2, CFL 0.3, \(c_h=2\), \(\kappa=2\), and periodic boundaries.

At \(t=0.5\), the 128×128 calculation produced interacting shocks and thin
current sheets while remaining physical:

| Diagnostic | Measured value |
|---|---:|
| Density range | 0.10436–0.41824 |
| Gas-pressure range | 0.03999–0.44580 |
| Maximum velocity magnitude | 1.53104 |
| Maximum magnetic pressure | 0.23106 |
| Maximum \(|J_z|\) | 25.9751 |
| Final kinetic energy | 0.0430150 |
| Final magnetic energy | 0.0527245 |
| Final \(\|\nabla\cdot\mathbf B\|_2\) | 0.0709798 |
| \(\Delta x\|\nabla\cdot B\|_2/B_{rms}\) | 1.70767e-3 |

Mass and both momentum components changed by no more than \(8.33\times10^{-17}\).
Augmented total energy decreased by \(5.05\times10^{-7}\) through the documented
GLM damping sink.

The 64×64 solution was compared with the volume-restricted 128×128 result:

| Quantity | Coarse/fine \(L_1\) difference |
|---|---:|
| Density | 1.33929e-2 |
| Pressure | 1.71532e-2 |
| \(B_x\) | 2.83308e-2 |
| \(B_y\) | 3.28702e-2 |

The shock topology and large-scale structures are consistent between the two
resolutions, but these differences show that discontinuities and current sheets
are not grid-converged. The result is therefore a qualitative nonlinear
validation and resolution-sensitivity study, not a converged reference
solution. Full diagnostics and the six-panel figure are stored in
`figures/orszag_tang.json` and `figures/orszag_tang.png`.

## Fast magnetic rotor

The Balsara–Spicer rotor tests the propagation of strong torsional Alfvén waves
from a rapidly spinning dense core. The setup follows the standard benchmark
introduced by Balsara and Spicer, *Journal of Computational Physics* 149 (1999),
270–292, [doi:10.1006/jcph.1998.6153](https://doi.org/10.1006/jcph.1998.6153).

On \([0,1]^2\), the rotor is centred at \((0.5,0.5)\), with inner and outer
radii \(r_0=0.1\) and \(r_1=0.115\). Inside the rotor, \(\rho=10\) and angular
velocity is 20; density and velocity taper linearly to the ambient stationary
state with \(\rho=1\). Pressure is one everywhere, \(B_x=5/\sqrt{4\pi}\),
\(B_y=B_z=0\), and \(\gamma=1.4\). The field is initially divergence-free.

The committed calculation used HLL, MUSCL-minmod, SSP-RK2, outflow boundaries,
CFL 0.25, \(c_h=4\), and \(\kappa=4\). At \(t=0.15\), the 128² run measured:

| Diagnostic | Measured value |
|---|---:|
| Density range | 0.90416–5.47109 |
| Pressure range | 0.03706–1.79325 |
| Maximum velocity magnitude | 0.89534 |
| Maximum magnetic pressure | 2.11482 |
| Maximum \(|J_z|\) | 51.8341 |
| Final kinetic energy | 0.111291 |
| Final magnetic energy | 1.104737 |
| Final \(\|\nabla\cdot\mathbf B\|_2\) | 0.0818575 |
| \(\Delta x\|\nabla\cdot B\|_2/B_{rms}\) | 4.30233e-4 |

The initially horizontal field is wound around the rotor and launches outward
torsional structures. Density forms an elongated core and current concentrates
along the twisted-field fronts, qualitatively matching the standard benchmark.
Density and pressure remain positive.

Because the benchmark uses outflow boundaries and GLM waves reach those
boundaries, changes in global integrals include physical/numerical boundary
fluxes: mass changed by \(-4.35\times10^{-6}\) and total energy by
\(3.56\times10^{-5}\). They are not interpreted as closed-domain conservation
errors.

The 64² result was compared with the restricted 128² solution:

| Quantity | Coarse/fine \(L_1\) difference |
|---|---:|
| Density | 1.13636e-1 |
| Pressure | 4.63368e-2 |
| \(B_x\) | 3.28108e-2 |
| \(B_y\) | 2.64766e-2 |

The large-scale morphology agrees, but the density difference confirms that the
rotor edge and torsional fronts are not grid-converged. This is a qualitative
benchmark and resolution-sensitivity result. Full diagnostics and the generated
figure are `figures/magnetic_rotor.json` and `figures/magnetic_rotor.png`.

## Uniform magnetic diffusion

The resistive implementation was tested using the analytical decay

\[
B_y(x,t)=10^{-3}\exp[-\eta(2\pi)^2t]\sin(2\pi x)
\]

on a periodic unit square with \(\eta=0.1\), \(\rho=p=1\), \(\gamma=5/3\), and
zero velocity. The weak amplitude isolates the diffusion operator from
Lorentz-force feedback while still evolving the complete MHD state. Runs reached
\(t=0.05\) using MUSCL, SSP-RK2, HLL, and the explicit diffusion timestep limit.

| Cells per direction | \(B_y\) mean absolute error | Order |
|---:|---:|---:|
| 16  | 1.10537e-5 | — |
| 32  | 2.99154e-6 | 1.886 |
| 64  | 7.57915e-7 | 1.981 |
| 128 | 1.89362e-7 | 2.001 |

The measured rate approaches the expected second order. At 128², numerical
magnetic energy was \(1.6838526\times10^{-7}\), compared with the analytical
\(1.6845636\times10^{-7}\). Total-energy changes across all resolutions were no
larger than \(1.33\times10^{-15}\). The finest pressure remained between
0.999999939 and 1.000000170 as magnetic energy was redistributed into heat.

This validates uniform resistive diffusion and its energy coupling in a smooth,
weak-field regime. It does not yet validate spatially varying resistivity,
current-sheet evolution, or reconnection. The measured data and figure are
`benchmarks/convergence/magnetic_diffusion_convergence.json` and
`figures/magnetic_diffusion_convergence.png`.

## Unperturbed Harris current sheet

Before introducing a reconnection seed, the ideal equilibrium was evolved on
\([0,1]\times[-0.5,0.5]\) with periodic x boundaries and outflow y boundaries.
The profiles use \(B_0=1\), half-width \(L=0.05\), background pressure 0.2,
background density 1, and therefore constant temperature \(p/\rho=0.2\). Gas
and magnetic pressure sum exactly to 0.7 in the initial cell-centred state.
Runs used \(\gamma=5/3\), HLL, MUSCL-minmod, SSP-RK2, CFL 0.3, \(c_h=2\),
\(\kappa=2\), and \(\eta=0\), reaching \(t=0.1\).

| Diagnostic | 64² | 128² |
|---|---:|---:|
| Maximum velocity magnitude | 1.61973e-2 | 7.31778e-3 |
| Density mean absolute drift | 4.16113e-2 | 1.22978e-2 |
| Pressure mean absolute drift | 6.48233e-3 | 1.68035e-3 |
| \(B_x\) mean absolute drift | 7.53765e-3 | 1.96195e-3 |
| Maximum total-pressure error | 4.28531e-2 | 1.38532e-2 |
| Maximum \(|J_z|\) | 17.8200 | 20.1640 |
| Final \(\|\nabla\cdot\mathbf B\|_2\) | 0 | 0 |

Refinement reduces every drift measure substantially: the density, pressure,
and field drifts fall by factors 3.38, 3.86, and 3.84, respectively, while the
maximum spurious velocity falls by a factor 2.21. The initial analytic peak
current magnitude is \(B_0/L=20\), approached by the 128² discrete result.
All nine evolved domain integrals are conserved to roundoff at 128², and
divergence remains exactly zero because the solution preserves its x-invariant
structure with \(B_y=0\).

This demonstrates a stable, resolution-improving numerical equilibrium, not an
exactly well-balanced discretization. HLL diffusion broadens the sheet and the
discrete pressure/magnetic-stress balance launches small y-directed motions.
Those measured errors set the baseline that a future seeded resistive
reconnection run must exceed. The figure and full record are
`figures/harris_equilibrium.png` and `figures/harris_equilibrium.json`.

## Seeded Harris-sheet reconnection

A divergence-free perturbation with vector-potential amplitude 0.01 and
vertical width 0.1 was applied to the validated sheet. The X-point is at
\((0.5,0)\), while the magnetic island closes through the periodic x boundary.
Runs used the equilibrium configuration above, \(\eta=10^{-3}\), and reached
\(t=1\). The nominal Lundquist numbers are 50 based on sheet half-width 0.05
and 500 based on sheet half-length 0.5. Each resistive calculation was paired
with an otherwise identical ideal control.

| Diagnostic | 64² resistive | 128² resistive |
|---|---:|---:|
| Initial reconnected flux | 1.98544e-2 | 1.99635e-2 |
| Final reconnected flux | 2.25983e-2 | 2.27008e-2 |
| Flux increase | 2.74391e-3 | 2.73726e-3 |
| Final \(E_z(X)\) | -1.29767e-2 | -1.70845e-2 |
| Final kinetic energy | 1.74299e-4 | 3.65667e-4 |
| Final magnetic energy | 4.07101e-1 | 4.24546e-1 |
| Maximum speed | 4.21745e-2 | 7.77314e-2 |
| Maximum \(|J_z|\) | 12.0737 | 16.7286 |
| Dimensionless divergence | 1.61114e-4 | 1.18009e-4 |

The flux increase agrees to 0.24% between the two grids. The 128² field-line
plot shows the seeded O-point, a central X-point, inflow toward the sheet, and
oppositely directed outflow jets. Magnetic energy falls from 0.450435 to
0.424546 while kinetic energy grows from zero. The remaining energy is mainly
redistributed into internal energy, with additional exchange through the open
y boundaries. At 128², mass and total energy change by 0.00526 and 0.00785;
these are boundary-flux balances, not closed-domain conservation errors.

Detailed structures are not grid-converged: the restricted 128² and 64²
resistive states differ in density, pressure, \(B_x\), and \(B_y\) by mean
absolute values 0.0663, 0.0175, 0.0224, and 0.00514. The local X-point electric
field also remains resolution-sensitive.

The ideal controls are decisive. Their flux increases are 0.00516 and 0.00584
at 64² and 128², larger than the explicit-resistive increases. Thus the run
demonstrates topology change and exercises the reconnection diagnostics, but
does not isolate a reconnection rate controlled solely by physical resistivity.
No Sweet–Parker scaling is claimed. Establishing such scaling requires a
numerical-resistivity study, additional resolutions and resistivities, and
likely a less diffusive Riemann solver or constrained transport. Full results
are stored in `figures/harris_reconnection.json` and
`figures/harris_reconnection.png`.
