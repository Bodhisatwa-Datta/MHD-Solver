# Numerical method: Phase 1

## Finite-volume update

The domain is divided into uniform cells of width \(\Delta x\), storing
cell-averaged conserved states. The semi-discrete update is

\[
\frac{d\mathbf{U}_i}{dt}=-\frac{
\widehat{\mathbf{F}}_{i+1/2}-\widehat{\mathbf{F}}_{i-1/2}}{\Delta x}.
\]

Two ghost cells are used on each side. They are filled by constant extrapolation
for transmissive/outflow boundaries or wrapped from the opposite edge for
periodic boundaries. The latter makes the finite-volume flux telescope exactly,
up to floating-point roundoff.

## HLL flux

For left and right primitive states, the signal estimates are

\[
S_L=\min(v_L-c_L,v_R-c_R),\qquad
S_R=\max(v_L+c_L,v_R+c_R).
\]

The HLL flux uses the left physical flux for \(S_L\geq0\), the right physical
flux for \(S_R\leq0\), and otherwise

\[
\widehat{\mathbf{F}}_{HLL}=\frac{S_R\mathbf{F}_L-S_L\mathbf{F}_R+
S_LS_R(\mathbf{U}_R-\mathbf{U}_L)}{S_R-S_L}.
\]

HLL is simple and robust, but its two-wave model diffuses the contact wave.

## Reconstruction and time integration

The first-order method passes piecewise-constant cell states to HLL and advances
with forward Euler. The second-order method reconstructs density, velocity, and
pressure using component-wise minmod slopes,

\[
\operatorname{minmod}(a,b)=
\begin{cases}\operatorname{sign}(a)\min(|a|,|b|),&ab>0,\\0,&ab\leq0.
\end{cases}
\]

If an interface reconstruction produces non-positive density or pressure, that
one side reverts to its cell average by suppressing its slope. The physical state
is not clipped. Time integration uses the two-stage strong-stability-preserving
Runge-Kutta method (SSP-RK2).

## CFL condition

Each timestep obeys

\[
\Delta t=C_{CFL}\frac{\Delta x}{\max_i(|v_i|+c_i)},
\]

and the last step is shortened to reach the requested final time exactly. The
default CFL number is 0.8.
