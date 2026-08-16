# Numerical method

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

# One-dimensional ideal MHD

The MHD solver reuses the finite-volume mesh, boundary fills, minmod
reconstruction, and SSP-RK2 method. Primitive reconstruction acts on
\((\rho,v_x,v_y,v_z,p,B_y,B_z)\), with density and gas-pressure checks before
interface states enter the Riemann solver.

For propagation along x, define

\[
a^2=\frac{\gamma p}{\rho},\quad
v_A^2=\frac{|\mathbf B|^2}{\rho},\quad
v_{Ax}^2=\frac{B_x^2}{\rho}.
\]

The fast magnetosonic speed is

\[
c_f^2=\frac{1}{2}\left[a^2+v_A^2+
\sqrt{(a^2+v_A^2)^2-4a^2v_{Ax}^2}\right].
\]

MHD HLL uses \(S_L=\min(v_{x,L}-c_{f,L},v_{x,R}-c_{f,R})\) and
\(S_R=\max(v_{x,L}+c_{f,L},v_{x,R}+c_{f,R})\). The CFL speed is
\(\max_i(|v_{x,i}|+c_{f,i})\). HLL is deliberately used before HLLD: it is
robust and auditable, but diffuses contact and Alfvénic waves because it retains
only the two outer signal speeds.
