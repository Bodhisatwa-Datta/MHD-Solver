# Numerical method

## Finite-volume update

The domain is divided into uniform cells of width $\Delta x$, storing
cell-averaged conserved states. The semi-discrete update is

$$
\frac{d\mathbf{U}_i}{dt}=-\frac{
\widehat{\mathbf{F}}_{i+1/2}-\widehat{\mathbf{F}}_{i-1/2}}{\Delta x}.
$$

Two ghost cells are used on each side. They are filled by constant extrapolation
for transmissive/outflow boundaries or wrapped from the opposite edge for
periodic boundaries. The latter makes the finite-volume flux telescope exactly,
up to floating-point roundoff.

## HLL flux

For left and right primitive states, the signal estimates are

$$
S_L=\min(v_L-c_L,v_R-c_R),\qquad
S_R=\max(v_L+c_L,v_R+c_R).
$$

The HLL flux uses the left physical flux for $S_L\geq0$, the right physical
flux for $S_R\leq0$, and otherwise

$$
\widehat{\mathbf{F}}_{HLL}=\frac{S_R\mathbf{F}_L-S_L\mathbf{F}_R+
S_LS_R(\mathbf{U}_R-\mathbf{U}_L)}{S_R-S_L}.
$$

HLL is simple and robust, but its two-wave model diffuses the contact wave.

## Reconstruction and time integration

The first-order method passes piecewise-constant cell states to HLL and advances
with forward Euler. The second-order method reconstructs density, velocity, and
pressure using component-wise minmod slopes,

$$
\operatorname{minmod}(a,b)=
\begin{cases}\operatorname{sign}(a)\min(|a|,|b|),&ab>0,\\0,&ab\leq0.
\end{cases}
$$

If an interface reconstruction produces non-positive density or pressure, that
one side reverts to its cell average by suppressing its slope. The physical state
is not clipped. Time integration uses the two-stage strong-stability-preserving
Runge-Kutta method (SSP-RK2).

## CFL condition

Each timestep obeys

$$
\Delta t=C_{CFL}\frac{\Delta x}{\max_i(|v_i|+c_i)},
$$

and the last step is shortened to reach the requested final time exactly. The
default CFL number is 0.8.

# One-dimensional ideal MHD

The MHD solver reuses the finite-volume mesh, boundary fills, minmod
reconstruction, and SSP-RK2 method. Primitive reconstruction acts on
$(\rho,v_x,v_y,v_z,p,B_y,B_z)$, with density and gas-pressure checks before
interface states enter the Riemann solver.

For propagation along x, define

$$
a^2=\frac{\gamma p}{\rho},\quad
v_A^2=\frac{|\mathbf B|^2}{\rho},\quad
v_{Ax}^2=\frac{B_x^2}{\rho}.
$$

The fast magnetosonic speed is

$$
c_f^2=\frac{1}{2}\left[a^2+v_A^2+
\sqrt{(a^2+v_A^2)^2-4a^2v_{Ax}^2}\right].
$$

MHD HLL uses $S_L=\min(v_{x,L}-c_{f,L},v_{x,R}-c_{f,R})$ and
$S_R=\max(v_{x,L}+c_{f,L},v_{x,R}+c_{f,R})$. The CFL speed is
$\max_i(|v_{x,i}|+c_{f,i})$. HLL is deliberately used before HLLD: it is
robust and auditable, but diffuses contact and Alfvénic waves because it retains
only the two outer signal speeds.

# Two-dimensional GLM-MHD

The 2D solver stores cell averages on a uniform Cartesian mesh and applies the
unsplit method-of-lines operator

$$
\mathcal L(\mathbf U)=-\frac{\widehat{\mathbf F}_{i+1/2,j}-
\widehat{\mathbf F}_{i-1/2,j}}{\Delta x}
-\frac{\widehat{\mathbf G}_{i,j+1/2}-
\widehat{\mathbf G}_{i,j-1/2}}{\Delta y}.
$$

Primitive variables are independently reconstructed in x and y with minmod
slopes. Directional HLL signal bounds include both the fast magnetosonic speed
and the GLM cleaning speed. The multidimensional timestep is

$$
\Delta t=C_{CFL}\left(\frac{a_x}{\Delta x}+
\frac{a_y}{\Delta y}\right)^{-1},
$$

where $a_x$ and $a_y$ are the maximum directional transport speeds including
$c_h$. SSP-RK2 advances the hyperbolic system. The exactly integrated damping
source $\psi\leftarrow\psi e^{-\kappa\Delta t}$ is Strang-split into half
steps. Removed $\psi$-energy is also removed from augmented total energy, so
damping does not spuriously heat the gas.

The diagnostic $\|\nabla\cdot\mathbf B\|_2$ uses centred differences. GLM is
simple and compatible with cell-centred fields, but it controls rather than
eliminates divergence. Constrained transport is intentionally deferred.

Periodic wrapping and zero-gradient outflow fills can be chosen independently
in x and y. The Harris-sheet validation uses periodic x boundaries along the
sheet and outflow y boundaries across it.

# Explicit resistive update

Uniform resistivity adds centred second-order Laplacians to all three magnetic
components. Cell-centred currents are computed with centred differences, and
the conservative energy contribution is
$-\nabla\cdot[\eta(\mathbf J\times\mathbf B)]$. These terms are included in
both SSP-RK2 stages rather than applied as a first-order operator split.

In addition to the hyperbolic CFL condition, each timestep satisfies

$$
\Delta t_\eta \leq C_\eta
\left[2\eta\left(\Delta x^{-2}+\Delta y^{-2}\right)\right]^{-1},
$$

with $C_\eta=0.8$ by default. This restriction becomes expensive at high
resolution because $\Delta t_\eta\propto\Delta x^2$; super-time-stepping or an
implicit method may be justified later, but is deliberately not introduced
before the explicit implementation is validated.

## Effective numerical resistivity

Numerical magnetic diffusion is measured by evolving the weak mode
$B_y=A_0\sin(kx)$ with explicit $\eta=0$. Its final amplitude is obtained
from the discrete Fourier projection

$$
A(t)=2\langle B_y(x,t)\sin(kx)\rangle.
$$

An effective resistivity is then defined from the analytical diffusion law,

$$
\eta_{num}=-\frac{\ln[A(t)/A_0]}{k^2t}.
$$

The amplitude $A_0=10^{-3}$ makes Lorentz-force feedback quadratic and
negligible: measured maximum velocities remain below $3\times10^{-8}$. This
diagnostic characterizes smooth-mode damping for a chosen grid, CFL number, and
scheme. It is not a universal material coefficient and can underestimate
dissipation where the minmod limiter reduces the reconstruction to first order.
As a consistency check, the same mode is evolved with a known explicit
resistivity. Subtracting that value from the total inferred decay should recover
the independently measured zero-resistivity $\eta_{num}$ if the two diffusion
mechanisms are additive in this weak, smooth regime.

## Exact-time diagnostic snapshots

The 2D configuration may request a strictly increasing tuple of output times.
The timestep is shortened when necessary to land exactly on each requested
time, and a copy of the primitive state is stored. No temporal interpolation is
used. Empty output schedules retain the original low-memory behavior.

The reconnection sweep samples eleven states from $t=0$ through $t=1$.
Reconnected flux, X-point electric field, kinetic and magnetic energies,
primitive extrema, current, and divergence are evaluated from those synchronized
states. A descriptive late-time flux slope is obtained by least-squares fitting
$\Psi_{rec}(t)$ over $0.5\leq t\leq1$. It is a finite-window diagnostic, not
automatically a steady reconnection rate.
