# Governing equations

The one-dimensional compressible Euler equations are solved in conservation form,

\[
\frac{\partial \mathbf{U}}{\partial t} +
\frac{\partial \mathbf{F}(\mathbf{U})}{\partial x}=0,
\]

with conserved state and flux

\[
\mathbf{U}=
\begin{pmatrix}\rho\\ \rho v\\ E\end{pmatrix},\qquad
\mathbf{F}=
\begin{pmatrix}\rho v\\ \rho v^2+p\\ v(E+p)\end{pmatrix}.
\]

Here \(\rho\) is mass density, \(v\) is velocity, and \(E\) is total energy
density. The ideal-gas closure is

\[
p=(\gamma-1)\left(E-\frac{1}{2}\rho v^2\right),
\]

and the adiabatic sound speed is \(c=\sqrt{\gamma p/\rho}\). The code uses
\(\gamma=1.4\) by default. Valid states require finite values, \(\rho>0\), and
\(p>0\); violations raise an exception with the failing quantity.

There are no source terms, viscosity, thermal conduction, magnetic fields, or
multi-dimensional terms in Phase 1.

## Exact smooth advection solution

For constant velocity \(v_0\) and pressure \(p_0\), an arbitrary density profile
is an entropy/contact mode of the Euler system and is advected without changing
shape. The convergence benchmark uses

\[
\rho(x,t)=1+0.2\sin\left[2\pi(x-v_0t)\right],\qquad
v_0=1,\qquad p_0=1
\]

on a unit periodic domain. This is an exact nonlinear Euler solution, so its
error is not contaminated by a linearisation or separately computed reference.

# One-dimensional ideal MHD

In units with magnetic permeability equal to one, the evolved state is

\[
\mathbf U=(\rho,\rho v_x,\rho v_y,\rho v_z,E,B_y,B_z)^T.
\]

The longitudinal field \(B_x\) is spatially constant in 1D and is supplied as a
parameter. Total energy and total pressure are

\[
E=\frac{p}{\gamma-1}+\frac{1}{2}\rho|\mathbf v|^2+
\frac{1}{2}|\mathbf B|^2,\qquad
p_t=p+\frac{1}{2}|\mathbf B|^2.
\]

The x-directed flux is

\[
\mathbf F=\begin{pmatrix}
\rho v_x\\
\rho v_x^2+p_t-B_x^2\\
\rho v_xv_y-B_xB_y\\
\rho v_xv_z-B_xB_z\\
(E+p_t)v_x-B_x(\mathbf v\cdot\mathbf B)\\
v_xB_y-v_yB_x\\
v_xB_z-v_zB_x
\end{pmatrix}.
\]

Because \(B_x\) is constant, \(\nabla\cdot\mathbf B=\partial_xB_x=0\)
identically. This one-dimensional property does not replace divergence control
in a future multidimensional solver.

## Exact circularly polarized Alfvén wave

The smooth MHD convergence problem uses constant \(\rho_0\), \(p_0\), and
\(B_x\), with

\[
B_y=A\cos(k[x-v_At]),\qquad B_z=A\sin(k[x-v_At]),
\]

\[
v_y=-\frac{B_y}{\sqrt{\rho_0}},\qquad
v_z=-\frac{B_z}{\sqrt{\rho_0}},\qquad
v_A=\frac{B_x}{\sqrt{\rho_0}}.
\]

Because \(|\mathbf B_\perp|=A\) is constant, magnetic pressure is uniform. This
circular polarization is an exact finite-amplitude ideal-MHD solution rather
than a linearized wave.

# Two-dimensional GLM-MHD

The cell-centred 2D formulation augments ideal MHD with a generalized Lagrange
multiplier \(\psi\). Its conserved state is

\[
\mathbf U=(\rho,\rho v_x,\rho v_y,\rho v_z,E,B_x,B_y,B_z,\psi)^T,
\]

where the augmented total energy is

\[
E=\frac{p}{\gamma-1}+\frac{1}{2}\rho|\mathbf v|^2+
\frac{1}{2}|\mathbf B|^2+\frac{\psi^2}{2c_h^2}.
\]

The induction and cleaning subsystem is

\[
\frac{\partial\mathbf B}{\partial t}
+\nabla\cdot(\mathbf v\mathbf B-\mathbf B\mathbf v+\psi\mathbf I)=0,
\]

\[
\frac{\partial\psi}{\partial t}+c_h^2\nabla\cdot\mathbf B=-\kappa\psi.
\]

Thus divergence errors propagate at speed \(c_h\) and decay through the damping
rate \(\kappa\). The implementation is based on the hyperbolic/parabolic GLM
idea of Dedner et al., *Journal of Computational Physics* 175 (2002), 645–673,
[doi:10.1006/jcph.2001.6961](https://doi.org/10.1006/jcph.2001.6961).

# Uniform resistive MHD

For constant scalar resistivity \(\eta\), the induction equation becomes

\[
\frac{\partial\mathbf B}{\partial t}
=\nabla\times(\mathbf v\times\mathbf B)+\eta\nabla^2\mathbf B,
\]

where the Laplacian form assumes \(\nabla\cdot\mathbf B=0\). The implementation
also includes the resistive contribution to the conservative total-energy flux,

\[
\mathbf F_{E,\eta}=\eta\,\mathbf J\times\mathbf B,
\qquad \mathbf J=\nabla\times\mathbf B.
\]

This transports resistive electromagnetic energy and converts decaying magnetic
energy into internal energy while conserving periodic-domain total energy.

For a stationary, weak, divergence-free sinusoidal field,

\[
B_y(x,0)=A\sin(kx),
\]

the magnetic diffusion equation has the analytical solution

\[
B_y(x,t)=A\exp(-\eta k^2t)\sin(kx).
\]
