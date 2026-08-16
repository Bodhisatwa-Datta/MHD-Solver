# Governing equations: Phase 1

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
