---
title: 'SARAS: A general-purpose PDE solver for fluid dynamics'

tags:
  - C++
  - PDE
  - turbulence
  - fluid dynamics

authors:
  - name: Roshan Samuel
    orcid: 0000-0002-1280-9881
    affiliation: 1
  - name: Shashwat Bhattacharya
    orcid: 0000-0001-7462-7680
    affiliation: 1
  - name: Ali Asad
    orcid: 0000-0001-9704-6686
    affiliation: 2
  - name: Soumyadeep Chatterjee
    orcid: 0000-0001-7957-1727
    affiliation: 2
  - name: Mahendra K. Verma
    orcid: 0000-0002-3380-4561
    affiliation: 2

affiliations:
 - name: Department of Mechanical Engineering, Indian Institute of Technology - Kanpur
   index: 1
 - name: Department of Physics, Indian Institute of Technology - Kanpur
   index: 2

date: 15 January 2020

bibliography: resources/paper.bib

---

# Summary

The laws that govern natural systems can often be modelled mathematically using
partial differential equations (PDEs).
Usually the resultant PDEs are not solvable analytically,
hence numerical solutions become important in such cases.
As a result, efficient numerical solutions of PDEs are important
for understanding such systems.
In this paper we briefly describe the design and validation
of a finite difference solver ``SARAS``.

``SARAS`` is a general-purpose PDE solver based on finite difference
method [@Anderson:book:CFD; @Ferziger:book:CFD], written in an object-oriented structure in C++.
In ``SARAS``, the underlying mathematical constructs like vector and scalar fields
are defined as classes.
Moreover, vector calculus operations associated with such fields, like gradient and
divergence, are defined using these classes.
This design makes the code intuitive, allowing users to quickly cast PDEs into
readable codes.
The initial conditions, boundary conditions, and source/forcing terms
are implemented using ``initial``, ``boundary`` and ``force`` classes.
These classes are readily extensible so that users can add custom
initial conditions, source terms, and so on.

``SARAS`` includes solvers for hydrodynamic flows, namely the incompressible
Navier-Stokes initial value problem (IVP), as well as for scalar convection,
like Rayleigh Benard Convection.
Presently, we use semi-implicit Crank-Nicholson [@Crank:1947] method for time-advancing
the IVP.
The solver uses Marker and Cell method (MAC) [@Harlow:PF1965] for discretizing
the velocity and pressure fields.

# Mathematics

The Navier-Stokes equations, which govern the dynamics of fluid flow, can be written as
$$
\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u} = -\nabla p + \mathbf{f} + \nu\nabla^2\mathbf{u},
$$
where $\mathbf{u}$ is the velocity field, $p$ is the pressure field, $\mathbf{f}$ is
the forcing term, and $\nu$ is the kinematic viscosity of the fluid.
The fluid is assumed to be incompressible. Hence $\nabla\cdot\mathbf{u} = 0$, and density is
constant (chosen to be unity).

If the velocity and pressure field at time $t = t_n$ are denoted as $\mathbf{u}_n$ and $p_n$
respectively, then the corresponding fields at the next time-step, $t = t_{n+1}$, namely
$\mathbf{u}_{n+1}$ and $p_{n+1}$, can be calculated as described below
[@Patankar:1972IJHMT; @Anderson:book:CFD; @Ferziger:book:CFD].
We compute an intermediate velocity field using the known values,
$\mathbf{u}_n$ and $p_n$, as
$$
\mathbf{u}^* = \mathbf{u}_{n} + \Delta t\left[\nu\nabla^2 \left( \frac{\mathbf{u}_n + \mathbf{u}^*}{2}\right) - \mathbf{u}_n.\nabla\mathbf{u}_n - \nabla p_n\right].
$$
The forcing term has been neglected here for simplicity.
Note that the diffusion term (also called the viscous term) is handled semi-implicitly,
with equal contribution from $\mathbf{u}_n$ and $\mathbf{u}^*$, as given below:
$$
\mathbf{u}^* - \Delta t\left[\frac{\nu\nabla^2\mathbf{u}^*}{2} \right ] = \mathbf{u}_{n} + \Delta t\left[\frac{\nu\nabla^2\mathbf{u}_n}{2} - \mathbf{u}_n.\nabla\mathbf{u}_n - \nabla p_n\right].
$$
The above equation has to be solved iteratively, and this is achieved through
OpenMP-parallelized Jacobi iterations.
The intermediate velocity field, $\mathbf{u}^*$, does not satisfy the continuity equation,
and requires appropriate correction.
This correction is obtained from the pressure correction term, which is calculated using the pressure Poisson equation,
$$
\nabla^2 p^* = \frac{\nabla.\mathbf{u}^*}{\Delta t}.
$$
``SARAS`` uses a Geometric Multigrid library to solve the above equation [@Wesseling:MG2004; @Briggs:MG2000].
Presently the library employs the Full Multigrid (FMG) V-Cycle to solve the Poisson equation.
Other methods like F-Cycle and W-Cycle are planned updates to the library in future.

Finally, using the above pressure correction, the velocity and pressure fields
corresponding to the next time-step are obtained as
$$
p_{n+1} = p_n + p^*,
$$
$$
\mathbf{u}_{n+1} = \mathbf{u}^* - \Delta t(\nabla p^*).
$$
The numerical implementation of the above procedure will be discussed in the next section.

# Numerical Method and Implementation
``SARAS`` uses finite-difference method [@Ferziger:book:CFD] to calculate derivatives of the field variables.
Presently the solver uses second-order central difference stencils to compute first and second derivatives.
We use a highly optimized and fast array manipulation library, Blitz++ [@Veldhuizen:CP1998], for all array operations.
Blitz++ also offers finite-difference stencils on uniformly spaced grids.
``SARAS`` augments this feature with grid-transformation terms [@Anderson:book:CFD] to perform simulations on grids with non-uniform spacing.

``SARAS`` offers a set of extensible boundary and initial conditions, as well as source terms.
Presently the solver can switch between periodic, Neumann, and Dirichlet boundary conditions.
There also exists a mixed boundary condition class that can simulate many practical applications,
such as a conducting heating plate on an adiabatic wall [@Teimurazov:2017].
Currently the solver is restricted to Cartesian grids, but it will be extended to cylindrical, toroidal, and spherical grids in the future.
The solver also supports adaptive time-stepping, where the Courant-Friedrichs-Lewy (CFL) condition
[@Courant:1928CFL] is used to dynamically compute the appropriate time-step.

# Results
We validate our code using two very well-known test problems: lid-driven cavity and decaying turbulence.
We simulate these problems using ``SARAS`` and compare the results with standard and validated solutions.

## Problem 1
We solve the two-dimensional lid-driven cavity (LDC) problem using ``SARAS``,
and compare the results with those of @Ghia:JCP1982.
LDC is an important fluid system and it serves as a benchmark for testing numerical methods.
This system consists of a square cavity of dimension $1 \times 1$
with no-slip boundary conditions on all the four walls.
However, the top wall moves laterally to the right with a constant velocity of $U = 1.0$
that serves as the reference velocity for non-dimensionalization of the problem.
The length of the side of the cavity, $L = 1.0$, is the reference length.
The Reynolds number of the flow is approximately $\mathrm{Re} \approx 1000$.

At the start of the simulation, the fluid, which is at rest, is driven impulsively by the top lid.
This results in the formation of a vortex at the upper-right corner of the cavity,
and this vortex rapidly grows in size and occupies the entire region of the cavity.
For the simulation we employ a $129 \times 129$ grid, and carry it out till $t = 30$.
The solver computes this solution using 4 MPI processes in approximately 12 minutes on an Intel workstation.
The output by ``SARAS`` at the final time is used for comparison with the results of @Ghia:JCP1982

In the website we supply the Bash and Python scripts to compile and run this test case.
After automatic execution of ``SARAS``, a Python script reads the output
from ``SARAS`` and compares the horizontal and vertical velocity profiles across
the geometric center of the square cavity.
The corresponding results from @Ghia:JCP1982 are also available with the installation.
We compare the two results and exhibit them in Figure \ref{figure1}.
We observe that the profiles computed by ``SARAS`` match very well with those by @Ghia:JCP1982, thus providing a strong validation for our code.
This quick validation of the solver, which can be done as a part of its installation, is one of the strengths of this package.

![Velocity profiles from the simulation of lid-driven cavity on a $129^2$ grid with ``SARAS`` (orange lines),
  plotted along with the data from @Ghia:JCP1982 (blue stars):
  (a) The vertical profile of the x-component of velocity, $u_x$, along the line across the geometric center of the cavity
  (b) The horizontal profile of the z-component of velocity, $u_z$, along the line across the geometric center of the cavity.
  \label{figure1}](resources/ldc_profiles.png)

## Problem 2
We simulate decaying turbulence using ``SARAS`` with Taylor-Green vortex [@Taylor:RSPA1937] as the initial condition.
That is,
$$
\mathbf{u}(x,y,z, t=0) = u_{0} \begin{bmatrix}
        \sin(2 \pi k_{0} x) \cos(2 \pi k_{0}y) \cos(2 \pi k_{0}z) \\
       -\cos(2 \pi k_{0} x) \sin(2 \pi k_{0}y) \cos(2 \pi k_{0}z) \\
        0
\end{bmatrix},
$$
where $u_0 = 1$ and $k_0 = 1$.
We perform our simulation in a periodic box of size $1 \times 1 \times 1$ ($L=1$) with a grid resolution of $257^{3}$ up to $t = 3.0$.
Here, we nondimensionalize time using $L/u_0$.
The initial Reynolds number of the flow is $\mathrm{Re} = 1000$.
We choose a constant $dt = 0.001$ for time-integration.
Besides, we use a uniform mesh along all the three directions. 

To validate the accuracy of the finite difference scheme of ``SARAS``,
we compare the results of ``SARAS`` with those of a pseudo-spectral code ``TARANG``,
which has been benchmarked and scaled up to 196608 cores of Cray XC40,
Shaheen II of KAUST [@Chatterjee:JPDC2018; @Verma:Pramana2013tarang].
Note that pseudo-spectral method yields very accurate derivatives [@Canuto:book:SpectralFluid].
We perform our spectral simulation using the same initial condition and grid resolution as above, up to $t = 3.0$.
The comparison of the results from the two codes is discussed below.
Comparison of ``SARAS`` results with those of ``TARANG`` provides another validation of ``SARAS``.

The results of ``TARANG`` and ``SARAS`` are quite similar, as exhibited in Figures \ref{figure2} and \ref{figure3}.
In Figure \ref{figure2}(a) we exhibit the plots of the total energy ($\int d{\bf r} u^2/2$) vs. time for both the runs.
As is evident from the plots, both the codes yield very similar evolution profiles for the total energy.
In addition, we also compute the energy spectra for the results at $t=1$.
As shown in Figure \ref{figure2}(b), both the energy spectra are very similar.
Interestingly, in both the plots, the energy spectra in the inertial range are close to the $k^{-5/3}$
power law [@Kolmogorov:DANS1941Dissipation; @Verma:book:ET].

![For the simulation of decaying turbulence on a $257^3$ grid with ``TARANG`` (thick blue lines) and ``SARAS`` (red-dashed lines):
  (a) plot of the total energy $E_u= \int d{\bf r} u^2/2$ vs $t$,
  (b) plot of $E_u(k)$ vs $k$ at $t =1$.
  \label{figure2}](resources/tgv_spectrum.png)

Figure \ref{figure3} exhibits the magnitude of velocity, $|\mathbf{u}|$, on the horizontal mid-plane ($z=1/2$) at $t=1$ and $t=3$.
Clearly, the results of ``TARANG`` and ``SARAS`` are very similar.

![For the simulation of decaying turbulence on a $257^3$ grid, vector plots of the velocity field,
  and the density plots of the magnitude of velocity ($|\mathbf{u}|$) computed at the horizontal mid-plane:
  for the data from ``TARANG``(a, c), and ``SARAS``(b, d) at $t = 1$ (top row) and $t = 3$ (bottom row).
  \label{figure3}](resources/tgv_velocity.png)


# Conclusions

This paper provides a brief description and validations tests of a new parallel finite-difference code ``SARAS``.
``SARAS`` has been designed as a general PDE solver using the object-oriented features of C++.
We validate ``SARAS`` using two test cases: lid-driven cavity and decaying turbulence.


# Acknowledgements

We gratefully acknowledge the contributions from Gaurav Gautham, Saurav Bhattacharjee, Rishabh Sahu,
and Fahad Anwer during the development of SARAS.

---

# References

