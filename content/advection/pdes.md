# Partial Differential Equations

We will now look at partial differential equations.  There are [3
broad classes of
PDEs](https://en.wikiversity.org/wiki/Partial_differential_equations#Elliptic,_Hyperbolic,_and_Parabolic_PDEs)
that we will study: hyperbolic, elliptic, and parabolic PDEs.  These
will each require different solution methods.

```{seealso}
For more in-depth discussions of the methods used for PDEs in astrophysics,
you can read my notes here: [_Introduction to Computational Astrophysical Hydrodynamics_](http://bender.astro.sunysb.edu/hydro_by_example/CompHydroTutorial.pdf).
```

## Hyperbolic PDEs

The canonical hyperbolic PDE is the wave equation:

$$\frac{\partial^2 \phi}{\partial t^2} = c^2 \frac{\partial^2 \phi}{\partial x^2}$$

The solution to this is traveling waves moving at speed $c$ in both directions:

$$\phi(x, t) = A f_0(x - ct) + B g_0(x + ct)$$

A simpler example of this that we will look at is the _linear advection equation_:

$$\frac{\partial a}{\partial t} + u \frac{\partial a}{\partial x} = 0$$

The defining feature of hyperbolic PDEs is that there is a finite, real speed at which the solution changes.

In astrophysics, advection-like equations describe fluid dynamics.

## Elliptic PDEs

The Poisson equation is the canonical elliptic PDE:

$$\nabla^2 \phi = f$$

Notice that there is no time variable here.  Information about the solution
moves at infinite speed, and the solution depends only on the boundary conditions and the source, $f$.

In astrophysics, the Poisson equation is used to obtain the gravitational
potential, via [_Poisson's equation_](https://en.wikipedia.org/wiki/Poisson%27s_equation#Newtonian_gravity)

$$\nabla^2 \Phi = 4 \pi G \rho$$

## Parabolic PDEs

The heat equation is the canonical parabolic PDE:

$$\frac{\partial \phi}{\partial t} = k \frac{\partial^2 \phi}{\partial x^2}$$

As with elliptic PDEs, information about the solution is communicated instantaneously&mdash;there is no finite propagate speed like the hyperbolic case.
But this is time-dependent.  

Diffusion of photons or heat in stars is perhaps the main area where parabolic PDEs
arise in astrophysics.

