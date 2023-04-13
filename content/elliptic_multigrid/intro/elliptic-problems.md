# Poisson problem

Consider the Poisson equation

$$\nabla^2 \phi = f$$

This is a second-order elliptic equation, and therefore requires 2 boundary conditions.

There is no time-dependence in this equation.  The potential $\phi$
responds instantaneously to the source $f$ and the boundary
conditions.

Consider the one-dimensional Poisson equation:

$$\phi^{\prime\prime} = f$$

on the domain $[a, b]$.

At each boundary. we can specify the boundary conditions as

* [Dirichlet](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition): $\phi(a) = A$
* [Neumann](https://en.wikipedia.org/wiki/Neumann_boundary_condition): $\phi^\prime(a) = C$
* [Robin](https://en.wikipedia.org/wiki/Robin_boundary_condition): $\alpha \phi(a) + \beta \phi^\prime(a) = \gamma$ (a mix of the above two)
* periodic

We can also have diffrent boundary conditions on each boundary.
If the values are set to 0, we call the
conditions homogeneous, otherwise we call them inhomogeneous.

```{note}
Not any set of boundary conditions is allowed.
Consider $f = 0$, so our Poisson equation (is the Laplace equation):

$$\phi^{\prime\prime} = 0$$

and the solution is just

$$\phi(x) = a x + b$$

if we try to enforce different inhomogeneous Neumann boundary
conditions on each end, then we get conflicting values for the
slope&mdash;this is unsolvable.
```

To understand solvable boundary conditions for the general case, 

$$\nabla^2 \phi = f$$

we start by integrating over the domain

$$\int_\Omega \nabla^2 \phi d\Omega = \int_{\partial \Omega} \nabla \phi \cdot {\bf n} dS = \int_\Omega f d\Omega$$

If we have homogeneous Neumann BCs on all sides, $\nabla \phi \cdot
{\bf n} = 0$, then the source, $f$, must satisfy

$$\int_\Omega f d\Omega = 0$$

The same condition will apply, if the boundary conditions are periodic.
