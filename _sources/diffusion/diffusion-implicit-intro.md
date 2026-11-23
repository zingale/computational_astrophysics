# Implicit Methods for Diffusion

One of the issues we have with our diffusion implementation is that the timestep constraint:

$$\Delta t = \frac{C}{2} \frac{\Delta x^2}{k}$$

can be very restrictive as $\Delta x$ gets small.  

```{note}
In astrophysical simulations, we are often combining advection and diffusion,
and the diffusive timestep restriction can dominate, making our overall simulation
computationally expensive.

For this reason, we often wish to use an _implicit-in-time_
integration method for diffusion.
```

Going back to our original discretization, we would write the implicit update as:

$$\frac{\phi_i^{n+1} - \phi_i^n}{\Delta t}
 = k \frac{\phi_{i+1}^{n+1} - 2 \phi_i^{n+1} + \phi_{i-1}^{n+1}}{\Delta x^2}$$
 
If we define

$$\alpha = \frac{k\Delta t}{\Delta x^2}$$

then we can write this as:

$$-\alpha \phi^{n+1}_{i-1} + (1 + 2\alpha) \phi_i^{n+1}
- \alpha \phi_{i+1}^{n+1} = \phi_i^n$$

This is a coupled set of algebraic equations for the new time state of $\phi$&mdash;a linear system.

```{important}
The implicit discretization is stable for all values of $C$.  But keep in mind that
stability and accuracy are not the same thing.
```

There are several ways we can solve this:

* _direct solve_ : create a linear system of the form ${\bf A}{\bf x} = {\bf b}$
  and solve for the new state using Gaussian elimination or something similar
  
* _iterative solve_ : apply relaxation to the system, just like we did with our
  elliptic equations.
  

## Direct solve

In order to write down the linear system, we need to understand our boundary
conditions.  Consider the left boundary.  We'll denote the cell just inside
the domain as $\phi_\mathrm{lo}$ and the ghost cell just outside of the domain
as $\phi_{\mathrm{lo}-1}$.

* _Homeogeneous Neumann_ : for Neumann boundary conditions, we want the
  gradient at the boundary to be zero.  This means:
  
  $$\frac{\phi_\mathrm{lo} - \phi_{\mathrm{lo}-1}}{\Delta x} = 0$$
  
  This is second-order accurate at the physical boundary, $x_{\mathrm{lo}-1/2}$.
  
  Inserting this into our implicit system above, we find that the update for
  the first zone inside the boundary is:
  
  $$(1 + \alpha) \phi_\mathrm{lo}^{n+1} - \alpha \phi_{\mathrm{lo}+1}^{n+1} = \phi_\mathrm{lo}^n$$
  
* _Dirichlet_ : for Dirichlet boundary conditions, we want the value at the
  physical boundary to be specified&mdash;we'll call it $\phi_A$.  To second-order, that means that the
  average of the two zones on either side of the boundary should give us
  the desired value:
  
  $$\frac{\phi_\mathrm{lo} + \phi_{\mathrm{lo}-1}}{2} = \phi_A$$
  
  Inserting this into our implicit update, we find that the update for the first
  zone inside the boundary is:
  
  $$(1 + 3\alpha) \phi_\mathrm{lo}^{n+1} - \alpha \phi_{\mathrm{lo}+1}^{n+1} = \phi_\mathrm{lo}^n + 2 \alpha \phi_A$$
  
Similar constructions are used for the right boundary.

Putting this together, our linear system for implicit diffusion (assuming Neumann boundary conditions) appears as:

$$
\left (
\begin{array}{ccccccc}
1+\alpha &   -\alpha &           &        &         &           &          \\
-\alpha  & 1+2\alpha & -\alpha   &        &         &           &          \\
         & -\alpha   & 1+2\alpha & -\alpha&         &           &          \\
         &           & \ddots    & \ddots & \ddots  &           &          \\
         &           &           & \ddots & \ddots  & \ddots    &          \\
         &           &           &        & -\alpha & 1+2\alpha &-\alpha   \\
         &           &           &        &         & -\alpha   &1+\alpha \\
\end{array}
\right )
\left (
\begin{array}{c}
\phi_\mathrm{lo}^{n+1} \\
\phi_\mathrm{lo+1}^{n+1} \\
\phi_\mathrm{lo+2}^{n+1} \\
\vdots \\
\vdots \\
\phi_\mathrm{hi-1}^{n+1} \\
\phi_\mathrm{hi}^{n+1} \\
\end{array}
\right )
=
\left (
\begin{array}{c}
\phi_\mathrm{lo}^{n} \\
\phi_\mathrm{lo+1}^{n} \\
\phi_\mathrm{lo+2}^{n} \\
\vdots \\
\vdots \\
\phi_\mathrm{hi-1}^{n} \\
\phi_\mathrm{hi}^{n} \\
\end{array}
\right )
$$

We note a few things:

* This is a [tridiagonal matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix).
  Efficient solution methods exist for this that take advantage of the [sparsity pattern](https://en.wikipedia.org/wiki/Sparse_matrix).
  
* The ghost cells in our grid do not explicitly appear in the system

* The matrix is [diagonally dominant](https://en.wikipedia.org/wiki/Diagonally_dominant_matrix).  This means that
  it will work well with iterative methods.