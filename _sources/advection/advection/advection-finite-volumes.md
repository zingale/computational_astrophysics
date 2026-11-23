# Finite-Volume Discretization

Ultimately we are interested in conservation laws, so we want a method
that respects conservation.  We'll look at this with advection first.

Here we see 3 different types of grids:

![discretization types](grids.png)

The first 2 are *finite-difference* grids&mdash;the data is
represented at a specific point in the domain.  The two differ in that
the second case is a cell-centered finite-difference grid.  In this
case, there is not a point on the physical boundaries, while with the
first example (what we were using previously), there is a data point
on each boundary.  For the cell-centered case, we typically talk about
dividing the data into cells or zones and the data living at the
center of the zone.

The last example is a *finite-volume* grid.  Here, we don't think of
the data living at a specific point, but instead we keep track of the
total amount of a quantity (or its average) in a volume.  This is
represented above as the shaded region inside our zone.

```{tip}
We often use the terms _zone_, _cell_, and _volume_ interchangeably.
```

```{note}
We use half-integer indices to denote the location of the boundaries of
a volume, so:

$$x_{i-1/2} = x_i - \frac{\Delta x}{2}$$

and the volume represented by $i$ extends from $[x_{i-1/2}, x_{i+1/2}]$.
```

## Cell-averages

Imagine we have a function $f(x)$ that we wish to represent on our grid, using a finite-volume discretization.

![finite-volume grid](fv_grid.png)

We can define the average of $f(x)$ in a zone that goes from
$x_{i-1/2}$ to $x_{i+1/2}$ as:

$$\langle f\rangle_i = \frac{1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}}f(x) dx$$

We use the angle brackets to indicate that this is an average, and use
an integer index, $i$, to indicate that this data lives in zone $i$.

Now consider our linear advection equation,

$$\frac{\partial a}{\partial t} = -u\frac{\partial a}{\partial x}$$

written in conservative form:

$$\frac{\partial a}{\partial t} = -\frac{\partial F(a)}{\partial x}$$

where $F(a) = ua$ is the flux.

Let's average (integrate) this equation over a zone (from $[x_{i-1/2},
x_{i+1/2}]$):

$$\frac{1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} \frac{\partial a}{\partial t} dx =
   -\frac{1}{\Delta x} \int_{x_{i-1/2}}^{x_{i+1/2}} \frac{\partial F}{\partial x} dx$$

we can recognize the left is the time-derivative of the average of $a$
and the right, via the divergence theorem is just the flux different
through the boundary of the zone:

$$\frac{\partial}{\partial t} \langle a\rangle_i = -\frac{1}{\Delta x}
   \left [ F_{i+1/2} - F_{i-1/2} \right ]$$

Here, $F_{i-1/2}$ indicates the flux moving through the interface at location $x_{i-1/2}$.

```{note}
Finite-volume method are conservative because the flux leaving one volume enters
the adjacent zone.

For example, in the updates for zones $i-1$ and $i$, we have:

$$\frac{\partial}{\partial t} \langle a\rangle_{i-1} = -\frac{1}{\Delta x}
   \left [ F_{i-1/2} - F_{i-3/2} \right ]$$

$$\frac{\partial}{\partial t} \langle a\rangle_i = -\frac{1}{\Delta x}
   \left [ F_{i+1/2} - F_{i-1/2} \right ]$$

So any amount of $\langle a\rangle_{i-1}$ that leaves zone $i-1$ through the right interface is gained
by $\langle a\rangle_i$ by that same flux entering through its left interface.  This is an
example of a [telescoping sum](https://en.wikipedia.org/wiki/Telescoping_series).
```

## Time update

Let's start with a first-order update in time.  We can use the Euler
method we already saw,

$$\frac{\partial \langle a \rangle}{\partial t} \approx \frac{\langle a\rangle^{n+1} - \langle a\rangle^{n}}{\Delta t}$$

in which case our update appears as:

$$\langle a\rangle_i^{n+1} = \langle a\rangle_i^{n} - \frac{\Delta t}{\Delta x} (F_{i+1/2}^n - F_{i-1/2}^n)$$

We expect this to be first order accurate in time.  Notice that the
fluxes are evaluated at the old time.


## Reconstruction and the Riemann problem

We need the value of the flux on the interface, $F_{i+1/2}$, we find this using the analytic expression
for the flux, $F(a)$, as:

$$[F(a)]_{i+1/2} = F(a_{i+1/2})$$

So now we need to find the interface state, $a_{i+1/2}$.  

Getting the value of $a$ on the interface from the average, $\langle a
\rangle$ is called *reconstruction*.  It means we need to infer how
$a$ actually varies throughout the cell just from the information we
have about its average value.  There are a variety of methods we can
use (some of which we will explore later).  For now, we will do the
simplest, and assume that $a(x)$ is constant in each cell:

$$a_{i+1/2,L} = \langle a\rangle_i^n$$
$$a_{i+1/2,R} = \langle a\rangle_{i+1}^n$$

Notice that we actually have 2 values for the interface state when we
do this, one coming from each side of the interface&mdash;we label
these as the left and right states:

![left and right states for advection](riemann-adv-mol.png)

Now we need to resolve this degeneracy&mdash;which of the two states
(or what combination of them) is the correct interface state?  This is called the
*Riemann problem*.

```{note}
The Riemann problem is where the physics of the actual set of equations
comes into play.  Thus far, everything we've done has been general.
```

For advection, the Riemann problem is easy.  We know that for 
$u > 0$ that $a(x)$ moves
from left to right, so the correct state on the interface is the left
state&mdash;this is upwinding.

$$
a_{i+1/2} = \mathcal{R}(a_{i+1/2,L}, a_{i+1/2,R}) = 
\left \{ 
\begin{array}{c}
 a_{i+1/2,L} \quad u > 0 \\
 a_{i+1/2,R} \quad u < 0 
\end{array}
\right .
$$

where we indicate the Riemann problem as
$\mathcal{R}(a_{i+1/2,L},a_{i+1/2,R})$.

For the advection equation, with $u > 0$, our interface states are 

$$F_{i+1/2} = F(a_{i+1/2}) = u a_{i+1/2} = u \langle a \rangle_i$$

$$F_{i-1/2} = F(a_{i-1/2}) = u a_{i-1/2} = u \langle a \rangle_{i-1}$$

Inserting these into our difference equation, we have:

$$\langle a \rangle_i^{n+1} = \langle a \rangle_i^n - \frac{u \Delta t}{\Delta x} (\langle a \rangle_i^n - \langle a \rangle_{i-1}^n)$$

This is precisely the upwind finite-difference scheme we saw earlier.
