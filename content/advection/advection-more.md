# More on Advection and Burgers' Equation

## Piecewise Parabolic Reconstruction

In astrophysics, the [piecewise parabolic method (PPM)](https://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/A141984.pdf) is widely used.  This
reconstructs the data in each cell as parabola (instead of piecewise constant for our first-order method or piecewise linear for our second order method).  As a result, this is much more accurate.

Here's an example code that implements PPM for advection and shows some visualizations: [PPMpy](https://python-hydro.github.io/ppmpy/intro.html)

````{tip}
You can install PPMpy via `pip` as:

```
pip install ppmpy
```
````

An advection solver can be created using this, as shown here:
https://python-hydro.github.io/ppmpy/advection.html

## Multi-dimensions

In two-dimensions, the advection equation is:

$$\frac{\partial a}{\partial t} + u \frac{\partial a}{\partial x} + v \frac{\partial a}{\partial y} = 0$$

or in conservative form:

$$\frac{\partial a}{\partial t} + \frac{\partial (ua)}{\partial x} +
                                  \frac{\partial (va)}{\partial y} = 0$$

defining $F^{(x)}(a) = (ua)$ and $F^{(y)}(a) = (va)$, this is:

$$\frac{\partial a}{\partial t} + \frac{\partial F^{(x)}(a)}{\partial x} +
                                  \frac{\partial F^{(y)}(a)}{\partial y} = 0$$

Now defining the average in two-dimensions as:

$$\langle a \rangle_{i,j} = \int_{x_{i-1/2}}^{x_{i+1/2}} \int_{y_{j-1/2}}^{y_{j+1/2}} a(x, y) dx dy$$

and integrating the conservative form of the advection equation over a cell,
we have:

$$
\begin{align*}\frac{\partial \langle a\rangle_{i,j}}{\partial t} =
   &-\frac{1}{\Delta x \Delta y} \int_{x_{i-1/2}}^{x_{i+1/2}} \int_{y_{j-1/2}} \frac{\partial F^{(x)}(a)}{\partial x} dx dy \\
   &-\frac{1}{\Delta x \Delta y} \int_{x_{i-1/2}}^{x_{i+1/2}} \int_{y_{j-1/2}} \frac{\partial F^{(y)}(a)}{\partial y} dx dy
\end{align*}
$$

Now for each of those double integral terms, the integral in the normal direction
is just the divergence theorem, and means that we just evaluate the flux at the 
limits of integration in that direction.  The other integral represents an
average of the flux over the interface in the transverse direction.  We can
just represent that average as the flux evaluated in the center (in the transverse direction), giving:

$$\begin{align*}
\frac{\partial \langle a\rangle_{i,j}}{\partial t} =
   &- \frac{[F^{(x)}(a)]_{i+1/2,j} - [F^{(x)}(a)]_{i-1/2,j}}{\Delta x} \\
   &- \frac{[F^{(y)}(a)]_{i,j+1/2} - [F^{(y)}(a)]_{i,j-1/2}}{\Delta y} 
\end{align*}$$

This can be solved in pretty much the same fashion as we did in one-dimension,
but now our grid will need to be two-dimensional and we would have to compute
fluxes and upwind in both directions.

The [pyro code](https://python-hydro.github.io/pyro2/) solves advection
(and more) in two-dimensions.

````{tip}
You can install pyro via `pip` as:

```
pip install pyro-hydro
```
````

Some discussion of the advection solver is here:
https://python-hydro.github.io/pyro2/advection_basics.html


