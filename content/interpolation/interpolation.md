Interpolation
=============

We frequently encounter situations where we only have data at a
discrete number of points, or a function is very to expensive to
compute _on demand_, so we tabulate it at several points.

[Interpolation](https://en.wikipedia.org/wiki/Interpolation) fills in the gaps in our data by making an assumption about the
behavior of the underlying functional form of the data.

Many types of interpolation exists:

* Some ensure that no new extrema a introduced.

* Some ensure a physical constraint (like thermodynamic consistency) is satisfied/

* Some ensure smoothness of the fit.

* Some ensure the quantity being interpolated is conserved.

Generally speaking, you need to balance the number of points used in
constructing an interpolant (which _can_ increase accuracy) against
pathologies (like oscillations).

```{important}
_Interpolation_ and _Fitting_ are different operations.
Fitting seeks to produce a simple functional model that represents the
entire dataset.  Interpolation looks to fill in the gaps in some
region of your dataset.
```

```{admonition} Reference
The discussion in Pang, _An Introduction to Computational Physics_ provides a
good introduction.
```


