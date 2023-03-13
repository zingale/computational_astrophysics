Interpolation
=============

We frequently encounter situations where we only have data at a
discrete number of points (or it is to expensive to compute it _on
demand_, so we tabulate it).

Interpolation fills in the gaps by making an assumption about the
behavior of the underlying functional form of the data.

Many types of interpolation exists:

* Some ensure that no new extrema a introduced.

* Some ensure a physical constraint (like thermodynamic consistency) is satisfied/

* Some ensure smoothness of the fit.

* Some ensure the quantity being interpolated is conserved.

Generally speaking, you need to balance the number of points used in
constructing an interpolant (which _can_ increase accuracy) against
pathologies (like oscillations).

```{note}
_Interpolation_ and _Fitting_ are different operations.
Fitting seeks to produce a simple functional model that represents the
entire dataset.  Interpolation looks to fill in the gaps in some
region of your dataset.
```

References:

* Pang, _An Introduction to Computational Physics_ : this provides a
  nice discussion.

