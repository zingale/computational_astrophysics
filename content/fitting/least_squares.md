Least Squares Regression
========================

Imagine a dataset where we have a set of points, $(x_i, y_i)$ with
associated errors in $y$, $\sigma_i$, and we want to fit a line or
curve to them.

Our data looks like:

```{figure} ./fitting-illustration.png
---
width: 100%
---
Set of points and associated errors
```

We want to create a model, $Y(x; {a_j})$ that best fits these points.
Here the values ${a_j}$ are free parameters that we will tune to get a
good fit.

For each point, we will measure the _vertical_ distance to the model:

$$\Delta_i = Y(x_i; {a_j}) - y_i$$

and we want to minimize the distance.  We do this by defining

$$\chi^2({a_j}) \equiv \sum_i \left(\frac{\Delta_i}{\sigma_i}\right )^2$$

Here, each point's distance, $\Delta_i$ is weighted by its error,
$\sigma_i$ &mdash; this ensures that points that we are least certain
of has less influence in the fit.  Minimizing $\chi^2$ with respect to
the fit parameters, ${a_j}$ is called _least squares minimization_.


