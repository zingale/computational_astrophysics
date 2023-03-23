# Going Further

## Cubic splines

Cubic splines are a popular way of interpolating data.  A cubic
polynomial is fit to 2 adjacent data samples.  Since we need 4
parameters to fit a cubic, we use the values of the function at the
two points as well as the derivatives.  And the derivatives of the
cubic at the end point are matched to those of the next cubic spline.

This means that all of the cubics are coupled together, and a linear
system needs to be solved for find all of the coefficients simultaneously.
