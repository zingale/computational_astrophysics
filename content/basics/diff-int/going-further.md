Going Further
=============

There are many other types of integration methods that are useful.

Gaussian Quadrature
-------------------

A very powerful technique is [Gaussian
quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature) of
which Gauss-Legendre quadrature is one of the most popular techniques.
By using $N$ points, Gauss-Legendre can exactly integrate a polynomial
of order $2N - 1$!


Monte Carlo Integration
-----------------------

For multi-dimensional integration, we can apply our composite
integration rules in succession over each dimension.  But when the
dimension is really high ($N > 3$), this becomes computationally
expensive and randomly sampling the integrand becomes
better&mdash;this is [Monte Carlo
integration](https://en.wikipedia.org/wiki/Monte_Carlo_integration).
