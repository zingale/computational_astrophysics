# Homework 1

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```


1. In class, we derived a finite-difference approximation to the
   first-derivative that was second-order accurate:

   $$f^\prime |_i = \frac{f_{i+1} - f_{i-1}}{2\Delta x}$$

   This uses one point to the left and right of point $i$.

   Now, instead, let's construct a *one-sided* difference for the
   first-derivative that is second-order accurate.  To do this,
   consider the points $f_i$, $f_{i+1}$, and $f_{i+2}$, and Taylor-expand
   them about $x_i$ and solve for $df/dx |_i$.

   Demonstrate that this convergeces second-order accurate by
   computing the derivative of $f(x) = \sin(x)$ at $x = 1$ for
   several values of $\Delta x$.
   
