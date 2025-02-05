# Homework 1

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```


1. Consider the function:

   $$f(x) = \frac{1}{\sqrt{x^3 + 1} - 1}$$

   defined for $x > 0$.  When $x$ is small, roundoff error
   from the subtraction will dominate the answer.  We want
   to see if we can get a better result.

   * Show that by multiplying and dividing $f(x)$ by
     $\sqrt{x^3 + 1} + 1$, you get an analytically
     equivalent expression without any subtractions.
     Call this new expression $g(x)$.

   * Using the fact that $\sqrt{x^3 + 1} \approx 1 + x^3/2$ when $|x| \ll 1$,
     show that you get:

     $$f(x) \approx \frac{2}{x^3}$$

   * Write a program to evaluate $f(x)$, $g(x)$, and the approximation
     to $f(x)$ for $|x| \ll 1$ for the values of $x$: `1.e-4`, `1.e-5`, and `1.e-6`.

     What do you observe?

2. In class, we derived a finite-difference approximation to the
   first-derivative that was second-order accurate:

   $$f^\prime |_i = \frac{f_{i+1} - f_{i-1}}{2\Delta x}$$

   This uses one point to the left and right of point $i$.

   Now, instead, let's construct a *one-sided* difference for the
   first-derivative that is second-order accurate.  To do this,
   consider the points $f_i$, $f_{i+1}$, and $f_{i+2}$, and Taylor-expand
   them about $x_i$ and solve for $df/dx |_i$.

   Demonstrate that this converges second-order accurate by
   computing the derivative of $f(x) = \sin(x)$ at $x = 1$ for
   several values of $\Delta x$.
