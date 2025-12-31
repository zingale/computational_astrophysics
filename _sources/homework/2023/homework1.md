# Homework 1

```{note}
You are free to discuss these questions with your classmates and on
our class slack, but you must write your own solutions, including your
own source code.

All code should be uploaded to Brightspace along with any analytic
derivations, notes, etc.
```


1. Following the example we did in class to find a difference
   approximation to the first derivative of $f(x)$, find an
   approximation to the _second derivative_ of $f(x)$.

   You should consider $f(x-h)$, $f(x)$, and $f(x+h)$.

   * What is the truncation error of this approximation?

   * Write a program to compute the second derivative of $f(x) =
     \sin(x)$ as a function of $h$ and print the absolute value of the
     error in your estimate compared to the analytic solution for
     several values of $h$.

   * Does the error converge as you expect?

2. Consider the function:

   $$f(x) = \sqrt{x^2 + 1} - 1$$

   When $x$ is small, roundoff error from the subtraction will
   dominate the answer.  In the limit of $x \rightarrow 0$, this
   expression takes the form:

   $$f(x) \sim \frac{1}{2} x^2$$

   * Show that by multiplying and dividing $f(x)$ by $\sqrt{x^2 + 1} +
     1$ you get an analytically equivalent expression for $f(x)$
     without any subtractions.

   * Write a program to evaluate $f(x)$ and the new version of $f(x)$
     for $x$ with the values `1.e-6`, `1.e-7`, `1.e-8`.

   * Which version is more accurate?
