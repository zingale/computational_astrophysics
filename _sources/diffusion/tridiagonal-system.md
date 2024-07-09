# Solving a Tridiagonal System

We want to solve a tridiagonal system of the form:
    
$$
\left ( \begin{array}{ccccccc}
   b_0 & c_0 \\
   a_1 & b_1 & c_1 \\
       & a_2 & b_2 & c_2 \\
       &     & \ddots & \ddots & \ddots \\
       &     &        & \ddots & \ddots & \ddots \\
       &     &        &        & a_{N-2} & b_{N-2} & c_{N-2} \\
       &     &        &        &         & a_{N-1} & b_{N=1} \\
\end{array}
\right )
\left (
\begin{array}{c} x_0 \\ x_1 \\ x_2 \\ \vdots \\ \vdots \\ x_{N-2} \\ x_{N-1} \\\end{array} \right )
  =
  \left ( \begin{array}{c} d_0 \\ d_1 \\ d_2 \\ \vdots \\ \vdots \\ d_{N-2} \\ d_{N-1} \\\end{array} \right )
$$

Here, the main diagonal is represented by $b_i$, which runs from $i =
0, \ldots, N-1$.  The two off diagonals have one fewer entries.

We can solve this just using Gaussian elimination, as before, but we
won't need to pivot, since the structure is the same on all rows.  The
result is the [tridiagonal matrix
solver](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm).
Notice that we only need to update 2 elements in each row instead of
$N$, so by taking advantage of this pattern, this solver will be much
faster than the full Gaussian elimination.


## Solving in python

In python, we can use [`linalg.solve_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html).  This works
by packed the diagonals into a 2D array, as shown on the pages there.

## Solving in C++

For C++, I've written a header-only set of functions that take vectors
of length $N$ for the three diagonals and the righthand side and
return $x$.  You can find these here:
https://github.com/zingale/computational_astrophysics/tree/main/examples/tridiagonal

You can run the test by compiling as:

```bash
g++ -I. -o tridiag tridiag.cpp
```
