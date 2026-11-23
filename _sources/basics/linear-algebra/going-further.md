Going Further
=============

Numerical linear algebra is a big topic.


Determinant
-----------

For a square matrix, we can compute the
[determinant](https://en.wikipedia.org/wiki/Determinant) easily once
the matrix is in row-echelon form, since it is an upper triangular matrix.
It is simply:

$$\mbox{det}\{{\bf A}\} = \prod_{i=1}^N A_{i,j}^\mbox{row-echelon}$$

where the superscript $\mbox{row-echelon}$ here means that ${\bf A}$ is in
row-echelon form.

There is one caveat: each time we swap rows when we pivot, we pick up a minus sign,
so the determinant becomes:

$$\mbox{det}\{{\bf A}\} = (-1)^{n_s} \prod_{i=1}^N A_{i,j}^\mbox{row-echelon}$$

where $n_s$ is the number of swaps.


Matrix Inverse
--------------

We can use Gaussian elimination to find the inverse of our matrix ${\bf A}$ by starting with

$${\bf A}{\bf A}^{-1} = {\bf I}$$

where ${\bf I}$ is the identity matrix.  Consider a single column, $j$ of ${\bf A}^{-1}$ which we'll call
${\bf x}^{(j)}$ and the corresponding column of ${\bf I}$ which we'll call ${\bf e}^{(j)}$.  Then we can
solve for that column of the inverse by solving:

$${\bf A} {\bf x}^{(j)} = {\bf e}^{(j)}$$

We can actually work on all of the columns together
to save on the computaitonal expense.

Sparse Matrices
---------------

For a matrix with a pattern of zeros, like a [tridiagonal
matrix](https://en.wikipedia.org/wiki/Tridiagonal_matrix) a linear
system can be solved with many fewer operations.  At tridiagonal or
banded matrix solver can be used in these cases.

In NumPy, ``numpy.tri()`` can work for this.


LU Decomposition
----------------

An alternative to Gaussian elimination is [LU
decomposition](https://en.wikipedia.org/wiki/LU_decomposition).  This
works factoring the matrix ${\bf A}$ as:

$${\bf A} = {\bf L}{\bf U}$$

where ${\bf L}$ is a lower triangle matrix and ${\bf U}$ is an upper triangle matrix.

LU decomposition is used when we need to solve the system ${\bf A}{\bf
x}$ with multiple different righthand sides.  Since most of the
computational expense is in doing the factorization, this can be much
less expensive when solving the same system with different righthand
sides.



Iterative Methods
-----------------

The methods that we saw here are _direct methods_, which give exact
solutions up to round-off error.  _Iterative methods_ are approximate
methods that solve the system to some tolerance.

These methods work well for linear systems that arise from the discretization
of an operator (like the Laplacian) and have the advantage that they
can be used without having to store the matrix itself.

We will see applications of these when we discuss the Poisson equation
and diffusion.

