Gaussian Elimination
====================

We are now ready to write out the algorithm for Gaussian elimination.
The basic idea is to get the matrix ${\bf A}$ into row-echelon form.
If the matrix is square, and non-singular ($\mbox{det}\{\bf A\} \ne 0$), then a
row-echelon matrix is an upper-trangular matrix.  The process of
transforming the matrix this way is called _forward-elimination_.
Then we solve for ${\bf x}$ by doing _back substitution_.

Forward Elimination
-------------------

Consider the following state of our matrix ${\bf A}$ mid-way through the
row-echelon transformation:

$$
\left (
\begin{array}{ccccccc}
  a_{1,1} & a_{1,2} & \cdots & a_{1,k} & a_{1,k+1} & \cdots & a_{1,N} \\
          & a_{2,2} & \cdots & a_{2,k} & a_{2,k+1} & \cdots & a_{2,N} \\
          &         & \ddots  \\
          &         &        & a_{k,k} & a_{k,k+1} & \cdots & a_{k,N} \\
          &         &        & a_{k+1,k} & a_{k+1,k+1} & \cdots & a_{k+1,N} \\
          &         &        & \vdots    & \vdots      &        & \vdots \\
          &         &        & a_{N,k} & a_{N,k+1} & \cdots & a_{N,N} \\
\end{array}
\right )
$$

At this point, the matrix is in row-echelon form down to row $k$, but
the next row, $k+1$ is not, since there is a non-zero element below
the diagonal ($a_{k+1,k}$).

Our next step in forward-elimination
is to operate on all rows $k+1, \ldots, N$ using row $k$ to eliminate
any non-zero entries in column $k$.

The procedure is:

1. Define a coefficient

   $$f_{k+1} = \frac{a_{k+1,k}}{a_{k,k}}$$

2. subtract $f_{k+1} \times \{ \mbox{row}~k \}$
   from row $k+1$

3. correct the $k+1$ row of the righthand size vector, ${\bf b}$
   with the same factor

4. Repeat the procedure for all rows $k+2, \ldots, N$


Back Substitution
-----------------

At the end of forward elimination, our system looks like:

$$
\left (
\begin{array}{ccccccc}
  a_{1,1} & a_{1,2} & \cdots & a_{1,k} & \cdots   & a_{1,N-1} & a_{1,N} \\
          & a_{2,2} & \cdots & a_{2,k} & \cdots   & a_{2,N-1} & a_{2,N} \\
          &         & \ddots &         &          &           & \vdots  \\
          &         &        & a_{k,k} & \cdots   & a_{k,N-1} & a_{k,N} \\
          &         &        &         & \ddots   &           & \vdots \\
          &         &        &           &        & a_{N-1,N-1} & a_{N-1,N} \\
          &         &        &           &        &             & a_{N,N} \\
\end{array}
\right )
\left (
\begin{array}{c}
 x_1 \\ x_2 \\ \vdots \\ x_k \\ \vdots \\ x_{N-1} \\ x_N
\end{array}
\right )
=
\left (
\begin{array}{c}
 b_1 \\ b_2 \\ \vdots \\ b_k \\ \vdots \\ b_{N-1} \\ b_N
\end{array}
\right )
$$

The last element is simply

$$x_N = \frac{b_N}{a_{N,N}}$$

The second-to-last element is then found by solving:

$$a_{N-1,N-1} x_{N-1} + a_{N-1,N}\underbrace{x_N}_{\mbox{known}} = b_{N-1}$$

giving:

$$x_{N-1} = \frac{b_{N-1} - a_{N-1,N} x_N}{a_{N-1,N-1}}$$

and so forth, with the general case looking like:

$$x_k = \frac{b_k - \sum_{j=k+1}^N a_{k,j} x_j}{a_{k,k}}$$

