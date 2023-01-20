Gaussian Elimination
====================

We are now ready to write out the algorithm for Gaussian elimination.
The basic idea is to get the matrix ${\bf A}$ into row-echelon form.
If the matrix is square, and non-singular ($det{\bf A} \ne 0$), then a
row-echelon matrix is an upper-trangular matrix.  The process of
transforming the matrix this way is called _forward-elimination_.

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
