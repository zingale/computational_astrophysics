# routines for multiplying matrices

import numpy as np

def mult_Ax(A, x):
    """ return the product of matrix A and vector x: Ax = b """

    # x is a vector
    assert x.ndim == 1, "x should be a vector"

    N = len(x)

    # A is square, with each dimension of length N
    assert A.shape == (N, N), "A should be square with each dim of same length as x"

    # allocation the product array
    b = np.zeros((N), dtype=A.dtype)

    # each row of b is the product of the like row of A dotted with
    # the vector x
    for i in range(N):
        b[i] = A[i, :] @ x

    return b
