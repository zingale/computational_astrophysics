#ifndef INVERSE_H
#define INVERSE_H

#include <cassert>
#include <limits>

#include "matrix.H"

inline
Matrix inverse(const Matrix& A_in) {

    // copy our array
    auto A = A_in;

    const auto N = A.nrows();

    assert(A.nrows() == A.ncols());

    // allocate the inverse array
    Matrix Ainv(N, N);

    // create the identity
    Matrix I(N, N);
    for (std::size_t i = 0; i < N; ++i) {
        I(i, i) = 1.0;
    }

    // we will do Gaussian elimination but with A Ainv = I, so
    // each column of I represents a separate righthand side to
    // and Ax = b linear system.

    for (std::size_t krow = 0; krow < N; ++krow) {

        // find the pivot row based on the size of column k -- only consider the
        // rows beyond the current row

        std::size_t row_max{krow};
        double col_max = std::numeric_limits<double>::lowest();
        for (int kk = krow; kk < N; ++kk) {
            if (std::abs(A(kk, krow)) > col_max) {
                col_max = std::abs(A(kk, krow));
                row_max = kk;
            }
        }

        // swap the row with the largest element in the current column
        // with the current row (pivot) -- do this with all columns of
        // I too!

        if (row_max != krow) {
            for (int jcol = 0; jcol < N; ++jcol) {
                std::swap(A(krow, jcol), A(row_max, jcol));
                std::swap(I(krow, jcol), I(row_max, jcol));
            }
        }

        // do the forward-elimination for all rows below the current

        for (int irow = krow+1; irow < N; ++irow) {
            double coeff = A(irow, krow) / A(krow, krow);

            for (int jcol = krow+1; jcol < N; ++jcol) {
                A(irow, jcol) += -A(krow, jcol) * coeff;
            }

            A(irow, krow) = 0.0;
            for (int jcol = 0; jcol < N; ++jcol) {
                I(irow, jcol) += -coeff * I(krow, jcol);
            }
        }

    }

    // back substituion -- once for each column c in the matrix I
    for (std::size_t c = 0; c < N; ++c) {

    // back-substitution

        Ainv(N-1, c) = I(N-1, c) / A(N-1, N-1);

        for (std::size_t irow = N-2; irow != std::numeric_limits<std::size_t>::max(); --irow) {
            double bsum = I(irow, c);
            for (std::size_t jcol = irow+1; jcol < N; ++jcol) {
                bsum += -A(irow, jcol) * Ainv(jcol, c);
            }
            Ainv(irow, c) = bsum / A(irow, irow);
        }

    }

    return Ainv;
}


#endif
