#ifndef GAUSS_H
#define GAUSS_H

#include <algorithm>
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>

#include "matrix.H"

inline
std::vector<double> gauss_elim(Matrix& A, std::vector<double>& b) {

    // perform gaussian elimination with pivoting, solving A x = b.
    //
    //  A is an NxN matrix, x and b are an N-element vectors.  Note: A and b are
    //  changed upon exit to be in upper triangular (row echelon) form

    auto N = b.size();

    // A is square, with each dimension of length N

    assert(A.nrows() == N && A.ncols() == N);

    // allocate the solution array

    std::vector<double> x(N, 0.0);

    // main loop over rows

    for (std::size_t krow = 0; krow < N; ++krow) {

        // find the pivot row based on the size of column k -- only
        // consider the rows beyond the current row

        auto row_max{krow};
        double col_max = std::numeric_limits<double>::lowest();
        for (std::size_t kk = krow; kk < N; ++kk) {
            if (std::abs(A(kk, krow)) > col_max) {
                col_max = std::abs(A(kk, krow));
                row_max = kk;
            }
        }

        // swap the row with the largest element in the current column
        // with the current row (pivot) -- do this with b too!

        if (row_max != krow) {

            for (std::size_t jcol = 0; jcol < N; ++jcol) {
                std::swap(A(krow, jcol), A(row_max, jcol));
            }
            std::swap(b[krow], b[row_max]);
        }

        // do the forward-elimination for all rows below the current

        for (std::size_t irow = krow+1; irow < N; ++irow) {
            double coeff = A(irow, krow) / A(krow, krow);

            for (std::size_t jcol = krow+1; jcol < N; ++jcol) {
                A(irow, jcol) += -A(krow, jcol) * coeff;
            }

            A(irow, krow) = 0.0;
            b[irow] += -coeff * b[krow];
        }
    }

    // back-substitution

    x[N-1] = b[N-1] / A(N-1, N-1);

    // size_t is unsigned, so it is never < 0, it will
    // wrap to a large number
    for (std::size_t irow = N-2; irow != std::numeric_limits<std::size_t>::max(); --irow) {
        std::cout << irow << std::endl;
        double bsum = b[irow];
        for (std::size_t jcol = irow+1; jcol < N; ++jcol) {
            bsum += -A(irow, jcol) * x[jcol];
        }
        x[irow] = bsum / A(irow, irow);
    }

    return x;
}

#endif
