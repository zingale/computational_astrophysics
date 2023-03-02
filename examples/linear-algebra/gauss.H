#ifndef GAUSS_H
#define GAUSS_H

#include <array.H>
#include <vector>
#include <limits>
#include <iostream>
#include <iomanip>

inline
void print_Ab(const Array& A, const std::vector<double>& b) {

    // printout the matrix A and vector b in a pretty fashion.  We don't use the
    // numpy print here, because we want to make them side by side

    const int N = b.size();

    for (int irow = 0; irow < N; ++irow) {

        std::cout.precision(4);
        std::cout << std::fixed;

        if (irow == 0) {
            std::cout << "/ ";
        } else if (irow == N-1) {
            std::cout << "\\ ";
        } else {
            std::cout << "| ";
        }

        for (int jcol = 0; jcol < N; ++jcol) {
            std::cout << std::setw(8) << A(irow, jcol);
        }

        std::cout << " | ";
        std::cout << std::setw(8) << b[irow];

        if (irow == 0) {
            std::cout << " \\ ";
        } else if (irow == N-1) {
            std::cout << " / ";
        } else {
            std::cout << " | ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

inline
std::vector<double> gauss_elim(Array& A, std::vector<double>& b,
                               const bool quiet=false) {

    // perform gaussian elimination with pivoting, solving A x = b.
    //
    //  A is an NxN matrix, x and b are an N-element vectors.  Note: A and b are
    //  changed upon exit to be in upper triangular (row echelon) form

    const int N = b.size();

    // A is square, with each dimension of length N

    assert(A.nrows() == N && A.ncols() == N);

    // allocate the solution array

    std::vector<double> x(N, 0.0);

    // find the scale factors for each row -- this is used when pivoting

    std::vector<double> scales(N, 0.0);

    for (int irow = 0; irow < N; ++irow) {
        double maxval = std::numeric_limits<double>::lowest();
        for (int jcol = 0; jcol < N; ++jcol) {
            maxval = std::max(maxval, std::abs(A(irow, jcol)));
        }
        scales[irow] = maxval;
    }

    // keep track of the number of times we swapped rows

    int num_row_swap{};

    if (!quiet) {
        print_Ab(A, b);
    }

    // main loop over rows

    for (int krow = 0; krow < N; ++krow) {

        // find the pivot row based on the size of column k -- only consider the
        // rows beyond the current row

        int row_max{krow};
        double col_max = std::numeric_limits<double>::lowest();
        for (int kk = krow; kk < N; ++kk) {
            double scaled_col_val = std::abs(A(kk, krow)) / scales[kk];
            if (scaled_col_val > col_max) {
                col_max = scaled_col_val;
                row_max = kk;
            }
        }

        // swap the row with the largest scaled element in the current column
        // with the current row (pivot) -- do this with b too!

        if (row_max != krow) {

            // swap row krow with row row_max

            for (int jcol = 0; jcol < N; ++jcol) {
                double tmp_val = A(krow, jcol);
                A(krow, jcol) = A(row_max, jcol);
                A(row_max, jcol) = tmp_val;
            }

            // now swap the same elements in b

            double tmp_b = b[krow];
            b[krow] = b[row_max];
            b[row_max] = tmp_b;

            // finally swap the scales

            double tmp_scale = scales[krow];
            scales[krow] = scales[row_max];
            scales[row_max] = tmp_scale;

            if (!quiet) {
                std::cout << "pivoted" << std::endl;
            }

            ++num_row_swap;
        }

        // do the forward-elimination for all rows below the current

        for (int irow = krow+1; irow < N; ++irow) {
            double coeff = A(irow, krow) / A(krow, krow);

            for (int jcol = krow+1; jcol < N; ++jcol) {
                A(irow, jcol) += -A(krow, jcol) * coeff;
            }

            A(irow, krow) = 0.0;
            b[irow] += -coeff * b[krow];
        }

        if (!quiet) {
            print_Ab(A, b);
        }
    }

    // back-substitution

    // last solution is easy

    x[N-1] = b[N-1] / A(N-1, N-1);

    for (int irow = N-2; irow >= 0; --irow) {
        double bsum = b[irow];
        for (int jcol = irow+1; jcol < N; ++jcol) {
            bsum += -A(irow, jcol) * x[jcol];
        }
        x[irow] = bsum / A(irow, irow);
    }

    return x;
}

#endif
