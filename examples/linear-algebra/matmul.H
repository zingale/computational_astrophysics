#ifndef MATMUL_H
#define MATMUL_H

#include <vector>
#include <array.H>

inline
std::vector<double> matmul(const Array& A, const std::vector<double>& x) {

    assert (A.ncols() == static_cast<int>(x.size()));

    std::vector<double> b(A.nrows(), 0.0);

    for (int m = 0; m < A.nrows(); ++m) {
        for (int n = 0; n < A.ncols(); ++n) {
            b[m] += A(m, n) * x[n];
        }
    }

    return b;

}

#endif
