#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "matrix.H"
#include "gauss.H"

std::pair<std::vector<double>, double>
general_regression(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& yerr,
                   int M) {

    auto N = x.size();
    assert (y.size() == N && yerr.size() == N);

    // make the design matrix
    // each row corresponds to a single data point
    // each column corresponds to a particular basis function

    Matrix A(N, M);
    for (int n = 0; n < N; ++n) {
        A(n, 0) = 1.0 / yerr[n];
        for (int m = 1; m < M; ++m) {
            // x**m / yerr
            A(n, m) = A(n, m-1) * x[n];
        }
    }

    // build the linear system
    auto ATA = A.transpose() * A;

    std::vector<double> source(N, 0);
    for (int n = 0; n < N; ++n) {
        source[n] = y[n] / yerr[n];
    }

    auto b = A.transpose() * source;

    auto a = gauss_elim(ATA, b);

    double chisq{};
    for (int n = 0; n < N; ++n) {
        double y_model{};
        for (int m = 0; m < M; ++m) {
            y_model += a[m] * std::pow(x[n], m);
        }
        chisq += std::pow(y_model - y[n], 2) / std::pow(yerr[n], 2);
    }
    chisq /= static_cast<double>(N - M);

    return {a, chisq};
}
