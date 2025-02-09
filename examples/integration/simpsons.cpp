#include <iostream>
#include <iomanip>
#include <cmath>
#include <functional>
#include <cassert>
#include <numbers>

// An implementation of the composite Simpsons rule for integration.  This
// assumes that we will use an even number of intervals

// note: requires C++20 due to std::numbers

// the function we will integrate
double f(const double x) {
    return 1.0 + x * std::sin(std::numbers::pi * x) / 4.0;
}

// the analytic integral of f(x) over [1/2, 3/2]
double I_analytic() {
    return 1.0 - 0.5 / std::pow(std::numbers::pi, 2);
}

// implement the composite simpsons rule.  Here N is the number
// of intervals, so there will be N+1 points from xmin to xmax
double simpsons(const double xmin, const double xmax, const int N,
                std::function<double(const double)> func) {

    assert(N % 2 == 0);

    double I = 0.0;
    double dx = (xmax - xmin) / static_cast<double>(N);

    for (int n = 0; n < N; n += 2) {
        double x0 = xmin + static_cast<double>(n) * dx;
        I += dx / 3.0 * (func(x0) + 4.0 * func(x0 + dx) + func(x0 + 2.0*dx));
    }

    return I;
}

int main() {

    double xmin{0.5};
    double xmax{1.5};

    for (auto N : {2, 4, 8, 16, 32, 64}) {
        double I = simpsons(xmin, xmax, N, f);
        double err = std::abs(I - I_analytic());
        std::cout << std::setw(3) << N << ": " << err << std::endl;
    }
}
