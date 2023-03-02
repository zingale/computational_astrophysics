#include <iostream>
#include <limits>
#include <vector>

#include <array.H>
#include <gauss.H>

constexpr double sigma{10.0};
constexpr double rho{28.0};
constexpr double beta{8.0/3.0};

constexpr double tol{1.e-8};

std::vector<double> rhs(const std::vector<double>& x) {

    std::vector<double> f(x.size(), 0.0);

    f[0] = sigma * (x[1] - x[0]);
    f[1] = rho * x[0] - x[1] - x[0] * x[2];
    f[2] = x[0] * x[1] - beta * x[2];

    return f;

}

Array jacobian(const std::vector<double>& x) {

    Array J(x.size(), x.size());

    J(0, 0) = -sigma;
    J(0, 1) = sigma;
    J(0, 2) = 0.0;

    J(1, 0) = rho - x[2];
    J(1, 1) = -1.0;
    J(1, 2) = -x[0];

    J(2, 0) = x[1];
    J(2, 1) = x[0];
    J(2, 2) = -beta;

    return J;

}


std::vector<double> lorenz_roots(const std::vector<double>& x0) {

    std::vector<double> x(x0);

    double err = std::numeric_limits<double>::max();

    while (err > tol) {

        auto J = jacobian(x);
        auto f = rhs(x);

        for (auto& e : f) {
            e *= -1;
        }

        // solve J dx = -f
        auto dx = gauss_elim(J, f, true);

        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            x[i] += dx[i];
        }

        // we'll use the max norm
        err = std::numeric_limits<double>::lowest();
        for (auto e : dx) {
            err = std::max(err, std::abs(e));
        }
    }

    return x;
}

int main() {

    // pick an initial guess

    std::vector<double> x0{10.0, 10.0, 10.0};

    auto x = lorenz_roots(x0);

    for (auto e : x) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

}