#include <iostream>
#include <iomanip>
#include <cmath>
#include <print>

// compute the triple alpha reaction rate / rho^2 Y^3
double q3a(const double T) {
    double T8{T/1.e8};
    return 5.09e11 * std::pow(T8, -3.0) * std::exp(-44.027/T8);
}

int main() {

    for (auto T0 : {1.e8, 2.5e8, 5.e8, 1.e9, 2.5e9, 5e9}) {
        double dT{1.e-8 * T0};
        double nu = (T0 / q3a(T0)) * (q3a(T0 + dT) - q3a(T0)) / dT;

        std::println("{:8.5g} : {:5.2f}", T0, nu);
    }
}
