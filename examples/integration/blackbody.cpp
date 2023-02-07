#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <functional>
#include <vector>

// example of integrating to infinity.  We will integrate
//
// x**3 / (exp(x) - 1)
//
// from 0 to infinity.  This arises from integrating the Planck function over
// all wavelengths.

const double SMALL = 1.e-30;

inline double zv(const double x, const double c) {
    return x / (c + x);
}

inline double xv(const double z, const double c) {
    return c * z / (1.0 - z + SMALL);
}

double integrand(double x) {
    double I{0};

    if (x == 0) {
        I = 0.0;
    } else {
        I = std::pow(x, 3) / (std::exp(x) - 1.0);
    }
    return I;
}

double I_analytic() {
    return std::pow(M_PI, 4) / 15.0;
}

// implement the composite simpsons rule.  Here N is the number
// of intervals, so there will be N+1 points from xmin to xmax
// this version is written to integrate from [0, oo]
double simpsons_inf(const int N,
                    std::function<double(const double)> func) {

    assert(N % 2 == 0);

    double I = 0.0;
    double c = 3.0;
    
    // we will work in the transformed coords z = [0, 1]
    double dz = (1.0 - 0.0) / static_cast<double>(N);

    for (int n = 0; n < N; n += 2) {
        double z0 = static_cast<double>(n) * dz;
        double z1 = static_cast<double>(n+1) * dz;
        double z2 = static_cast<double>(n+2) * dz;

        double f0 = func(xv(z0, c));
        double f1 = func(xv(z1, c));
        double f2 = func(xv(z2, c));            

        I += dz / 3.0 * (f0 / std::pow(1.0 - z0 + SMALL, 2) +
                         4.0 * f1 / std::pow(1.0 - z1 + SMALL, 2) +
                         f2 / std::pow(1.0 - z2 + SMALL, 2));
    }

    I *= c;

    return I;

}


// implement the composite trapezoid rule.  Here N is the number
// of intervals, so there will be N+1 points from xmin to xmax
// this version is written to integrate from [0, oo]
double trap_inf(const int N,
                std::function<double(const double)> func) {

    double I = 0.0;
    double c = 3.0;
    
    // we will work in the transformed coords z = [0, 1]
    double dz = (1.0 - 0.0) / static_cast<double>(N);

    for (int n = 0; n < N; ++n) {
        double z0 = static_cast<double>(n) * dz;
        double z1 = static_cast<double>(n+1) * dz;

        double f0 = func(xv(z0, c));
        double f1 = func(xv(z1, c));

        I += 0.5 * dz * (f0 / std::pow(1.0 - z0 + SMALL, 2) +
                         f1 / std::pow(1.0 - z1 + SMALL, 2));
    }

    I *= c;

    return I;

}

int main() {

    std::vector<double> bins{2, 4, 8, 16, 32, 64};  

    for (auto N : bins) {
        auto I = trap_inf(N, integrand);
        double err = std::abs(I - I_analytic());
        std::cout << std::setw(3) << N << ": " << I << " " << err << std::endl;
    }
}