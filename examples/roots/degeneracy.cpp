// Given an electron number density, n_e, integrate the Fermi-Dirac
// distribution and find the degeneracy parameter that gives that
// number density.
//
// We assume that we are non-relativisitic and work in terms
// of dimensionless momentum, x = p / (m_e c).
//
// We are ultimately integrating
//
// n_e = (8π/h^3) (2m_e k T)^(3/2) ∫ x^2 dx / (exp(x^2 - psi) + 1)
//
// where psi is the degeneracy parameter.

#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>

namespace {
    constexpr double SMALL{1.e-30};
    constexpr double h_planck{6.63e-27};
    constexpr double k_B{1.38e-16};
    constexpr double m_e{9.11e-28};
    constexpr double m_u {1.67e-24};
}

inline
double integrand(const double x, const double psi) {
    return x * x / (std::exp(x*x - psi) + 1.0);
}

inline
double xv(const double z, const double c) {
    return c * z / (1.0 - z + SMALL);
}

inline
double fd_integral(const int N, const double psi) {

    assert(N % 2 == 0);

    // we will transform from integrating over x to integrating over
    // z = x / (c + x) with z = [0, 1]

    double c{2.0};

    // N is the number of intervals
    double dz = 1.0 / static_cast<double>(N);

    double I{0.0};
    for (int n = 0; n < N; n += 2) {
        double zl{n * dz};
        double zc{(n + 1) * dz};
        double zr{(n + 2) * dz};

        double fl = integrand(xv(zl, c), psi) / std::pow(1.0 - zl + SMALL, 2);
        double fc = integrand(xv(zc, c), psi) / std::pow(1.0 - zc + SMALL, 2);
        double fr = integrand(xv(zr, c), psi) / std::pow(1.0 - zr + SMALL, 2);

        I += (1.0/3.0) * dz * (fl + 4*fc + fr);
    }
    I *= c;
    return I;
}

inline
double compute_ne(const double T, const double psi,
                  double tol=1.e-8) {
    // given a temperature and degeneracy parameter, psi, compute the
    // number density of electrons by integrating the Fermi-Dirac
    // distribution over all momenta

    // we'll do the dimensionless integral of
    // x**2 / (exp(x**2 - psi) + 1) using Simpson's rule
    // we will pick a value of N and do the integration
    // and then change N until the error is small

    // we'll start with N = 50 and then keep doubling
    // until the difference in I is small

    int N{50};
    double err = std::numeric_limits<double>::max();
    double I_old = std::numeric_limits<double>::max();

    double I{};

    while (err > tol) {
        I = fd_integral(N, psi);
        err = std::abs(I - I_old) / std::abs(I_old);
        I_old = I;
        N *= 2;
    }

    double ne = 8.0 * M_PI / std::pow(h_planck, 3) *
        std::pow(2.0 * m_e * k_B * T, 1.5) * I;

    return ne;
}

inline
double find_psi(const double rho, const double T, const double Ye,
                double tol=1.e-6) {
    // use a secant method to find the psi that gives the correct
    // n_e / density (rho)

    double ne_input = rho * Ye / m_u;

    // initial guesses
    double psi_m1{10.0};
    double psi_0{20.0};

    // we want to zero ne(T, psi) - ne_input

    auto ne_m1 = compute_ne(T, psi_m1) - ne_input;

    double err = std::numeric_limits<double>::max();

    while (err > tol) {
        auto ne_0 = compute_ne(T, psi_0) - ne_input;
        double dne_dpsi = (ne_0 - ne_m1) / (psi_0 - psi_m1);

        psi_m1 = psi_0;
        ne_m1 = ne_0;

        psi_0 -= ne_0 / dne_dpsi;

        err = std::abs(ne_0 / dne_dpsi) / std::abs(psi_0);
    }

    return psi_0;
}

int main() {

    double T{1.e7};
    double Ye = 0.5;

    double rho_min{1.e-4};
    double rho_max{1.e6};

    int N{50};

    double dlogrho = (std::log10(rho_max) - std::log10(rho_min)) /
        static_cast<double>(N - 1);

    for (int nrho = 0; nrho < N; ++nrho) {
        double rho = std::pow(10.0, std::log10(rho_min) + nrho*dlogrho);

        auto psi = find_psi(rho, T, Ye);
        std::cout << rho << " " << psi << std::endl;
    }
}
