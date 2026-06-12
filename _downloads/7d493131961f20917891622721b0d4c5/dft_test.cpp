#include <iostream>
#include <iomanip>
#include <vector>

#include <dft.H>

int main() {

    // test our DFT implementation by taking the transform of a single-mode sine
    // wave

    double nu_0{0.2};
    double xmax{50.0};

    // number of points in our sample
    int N{128};

    std::vector<double> x(N, 0.0);
    std::vector<double> f(N, 0.0);

    double dx = xmax / static_cast<double>(N);

    // this will create evenly spaced bins from [0, xmax) -- i.e. xmax is
    // excluded.

    for (int n = 0; n < N; ++n) {
        x[n] = static_cast<double>(n) * dx;
        f[n] = std::sin(2.0 * M_PI * nu_0 * x[n]);
    }

    // take the FFT

    auto f_k = dft(f);

    // compute frequency

    std::vector<double> nu_k(N, 0.0);

    for (int n = 0; n < N; ++n) {
        nu_k[n] = static_cast<double>(n) / xmax;
    }

    // output

    for (int n = 0; n < N; ++n) {
        std::cout << std::setw(5) << n
                  << std::setw(12) << x[n]
                  << std::setw(12) << nu_k[n]
                  << std::setw(12) << f[n]
                  << std::setw(14) << std::real(f_k[n])
                  << std::setw(14) << std::imag(f_k[n]) << std::endl;
    }


}