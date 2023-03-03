#ifndef DFT_H
#define DFT_H

#include <complex>
#include <vector>

std::vector<std::complex<double>> dft(std::vector<double> f) {
    // compute the discrete Fourier transform of the real-valued function f

    using namespace std::complex_literals;

    std::vector<std::complex<double>> f_k(f.size(), 0.0);

    int N = static_cast<int>(f.size());

     for (int k = 0; k < N; ++k) {
        for (int n = 0; n < N; ++ n) {
            f_k[k] += f[n] * std::exp(-2.0 * M_PI * 1i * static_cast<double>(n * k) / static_cast<double>(N));
        }
     }

    return f_k;
}

#endif
