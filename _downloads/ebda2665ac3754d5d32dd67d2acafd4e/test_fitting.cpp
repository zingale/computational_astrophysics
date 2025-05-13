#include <iostream>
#include <fstream>
#include <random>
#include <vector>

#include "fitting.H"


std::pair<std::vector<double>, std::vector<double>>
y_experiment(const std::vector<double>& av, double sigma,
             const std::vector<double>& x) {

    auto N = x.size();

    std::vector<double> y(N, 0.0);

    // compute Gaussian random numbers that will represent
    // the error for each point
    std::mt19937 generator(12345);
    std::normal_distribution<double> randn(0.0, sigma);

    std::vector<double> yerr(N, 0.0);

    for (std::size_t n = 0; n < N; ++n) {
        auto r = randn(generator);
        y[n] = av[0] + av[1] * x[n] + av[2] * x[n] * x[n] + r;
        yerr[n] = sigma;
    }

    return {y, yerr};
}

int main() {

    constexpr std::size_t N{40};
    constexpr double xmin{0.0};
    constexpr double xmax{100.0};

    std::vector<double> x(N, 0.0);
    double dx = (xmax - xmin) / static_cast<double>(N - 1);

    for (std::size_t n = 0; n < N; ++n) {
        x[n] = xmin + static_cast<double>(n) * dx;
    }

    double sigma{8.0};
    auto [y, yerr] = y_experiment({2.0, 1.50, -0.02}, sigma, x);

    std::ofstream of("experiment.dat");
    for (std::size_t n = 0; n < N; ++n) {
        of << x[n] << " " << y[n] << " " << yerr[n] << std::endl;
    }

    constexpr int M{3};
    auto [a, chisq] = general_regression(x, y, yerr, M);

    std::cout << "chisq = " << chisq << std::endl;

    std::cout << "a = ";
    for (auto e : a) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

}
