#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <functional>
#include <fstream>

class LagrangePoly {
    // create a lagrange polynomial that fits through N
    // points in [xmin, xmax] sampling the function func

public:

    // our data samples that we will interpolate through

    std::vector<double> xp;
    std::vector<double> fp;

    // constructor -- this takes a function to create the
    // data samples

    LagrangePoly(double xmin, double xmax, int N,
                 std::function<double(double)> func)
        : xp(N, 0.0), fp(N, 0.0)
    {
        double dx = (xmax - xmin) / static_cast<double>(N-1);

        for (int i = 0; i < N; ++i) {
            xp[i] = xmin + static_cast<double>(i) * dx;
            fp[i] = func(xp[i]);
        }
    }

    // construct and evaluate the Lagrange polynomial at x0

    double evalf(double x0) {

        double f{};

        for (std::size_t m = 0; m < xp.size(); ++m) {

            // create the Lagrange basis function for point m
            double l{1};

            for (std::size_t n = 0; n < xp.size(); ++n) {
                if (n == m) {
                    continue;
                }
                l *= (x0 - xp[n]) / (xp[m] - xp[n]);
            }
            f += fp[m] * l;
        }

        return f;
    }

};


// a test function

double f_test(double x) {
    // a simple tanh profile
    return 0.5 * (1.0 + std::tanh((x - 1.0) / 0.1));
}


// now a test driver

int main() {

    // setup the interpolant and store the points in a file for plotting

    int N{11};
    double xmin{0.0};
    double xmax{2.0};

    LagrangePoly l(xmin, xmax, N, f_test);

    std::ofstream of("data_samples.txt");
    for (std::size_t i = 0; i < l.xp.size(); ++i) {
        of << std::setw(14) << l.xp[i] << " " << std::setw(14) << l.fp[i] << std::endl;
    }

    // now finely sample the domain and interpolate our samples and
    // output the interpolation to another file

    std::ofstream of2("interpolation.txt");

    int N_fine{200};
    double dx_fine = (xmax - xmin) / static_cast<double>(N_fine-1);

    for (int i = 0; i < N_fine; ++i) {
        double x = xmin + static_cast<double>(i) * dx_fine;

        auto f_interp = l.evalf(x);
        auto f_true = f_test(x);

        of2 << std::setw(14) << x << " "
            << std::setw(14) << f_interp << " "
            << std::setw(14) << f_true << std::endl;

    }
}
