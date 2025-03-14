#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <print>

// our test function

double f(double x) {
    return x * std::sin(x);
}


// a function that will produce N samples of the input func(x)

std::pair<std::vector<double>, std::vector<double>>
sample_function(std::function<double(double)> func,
                int npts, double xmin, double xmax) {

    std::vector<double> xv(npts, 0.0);
    std::vector<double> fv(npts, 0.0);

    double dx = (xmax - xmin) / static_cast<double>(npts - 1);

    for (int i = 0; i < npts; ++i) {
        xv[i] = xmin + static_cast<double>(i) * dx;
        fv[i] = func(xv[i]);
    }

    return {xv, fv};
}


// the main linear intepolation function

double interpolate(double x0,
                   const std::vector<double>& xv,
                   const std::vector<double>& fv) {

    // find the index i where xv[i] > x0
    int idx = -1;
    for (std::size_t i = 0; i < xv.size()-1; ++i) {
        if (x0 > xv[i] && x0 < xv[i+1]) {
            idx = i;
            break;
        }
    }

    if (idx == -1) {
        std::cout << "x0 is not in range of xv" << std::endl;
    }

    double slope = (fv[idx+1] - fv[idx]) / (xv[idx+1] - xv[idx]);
    return slope * (x0 - xv[idx]) + fv[idx];
}


// test driver

int main() {

    // we want to find the N such that the error in our interpolation
    // is < 1%.  We'll just consider a bunch of different sizes for N

    double xmin{0.0};
    double xmax{10.0};
    double x0{4.5};

    for (auto N : {5, 10, 20, 40, 80, 160}) {
        const auto [xv, fv] = sample_function(f, N, xmin, xmax);

        auto finterp = interpolate(x0, xv, fv);
        auto fexact = f(x0);

        double err = std::abs(finterp - fexact);

        std::println("{:4} points, error = {:7.5f}", N, err);
    }
}
