#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <functional>
#include <limits>

// bisecton
double bisect(const double _xl, const double _xr,
              std::function<double(double)> f,
              const double tol=1.e-8) {

    double xl{_xl};
    double xr{_xr};

    // first check if there is a root in the interval
    double fl = f(xl);
    double fr = f(xr);

    if (fl * fr > 0.0) {
        std::cout << "no root in the interval" << std::endl;
        return std::numeric_limits<double>::infinity();
    }

    // main iteration loop

    double err = std::numeric_limits<double>::max();
    double xm = 0.5 * (xl + xr);

    while (err > tol) {
        double fm = f(xm);

        if (fl * fm > 0.0) {
            // the root is in the right half of the interval
            xl = xm;
            fl = fm;
        } else {
            // the root is in the right half of the interval
            xr = xm;
            fr = fm;
        }

        xm = 0.5 * (xl + xr);
        err = std::abs(xr - xl) / std::abs(xm);

    }

    return xm;
}

int main() {

    // let's try this out on a cubic
    auto r = bisect(-2.0, -1.0,
                    [](double x){return 0.5 * std::pow(x, 3) + M_PI/3.0 * x + 2;});
    std::cout << "the root is " << r << std::endl;
}