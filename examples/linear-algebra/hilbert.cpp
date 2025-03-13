#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "array.H"
#include "gauss.H"
#include "matmul.H"

int main() {

    // loop over N
    for (int n = 2; n < 100; ++n) {

        // create Hilbert matrix of size NxN
        Array H(n, n);

        for (int irow = 0; irow < n; ++irow) {
            for (int jcol = 0; jcol < n; ++jcol) {
                H(irow, jcol) = 1.0 / static_cast<double>(irow + jcol + 1);
            }
        }

        // create vector x
        std::vector x(n, 0.0);
        for (int i = 0; i < n; ++i) {
            x[i] = static_cast<double>(i);
        }

        // create b = H x
        auto b = matmul(H, x);

        // solve H xtilde = b
        auto x_tilde = gauss_elim(H, b, true);

        // compute error | x - xtilde |

        double error{-1.0};
        for (int i = 0; i < n; ++i) {
            error = std::max(error, std::abs(x[i] - x_tilde[i]));
        }

        //std::cout << std::setprecision(10);

        //for (int i = 0; i < n; ++i) {
        //    std::cout << i << " " << x[i] << " " << x_tilde[i] << std::endl;
        //}

        // output
        std::cout << "for n = " << n << " " << "error = " << error << std::endl;

        // if error ~ 1 then exit
        if (error > 1.0) {
            break;
        }
    }

}
