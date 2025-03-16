#include <iostream>
#include <vector>

#include "matrix.H"
#include "gauss.H"

int main() {

    Matrix A(3, 4);
    int ival = 0;
    for (std::size_t irow = 0; irow < A.nrows(); ++irow) {
        for (std::size_t jcol = 0; jcol < A.ncols(); ++jcol) {
            A(irow, jcol) = ival++;
        }
    }

    std::cout << A << std::endl;
    std::cout << A.transpose() << std::endl;

    auto AT = A.transpose();
    auto C = AT * A;
    std::cout << "A^T A =\n " << C << std::endl;

    std::vector<double> x{0, 1, 2, 3};

    auto b = A * x;
    for (auto e : b) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

    // try solvig A x = b
    Matrix A_new{{1, 4, -2},
                 {9, 12, 3},
                 {-1, 3, 8}};

    std::cout << "A_new = \n" << A_new <<std::endl;

    std::vector<double> x3{0, 1, 2};
    auto b_new = A_new * x3;

    std::cout << "b_new = ";
    for (auto e : b_new) {
        std::cout << e << " ";
    }
    std::cout << std::endl;

    auto x_new = gauss_elim(A_new, b_new);
    for (auto e : x_new) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
}
