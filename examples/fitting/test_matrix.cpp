#include <iostream>
#include <vector>

#include "matrix.H"


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
}
