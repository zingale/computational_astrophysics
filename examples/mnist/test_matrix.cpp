#include <iostream>
#include <vector>

#include "matrix.H"


int main() {

    Matrix A{{1, 2, 3},
             {4, 5, 6},
             {7, 8, 9},
             {10, 11, 12}};

    std::vector<double> x{0, 1, 2};

    auto b = A * x;
    std::cout << b << std::endl;

    std::cout << A.transpose() * b << std::endl;

    Matrix B{{1, 2}, {3, 4}};

    B.apply([](double x) -> double {return 2*x+1;});
    std::cout << B << std::endl;

    Matrix C{{1, 2, 3}, {4, 5, 6}};
    std::cout << C * 2 << std::endl;
    std::cout << 3 * C << std::endl;

    Matrix P{{1, 2, 3}, {4, 5, 6}};
    Matrix Q{{2, 4, 6}, {8, 10, 12}};

    std::cout << P % Q << std::endl;

    std::cout << P + Q << std::endl;
    std::cout << P - Q << std::endl;
}
