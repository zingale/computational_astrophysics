#include <iostream>

#include "matrix.H"

int main() {

    Matrix A{{{1, 2, 3},
              {4, 5, 6},
              {7, 8, 9}}};

    std::cout << "matrix A =\n" << A << std::endl;

    auto B = A * 2;

    auto C = 2 * A;


    std::cout << "B = A * 2 =\n" << B << std::endl;

    std::cout << "C = 2 * A =\n" << C << std::endl;

    std::cout << "A + B =\n" << A + B << std::endl;

    Matrix K{{{1}, {2}, {3}, {4}}};
    Matrix M{{{1, 2, 3, 4}}};

    std::cout << "matrix K =\n" << K << std::endl;
    std::cout << "matrix M =\n" << M << std::endl;

    std::cout << "K * M =\n" << K * M << std::endl;
    std::cout << "M * K =\n" << M * K << std::endl;

    // square each element
    std::cout << "squaring the elements of A:\n" << A.transform([] (double x) {return x*x;}) << std::endl;

}
