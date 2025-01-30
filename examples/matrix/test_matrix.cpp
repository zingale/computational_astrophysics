#include <iostream>

#include "matrix.H"

int main() {

    Matrix A{{{1, 2, 3},
              {4, 5, 6},
              {7, 8, 9}}};

    auto B = A * 2;

    std::cout << A << std::endl;

    std::cout << B << std::endl;

    std::cout << A + B << std::endl;
}
