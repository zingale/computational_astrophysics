#include <iostream>

#include "matrix.H"
#include "inverse.H"


int main() {

    Matrix A{{1.0, -2.0, 4.0, 3.5},
             {2.4, 10.2, -8.5, -12.0},
             {1.0, 0.0, -14.5, 2.1},
             {4.7, 10.3, 11.2, 50.0}};

    auto Ainv = inverse(A);

    std::cout << Ainv << std::endl;

    std::cout << A * Ainv << std::endl;
}
