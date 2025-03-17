#include <algorithm>
#include <iostream>
#include <random>

#include "matrix.H"

int main() {

    Matrix A(4, 5);

    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> randn(0.0, 1.0 / std::sqrt(4));
    std::generate(A._data.begin(), A._data.end(),
                  [&]() -> double {return randn(generator);});

    std::cout << A << std::endl;
}
