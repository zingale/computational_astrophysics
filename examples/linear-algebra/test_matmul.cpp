#include <vector>
#include <array.H>
#include <iostream>
#include <matmul.H>

int main() {

    Array A{{{1, 2, 3},
             {4, 5, 6}}};

    std::cout << A << std::endl;

    std::vector<double> x{1, 2, 3};

    auto b = matmul(A, x);

    for (auto element : b) {
        std::cout << element << std::endl;
    }
}
