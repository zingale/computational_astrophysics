#include <iostream>

#include <array.H>
#include <gauss.H>

int main() {

    Array A{{{1, 1, 1},
             {-1, 2, 0},
             {2, 0, 1}}};

    std::vector<double> b{6, 3, 5};

    auto x = gauss_elim(A, b);

    for (auto e : x) {
        std::cout << e << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;

    // next example

    Array A2{{{0, 0, 0, 4},
              {0, 0, 3, 0},
              {5, 6, 7, 8},
              {0, 4, 3, 2}}};

    std::vector<double> b2{5, 4, 9, 1};

    auto x2 = gauss_elim(A2, b2);

    for (auto e : x2) {
        std::cout << e << std::endl;
    }

    std::cout << std::endl;
    std::cout << std::endl;


}