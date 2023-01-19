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
}