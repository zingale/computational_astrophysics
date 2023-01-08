#include <iostream>
#include <limits>

int main() {

    std::cout << "maximum double = " << std::numeric_limits<double>::max() << std::endl;
    std::cout << "maximum double base-10 exponent = " << std::numeric_limits<double>::max_exponent10 << std::endl;
    std::cout << "smallest (abs) double = " << std::numeric_limits<double>::min() << std::endl;
    std::cout << "minimim double base-10 exponent = " << std::numeric_limits<double>::min_exponent10 << std::endl;
    std::cout << "machine epsilon (double) = " << std::numeric_limits<double>::epsilon() << std::endl;

}
