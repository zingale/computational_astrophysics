#include <iostream>

#include "vector2d.H"

int main() {

    // create 2 vectors

    Vector2d v1(1, 2);
    Vector2d v2(2, 4);

    // output our vectors

    std::cout << v1 << " " << v2 << std::endl;

    // output their sum

    std::cout << v1 + v2 << std::endl;

    // create a new vector from subracting our two vectors

    auto v3 = v1 - v2;
    std::cout << v3 << std::endl;

    // create a copy

    auto v4 = v3;
    std::cout << v3 << " " << v4 << std::endl;

    // change the data in the original
    v3.set_x(0.0);
    v3.set_y(0.0);

    // did both change? or just the original?

    std::cout << v3 << " " << v4 << std::endl;


    // vectors displancing N, E, S, W

    Vector2d v_north(0.0, 1.0);
    Vector2d v_east(1.0, 0.0);
    Vector2d v_south(0.0, -1.0);
    Vector2d v_west(-1.0, 0.0);

    // test compound addition -- if we sum up the four vectors above,
    // we should get zero

    Vector2d vv1;
    vv1 += v_north;
    vv1 += v_east;

    std::cout << "compound addition (intermediate result): " << vv1 << std::endl;
    vv1 += v_south;
    vv1 += v_west;

    std::cout << "compound addition (final result): " << vv1 << std::endl;

    // test compound subtraction

    Vector2d vv2(1.0, 1.0);
    vv2 -= v_north;

    std::cout << "compound subtraction: " << vv2 << std::endl;

    // test compound multiplication

    Vector2d vv3(2.0, 4.0);
    vv3 *= 4;

    std::cout << "compound multiplication: " << vv3 << std::endl;

    // test compound division

    Vector2d vv4(1.0, -5.0);
    vv4 /= 2;

    std::cout << "compound division: " << vv4 << std::endl;


}
