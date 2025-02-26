#ifndef ORBIT_STATE_H
#define ORBIT_STATE_H

#include <iostream>
#include <iomanip>

struct OrbitState {
    // a container to hold the star positions
    double mass{};
    double x{};
    double y{};
    double u{};
    double v{};

    OrbitState(double M,
               double x0, double y0, double u0, double v0)
        : mass(M), x(x0), y(y0), u(u0), v(v0)
    {}

    OrbitState() {}

};

inline
std::ostream& operator<< (std::ostream& os, const OrbitState& s) {
    os.precision(6);

    os << std::setw(14) << s.mass
       << std::setw(14) << s.x
       << std::setw(14) << s.y
       << std::setw(14) << s.u
       << std::setw(14) << s.v;

    return os;
}

#endif
