#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <limits>
#include <fstream>

const double GM = 4.0 * M_PI * M_PI;

struct OrbitState {
    // a container to hold the star positions
    double x{};
    double y{};
    double u{};
    double v{};

    OrbitState(double x0, double y0, double u0, double v0)
        : x{x0}, y{y0}, u{u0}, v{v0}
    {}

    OrbitState() {}

    OrbitState operator+(const OrbitState& other) const {
        return OrbitState(x + other.x, y + other.y, u + other.u, v + other.v);
    }

    OrbitState operator-(const OrbitState& other) const {
        return OrbitState(x - other.x, y - other.y, u - other.u, v - other.v);
    }

    // this handles OrbitState * a
    OrbitState operator*(double a) const {
        return OrbitState(a * x, a * y, a * u, a * v);
    }

};

inline
std::ostream& operator<< (std::ostream& os, const OrbitState& s) {
    os.precision(6);

    os << std::setw(14) << s.x
       << std::setw(14) << s.y
       << std::setw(14) << s.u
       << std::setw(14) << s.v;

    return os;
}

// this handles a * OrbitState
inline
OrbitState operator*(double a, const OrbitState& state) {
    return OrbitState(a * state.x, a * state.y, a * state.u, a * state.v);
}


class OrbitsRK4 {
    // model the evolution of a single planet around the Sun using
    // gravitational interaction of three stars

private:

    // a vector to store the history of our orbit
    std::vector<OrbitState> history;
    std::vector<double> time;

public:

    OrbitsRK4(const double a, const double e)
    {

        // put the planet at perihelion

        double x0 = 0.0;
        double y0 = a * (1.0 - e);
        double u0 = std::sqrt((GM / a) * (1.0 + e) / (1.0 - e));
        double v0 = 0.0;

        history.push_back(OrbitState(x0, y0, u0, v0));

        time.push_back(0.0);

    }

    int npts() {
        // return the number of integration points
        return time.size();
    }

    double get_time(const int n) {
        // return the physical time for step n
        return time[n];
    }

    OrbitState& get_state(const int n) {
        // return a reference to the state at time index n
        return history[n];
    }

    double energy(const int n) {
        // return the energy of the system for timestep n

        const auto& state = history[n];

        // kinetic energy

        double KE = 0.5 * (std::pow(state.u, 2) + std::pow(state.v, 2));
        double PE = -GM / std::sqrt(std::pow(state.x, 2) + std::pow(state.y, 2));

        return KE + PE;
    }

    OrbitState rhs(const OrbitState& state) {
        // compute the ydot terms

        double dxdt = state.u;
        double dydt = state.v;

        double dx = state.x;
        double dy = state.y;

        double r = std::sqrt(dx * dx + dy * dy);

        double dudt = -GM * dx / std::pow(r, 3);
        double dvdt = -GM * dy / std::pow(r, 3);

        return OrbitState(dxdt, dydt, dudt, dvdt);

    }

    void integrate(const double dt_in, const double tmax) {
        // integrate the equations of motion using 4th order R-K
        // method

        double t{0.0};
        double dt{dt_in};

        while (t < tmax) {

            const auto &state_old = history.back();

            auto ydot1 = rhs(state_old);

            OrbitState state_temp{};

            state_temp = state_old + 0.5 * dt * ydot1;
            auto ydot2 = rhs(state_temp);

            state_temp = state_old + 0.5 * dt * ydot2;
            auto ydot3 = rhs(state_temp);

            state_temp = state_old + dt * ydot3;
            auto ydot4 = rhs(state_temp);

            // final solution

            OrbitState state_new = state_old + (dt / 6.0) *
                (ydot1 + 2.0 * ydot2 + 2.0 * ydot3 + ydot4);

            // successful step
            t += dt;

            // store the new solution -- this will automatically be
            // used as the "old" state in the next step

            time.push_back(t);
            history.push_back(state_new);

        }

    }
};


int main() {

    OrbitsRK4 o(1.0, 0.3);
    o.integrate(0.05, 2.0);

    std::ofstream of("orbit_rk4.dat");

    for (int n = 0; n < o.npts(); ++n) {
        const auto& state = o.get_state(n);

        of << std::setw(14) << o.get_time(n);

        of << state;

        of.precision(10);

        of << std::setw(18) << o.energy(n);

        of << std::endl;

    }

}
