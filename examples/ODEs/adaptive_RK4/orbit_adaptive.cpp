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
    double x;
    double y;
    double u;
    double v;

    OrbitState(double x0, double y0, double u0, double v0)
        : x{x0}, y{y0}, u{u0}, v{v0}
    {}

    OrbitState()
        : x(0.0), y(0.0), u(0.0), v(0.0)
    {}

    OrbitState operator+(const OrbitState& other) const {
        return OrbitState(x + other.x, y + other.y, u + other.u, v + other.v);
    }

    OrbitState operator-(const OrbitState& other) const {
        return OrbitState(x - other.x, y - other.y, u - other.u, v - other.v);
    }

    OrbitState operator*(double a) const {
        return OrbitState(a * x, a * y, a * u, a * v);
    }

    // a * vec -- needs to be a friend
    friend OrbitState operator*(double a, const OrbitState& state);

};

std::ostream& operator<< (std::ostream& os, const OrbitState& s) {
    os.precision(6);

    os << std::setw(14) << s.x
       << std::setw(14) << s.y
       << std::setw(14) << s.u
       << std::setw(14) << s.v;

    return os;
}

OrbitState operator*(double a, const OrbitState& state) {
    return OrbitState(a * state.x, a * state.y, a * state.u, a * state.v);
}

class OrbitsRK4 {
    // model the evolution of a single planet around the Sun using gravitational
    // interaction of three stars

private:

    // a vector to store the history of our orbit
    std::vector<OrbitState> history;
    std::vector<double> time;
    int n_reset{0};

public:

    OrbitsRK4(const double a, const double e)
        {

            OrbitState initial_conditions;

            // put the planet at perihelion

            initial_conditions.x = 0.0;
            initial_conditions.y = a * (1.0 - e);
            initial_conditions.u = std::sqrt((GM / a) * (1.0 + e) / (1.0 - e));
            initial_conditions.v = 0.0;

            history.push_back(initial_conditions);

            time.push_back(0.0);

        }

    int npts() {
        // return the number of integration points
        return time.size();
    }

    int get_n_reset() {
        // return the number of times a step was reset
        return n_reset;
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

        double KE = 0.5 * (state.u * state.u + state.v * state.v);
        double PE = -GM / std::sqrt(state.x * state.x + state.y * state.y);

        return KE + PE;
    }

    OrbitState rhs(const OrbitState& state) {
        // compute the ydot terms

        OrbitState ydot;

        ydot.x = state.u;
        ydot.y = state.v;

        double dx = state.x;
        double dy = state.y;

        double r = std::sqrt(dx * dx + dy * dy);

        ydot.u = -GM * dx / std::pow(r, 3);
        ydot.v = -GM * dy / std::pow(r, 3);

        return ydot;

    }

    OrbitState single_step(const OrbitState& state_old, const double dt) {
        /// take a single RK-4 timestep through dt

        auto ydot1 = rhs(state_old);

        OrbitState state_temp;

        state_temp = state_old + 0.5 * dt * ydot1;
        auto ydot2 = rhs(state_temp);

        state_temp = state_old + 0.5 * dt * ydot2;
        auto ydot3 = rhs(state_temp);

        state_temp = state_old + dt * ydot3;
        auto ydot4 = rhs(state_temp);

        OrbitState s = state_old + (dt / 6.0) *
            (ydot1 + 2.0 * ydot2 + 2.0 * ydot3 + ydot4);

        return s;
    }


    void integrate(const double dt_in, const double err, const double tmax) {
        // integrate the equations of motion using 4th order R-K
        // method with an adaptive stepsize, to try to achieve the
        // relative error err.  dt here is the initial timestep

        // if err < 0, then we don't do adaptive stepping, but rather
        // we always walk at the input dt

        // start with the old timestep
        double dt_new{dt_in};
        double dt{dt_in};

        double t{0.0};

        while (t < tmax) {

            const auto &state_old = history.back();

            // adaptive stepping iteration loop -- keep trying to take
            // a step until we achieve our desired error
            double rel_error = std::numeric_limits<double>::max();

            int n_try{0};

            OrbitState state_new;

            while (rel_error > err) {

                dt = std::min(dt_new, tmax - t);

                // take 2 half steps

                auto state_temp = single_step(state_old, 0.5*dt);
                state_new = single_step(state_temp, 0.5*dt);

                // now take just a single step to cover dt

                auto state_single = single_step(state_old, dt);

                // state_new should be more accurate than state_single
                // since it used smaller steps

                // estimate the relative error now

                double pos_err =
                    std::max(std::abs((state_new.x - state_single.x) /
                                      state_single.x),
                             std::abs((state_new.y - state_single.y) /
                                      state_single.y));

                double vel_err =
                    std::max(std::abs((state_new.u - state_single.u) /
                                      state_single.u),
                             std::abs((state_new.v - state_single.v) /
                                      state_single.v));

                rel_error = std::max(pos_err, vel_err);

                // adaptive timestep algorithm from Garcia (Eqs. 3.30
                // and 3.31)

                double dt_est = dt * std::pow(std::abs(err/rel_error), 0.2);
                dt_new = std::min(std::max(0.9 * dt_est, 0.25 * dt), 4.0 * dt);

                n_try++;
            }

            if (n_try > 1) {
                // n_try = 1 if we took only a single try at the step
                n_reset += (n_try - 1);
            }

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

    OrbitsRK4 o(1.0, 0.95);
    o.integrate(0.05, 1.e-5, 1.0);

    std::cout << "number of points = " << o.npts() << std::endl;
    std::cout << "number of times a step was rejected = "
              << o.get_n_reset() << std::endl;

    std::ofstream of("orbit_adaptive.dat");

    for (int n = 0; n < o.npts(); ++n) {
        const auto& state = o.get_state(n);

        of << std::setw(14) << o.get_time(n);

        of << state;

        of.precision(10);

        of << std::setw(18) << o.energy(n);

        of << std::endl;

    }

}
