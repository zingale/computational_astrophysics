#include <algorithm>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <limits>
#include <fstream>
#include <utility>

#include "orbit.H"

const double GM = 4.0 * M_PI * M_PI;


class OrbitsRK4 {
    // model the evolution of a single planet around the Sun using gravitational
    // interaction of three stars

private:

    // a vector to store the history of our orbit
    std::vector<OrbitState> history;
    std::vector<double> time;

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

        OrbitState ydot;

        ydot.x = state.u;
        ydot.y = state.v;

        double dx = -state.x;
        double dy = -state.y;

        double r = std::sqrt(dx * dx + dy * dy);

        ydot.u = GM * dx / std::pow(r, 3);
        ydot.v = GM * dy / std::pow(r, 3);

        return ydot;

    }

    std::pair<OrbitState, OrbitState>
    single_step(OrbitState& state_old, const double dt) {
        /// take a single RK-4 timestep through dt

        auto f1 = dt * rhs(state_old);

        OrbitState state_temp;

        state_temp = state_old + 0.25 * f1;
        auto f2 = dt * rhs(state_temp);

        state_temp = state_old + (3./32.) * f1 + (9./32.) * f2;
        auto f3 = dt * rhs(state_temp);

        state_temp = state_old + (1932./2197.) * f1 - (7200./2197.) * f2  + (7296./2197.) * f3;
        auto f4 = dt * rhs(state_temp);

        state_temp = state_old + (439./216.) * f1 - 8. * f2 + (3680./513.) * f3 - (845./4104.) * f4;
        auto f5 = dt * rhs(state_temp);

        state_temp = state_old - (8./27.) * f1 + 2. * f2 - (3544./2565.) * f3 + (1859./4104.) * f4 - (11./40.) * f5;
        auto f6 = dt * rhs(state_temp);

        // final solution

        auto state4_new = state_old + (25./216.) * f1 + (1408./2565.) * f3 + (2197./4101.) * f4 - 0.2 * f5;
        auto state5_new = state_old + (16./135.) * f1 + (6656./12825.) * f3 + (28561./56430.) * f4 - (9./50.) * f5 + (2./55) * f6;

        return {state4_new, state5_new};
    }


    void integrate(const double dt_in, const double err, const double tmax) {
        // integrate the equations of motion using 4th order R-K method with an
        // adaptive stepsize, to try to achieve the relative error err.  dt here
        // is the initial timestep

        // if err < 0, then we don't do adaptive stepping, but rather
        // we always walk at the input dt

        const double atol{1.e-10};

        // start with the old timestep
        double dt_new{dt_in};
        double dt{dt_in};

        int n_reset{0};
        double t{0.0};

        while (t < tmax) {

            auto &state_old = history.back();

            // adaptive stepping iteration loop -- keep trying to take
            // a step until we achieve our desired error
            double rel_error = std::numeric_limits<double>::max();

            int n_try{0};

            OrbitState state_new;

            int iter{0};
            while (rel_error > 1) {
                dt = dt_new;
                if (t + dt > tmax) {
                    dt = tmax - t;
                }

                // do the embedded RK45

                OrbitState state4;
                std::tie(state4, state_new) = single_step(state_old, dt);

                // state5 should be more accurate than state4 since it
                // is 5th order

                // estimate the error now.
                rel_error = std::max({std::abs(state_new.x - state4.x) / (err * std::abs(state_old.x) + atol),
                                      std::abs(state_new.y - state4.y) / (err * std::abs(state_old.y) + atol),
                                      std::abs(state_new.u - state4.u) / (err * std::abs(state_old.u) + atol),
                                      std::abs(state_new.v - state4.v) / (err * std::abs(state_old.v) + atol)});
                std::cout << "relative error = " << rel_error << std::endl;

                // RK45 new timestep

                dt_new = dt * 0.84 * std::pow(std::abs(1.0/rel_error), 0.25);

                n_try++;
                iter++;
            }

            if (n_try > 1) {
                // n_try = 1 if we took only a single try at the step
                n_reset += (n_try - 1);
            }

            // successful step
            t += dt;

            std::cout << t << " " << dt << std::endl;

            // store the new solution -- this will automatically be used as the
            // "old" state in the next step

            time.push_back(t);
            history.push_back(state_new);

        }

    }
};


int main() {

    OrbitsRK4 o(1.0, 0.95);
    o.integrate(0.05, 1.e-5, 1.0);

    std::cout << "number of points = " << o.npts() << std::endl;

    std::ofstream of("orbit_adaptive.dat");

    for (int n = 0; n < o.npts(); ++n) {
        auto& state = o.get_state(n);

        of << std::setw(14) << o.get_time(n);

        of << state;

        of.precision(10);

        of << std::setw(18) << o.energy(n);

        of << std::endl;

    }

}
