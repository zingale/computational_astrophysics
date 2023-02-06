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
};

std::ostream& operator<< (std::ostream& os, const OrbitState& s) {
    os.precision(6);

    os << std::setw(14) << s.x
       << std::setw(14) << s.y
       << std::setw(14) << s.u
       << std::setw(14) << s.v;

    return os;
}

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
        // return a referece to the state at time index n
        return history[n];
    }

    double energy(const int n) {
        // return the energy of the system for timestep n

        auto& state = history[n];

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

    OrbitState update_state(const OrbitState& state_old, const OrbitState& ydot, const double dt) {

        OrbitState s;

        s.x = state_old.x + dt * ydot.x;
        s.y = state_old.y + dt * ydot.y;
        s.u = state_old.u + dt * ydot.u;
        s.v = state_old.v + dt * ydot.v;

        return s;
    }

    OrbitState single_step(const OrbitState& state_old, const double dt) {
        /// take a single RK-4 timestep through dt

        auto ydot1 = rhs(state_old);

        OrbitState state_temp;

        state_temp = update_state(state_old, ydot1, 0.5*dt);
        auto ydot2 = rhs(state_temp);

        state_temp = update_state(state_old, ydot2, 0.5*dt);
        auto ydot3 = rhs(state_temp);

        state_temp = update_state(state_old, ydot3, dt);
        auto ydot4 = rhs(state_temp);

        // final solution

        OrbitState s;

        s.x = state_old.x + (dt / 6.0) * (ydot1.x + 2.0 * ydot2.x + 2.0 * ydot3.x + ydot4.x);
        s.y = state_old.y + (dt / 6.0) * (ydot1.y + 2.0 * ydot2.y + 2.0 * ydot3.y + ydot4.y);
        s.u = state_old.u + (dt / 6.0) * (ydot1.u + 2.0 * ydot2.u + 2.0 * ydot3.u + ydot4.u);
        s.v = state_old.v + (dt / 6.0) * (ydot1.v + 2.0 * ydot2.v + 2.0 * ydot3.v + ydot4.v);

        return s;
    }


    void integrate(const double dt_in, const double err, const double tmax) {
        // integrate the equations of motion using 4th order R-K method with an
        // adaptive stepsize, to try to achieve the relative error err.  dt here
        // is the initial timestep

        // if err < 0, then we don't do adaptive stepping, but rather
        // we always walk at the input dt

        // safety parameters

        const double S1{0.9};
        const double S2{4.0};

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

            while (rel_error > err) {
                dt = dt_new;
                if (t + dt > tmax) {
                    dt = tmax - t;
                }

                // take 2 half steps

                auto state_temp = single_step(state_old, 0.5*dt);
                state_new = single_step(state_temp, 0.5*dt);

                // now take just a single step to cover dt

                auto state_single = single_step(state_old, dt);

                // state_new should be more accurate than state_single since it
                // used smaller steps

                // estimate the relative error now

                double pos_err = std::max(std::abs((state_new.x - state_single.x) / state_single.x),
                                          std::abs((state_new.y - state_single.y) / state_single.y));

                double vel_err = std::max(std::abs((state_new.u - state_single.u) / state_single.u),
                                          std::abs((state_new.v - state_single.v) / state_single.v));

                rel_error = std::max(pos_err, vel_err);

                // adaptive timestep algorithm from Garcia (Eqs. 3.30 and 3.31)

                double dt_est = dt * std::pow(std::abs(err/rel_error), 0.2);
                dt_new = std::min(std::max(S1*dt_est, dt/S2), S2*dt);

                n_try++;
            }

            if (n_try > 1) {
                // n_try = 1 if we took only a single try at the step
                n_reset += (n_try - 1);
            }

            // successful step
            t += dt;

            // store the new solution -- this will automatically be used as the
            // "old" state in the next step

            time.push_back(t);
            history.push_back(state_new);

        }

    }
};


int main() {

    OrbitsRK4 o(1.0, 0.95);
    o.integrate(0.05, 1.e-8, 1.0);

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
