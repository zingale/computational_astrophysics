#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <limits>
#include <fstream>

const double G = 1.0;
const int N = 3;

struct State {
    // a container to hold the star positions
    double x{};
    double y{};
    double u{};
    double v{};
};

std::ostream& operator<< (std::ostream& os, const State& s) {
    os.precision(6);

    os << std::setw(14) << s.x
       << std::setw(14) << s.y
       << std::setw(14) << s.u
       << std::setw(14) << s.v;

    return os;
}

class ThreeBody {
    /// model the gravitational interaction of three stars

private:

    double SMALL;
    std::vector<double> masses;
    std::vector<double> time;
    std::vector<std::vector<State>> stars;

public:

    ThreeBody (const double M0, const std::vector<double>& pos0,
               const double M1, const std::vector<double>& pos1,
               const double M2, const std::vector<double>& pos2,
               const double _SMALL=1.e-30)
               : SMALL{_SMALL},
                 masses{M0, M1, M2}
    {

        // we will organize the data as stars[n][istar].{t,x,y,u,v}, where i
        // istar the index of the star (0, 1, 2) and n is the timestep index

        std::vector<State> current(N);

        current[0].x = pos0[0];
        current[0].y = pos0[1];
        current[0].u = 0.0;
        current[0].v = 0.0;

        current[1].x = pos1[0];
        current[1].y = pos1[1];
        current[1].u = 0.0;
        current[1].v = 0.0;

        current[2].x = pos2[0];
        current[2].y = pos2[1];
        current[2].u = 0.0;
        current[2].v = 0.0;

        stars.push_back(current);

        time.push_back(0.0);

    }

    int npts() {
        // return the number of integration points
        return stars.size();
    }

    double get_time(const int n) {
        // return the physical time for step n
        return time[n];
    }

    std::vector<State>& get_state(const int n) {
        // return a reference to the state at time index n
        return stars[n];
    }

    std::pair<double, double> center_of_mass(const int n) {
        // return the center of mass of the system at time n

        auto Mtot = std::accumulate(masses.cbegin(), masses.cend(), 0.0);

        auto& state = stars[n];

        double xcm{};
        double ycm{};

        for (int istar = 0; istar < N; ++istar) {
            xcm += masses[istar] * state[istar].x;
            ycm += masses[istar] * state[istar].y;
        }
        xcm /= Mtot;
        ycm /= Mtot;

        return {xcm, ycm};
    }

    double energy(const int n) {
        // return the energy of the system for timestep n

        auto& state = stars[n];

        // kinetic energy

        double KE = 0.0;
        for (int istar = 0; istar < N; ++istar) {
            KE += 0.5 * masses[istar] *
                (std::pow(state[istar].u, 2) +
                 std::pow(state[istar].v, 2));
        }

        // potential energy -- we need all pairs, but order doesn't matter

        double PE = 0.0;
        for (int istar = 0; istar < N; ++istar) {
            for (int jstar = 0; jstar < istar; ++jstar) {
                PE += -G * masses[istar] * masses[jstar] /
                    std::sqrt(std::pow(state[istar].x - state[jstar].x, 2) +
                              std::pow(state[istar].y - state[jstar].y, 2));
            }
        }

        return KE + PE;
    }

    std::vector<State> rhs(const std::vector<State>& star_states) {
        // compute the ydot terms

        std::vector<State> ydot;

        for (int istar = 0; istar < N; ++istar) {
            State star_ydot;

            star_ydot.x = star_states[istar].u;
            star_ydot.y = star_states[istar].v;

            // for the gravitational acceleration, we need the force of this
            // star with all others

            double a_x{};
            double a_y{};

            for (int jstar = 0; jstar < N; ++jstar) {

                // no self-interaction

                if (istar == jstar) {
                    continue;
                }

                double dx = star_states[jstar].x - star_states[istar].x;
                double dy = star_states[jstar].y - star_states[istar].y;

                double r = std::sqrt(dx * dx + dy * dy) + SMALL;

                a_x += G * masses[jstar] * dx / std::pow(r, 3);
                a_y += G * masses[jstar] * dy / std::pow(r, 3);

            }

            star_ydot.u = a_x;
            star_ydot.v = a_y;

            ydot.push_back(star_ydot);

        }

        return ydot;

    }

    std::vector<State> update_state(const std::vector<State>& state_old,
                                    const std::vector<State>& ydot,
                                    const double dt) {

        std::vector<State> state_new;

        for (int istar = 0; istar < N; ++istar) {
            State s;
            s.x = state_old[istar].x + dt * ydot[istar].x;
            s.y = state_old[istar].y + dt * ydot[istar].y;
            s.u = state_old[istar].u + dt * ydot[istar].u;
            s.v = state_old[istar].v + dt * ydot[istar].v;

            state_new.push_back(s);
        }

        return state_new;
    }

    std::vector<State> single_step(const std::vector<State>& state_old, const double dt) {
        /// take a single RK-4 timestep through dt

        auto ydot1 = rhs(state_old);

        std::vector<State> state_temp;

        state_temp = update_state(state_old, ydot1, 0.5*dt);
        auto ydot2 = rhs(state_temp);

        state_temp = update_state(state_old, ydot2, 0.5*dt);
        auto ydot3 = rhs(state_temp);

        state_temp = update_state(state_old, ydot3, dt);
        auto ydot4 = rhs(state_temp);

        // final solution

        std::vector<State> state_new;

        for (int istar = 0; istar < N; ++istar) {
            State s;

            s.x = state_old[istar].x + (dt / 6.0) * (ydot1[istar].x + 2.0 * ydot2[istar].x + 2.0 * ydot3[istar].x + ydot4[istar].x);
            s.y = state_old[istar].y + (dt / 6.0) * (ydot1[istar].y + 2.0 * ydot2[istar].y + 2.0 * ydot3[istar].y + ydot4[istar].y);
            s.u = state_old[istar].u + (dt / 6.0) * (ydot1[istar].u + 2.0 * ydot2[istar].u + 2.0 * ydot3[istar].u + ydot4[istar].u);
            s.v = state_old[istar].v + (dt / 6.0) * (ydot1[istar].v + 2.0 * ydot2[istar].v + 2.0 * ydot3[istar].v + ydot4[istar].v);

            state_new.push_back(s);
        }

        return state_new;
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

            auto &state_old = stars.back();

            // adaptive stepping iteration loop -- keep trying to take
            // a step until we achieve our desired error
            double rel_error = std::numeric_limits<double>::max();

            int n_try{0};

            std::vector<State> state_new;

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

                double pos_err = std::numeric_limits<double>::lowest();
                double vel_err = std::numeric_limits<double>::lowest();

                for (int istar = 0; istar < N; ++istar) {
                    pos_err = std::max(pos_err,
                                       std::abs((state_new[istar].x - state_single[istar].x) / state_single[istar].x));
                    pos_err = std::max(pos_err,
                                       std::abs((state_new[istar].y - state_single[istar].y) / state_single[istar].y));

                    vel_err = std::max(vel_err,
                                       std::abs((state_new[istar].u - state_single[istar].u) / state_single[istar].u));
                    vel_err = std::max(vel_err,
                                       std::abs((state_new[istar].v - state_single[istar].v) / state_single[istar].v));
                }

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
            stars.push_back(state_new);

        }

    }
};


int main() {

    std::vector<double> pos0{3.0, 1.0};
    std::vector<double> pos1{-1.0, -2.0};
    std::vector<double> pos2{-1.0, 1.0};

    double M0 = 150.0;
    double M1 = 200.0;
    double M2 = 250.0;

    ThreeBody tb(M0, pos0, M1, pos1, M2, pos2);
    tb.integrate(0.05, 1.e-8, 2.0);


    std::cout << "number of points = " << tb.npts() << std::endl;

    std::ofstream of("three_body.dat");

    for (int n = 0; n < tb.npts(); ++n) {
        auto& state = tb.get_state(n);

        of << std::setw(14) << tb.get_time(n);

        // first the positions and velocities

        for (auto& star : state) {
            of << star;
        }

        // now the energy

        of.precision(10);

        of << std::setw(18) << tb.energy(n);

        // and the center of mass

        auto [xcm, ycm] = tb.center_of_mass(n);

        of << std::setw(18) << xcm << std::setw(18) << ycm;

        of << std::endl;
    }

}
