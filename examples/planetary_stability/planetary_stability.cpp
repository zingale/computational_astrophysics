#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <vector>
#include <random>
#include <cassert>
#include <deque>
#include <numeric>

constexpr double G_const{4.0 * M_PI * M_PI};

struct OrbitState {
    // a container to hold the star positions
    double mass;
    double x;
    double y;
    double u;
    double v;

    OrbitState(double M,
               double x0, double y0, double u0, double v0)
        : mass(M), x(x0), y(y0), u(u0), v(v0)
    {}

    OrbitState()
        : mass(0.0), x(0.0), y(0.0), u(0.0), v(0.0)
    {}
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

class SolarSystem {

    double M_star{};
    double escape_R{};

    // Yoshida integration constants
    std::array<double, 4> c{};
    std::array<double, 3> d{};

public:

    // for plotting, we'll store all the initial semi-major
    // axes, including the star's
    std::vector<double> a_initial;

    std::deque<std::vector<OrbitState>> history;
    std::deque<double> times;

    int n_objects{};

    const double SMALL{1.e-30};

    SolarSystem(double _M_star,
                const std::vector<double>& m_planets,
                const std::vector<double>& a_planets,
                double escape_radius=10000) :
        M_star(_M_star),
        escape_R(escape_radius)
    {

        assert (m_planets.size() == a_planets.size());

        // integration constants
        double w0 = -std::cbrt(2.0) / (2.0 - cbrt(2.0));
        double w1 = 1.0 / (2.0 - std::cbrt(2.0));

        c[0] = w1 / 2.0;
        c[1] = 0.5 * (w0 + w1);
        c[2] = c[1];
        c[3] = c[0];

        d[0] = w1;
        d[1] = w0;
        d[2] = w1;

        a_initial.push_back(0.0);
        std::copy(a_planets.cbegin(), a_planets.cend(),
                  std::back_inserter(a_initial));

        n_objects = a_initial.size();

        // create the initial solar system
        std::vector<OrbitState> system;

        // store the star
        system.emplace_back(OrbitState(M_star, 0.0, 0.0, 0.0, 0.0));

        // now the planets.  We will give them random phases

        // seed the random number generator
        std::mt19937 generator(123);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);

        for (std::size_t n = 0; n < m_planets.size(); ++n) {
            double phi = uniform(generator) * 2.0 * M_PI;
            double v_circ = std::sqrt(G_const * M_star / a_planets[n]);
            system.emplace_back(OrbitState(m_planets[n],
                                           a_planets[n] * std::cos(phi),
                                           a_planets[n] * std::sin(phi),
                                           -v_circ * std::sin(phi),
                                           v_circ * std::cos(phi)));

        }

        // now correct the star's state so the center of mass
        // is stationary and at the origin

        // note: since the star was initialized with all 0s,
        // we can include it in the sum here without issue
        double mx_sum{};
        double my_sum{};
        double mu_sum{};
        double mv_sum{};
        for (const auto& object : system) {
            mx_sum += object.mass * object.x;
            my_sum += object.mass * object.y;
            mu_sum += object.mass * object.u;
            mv_sum += object.mass * object.v;
        };
        system[0].x -= mx_sum / system[0].mass;
        system[0].y -= my_sum / system[0].mass;
        system[0].u -= mu_sum / system[0].mass;
        system[0].v -= mv_sum / system[0].mass;

        // store the initial conditions
        history.push_back(system);
        times.push_back(0.0);
    }

    double
    compute_hill_separation(const std::vector<OrbitState>& state) {
        double h_sep{std::numeric_limits<double>::max()};

        for (int p1 = 1; p1 < n_objects; ++p1) {
            double a1 = std::sqrt(std::pow(state[p1].x - state[0].x, 2) +
                                  std::pow(state[p1].y - state[0].y, 2));
            for (int p2 = 1; p2 < n_objects; ++p2) {
                if (p1 == p2) {
                    continue;
                }
                double a2 = std::sqrt(std::pow(state[p2].x - state[0].x, 2) +
                                      std::pow(state[p2].y - state[0].y, 2));

                double r_hill = 0.5 * (a1 + a2) *
                    std::cbrt((state[p1].mass + state[p2].mass) /
                              (3.0 * M_star));

                h_sep = std::min(h_sep, std::abs(a2 - a1) / r_hill);
            }
        }

        return h_sep;
    }

    std::vector<OrbitState>
    rhs(const std::vector<OrbitState>& states) {

        std::vector<OrbitState> ydots;

        for (int iobj = 0; iobj < n_objects; ++iobj) {

            double dxdt = states[iobj].u;
            double dydt = states[iobj].v;

            double dudt{};
            double dvdt{};

            for (int jobj = 0; jobj < n_objects; ++jobj) {
                if (iobj == jobj) {
                    continue;
                }
                double dx = states[jobj].x - states[iobj].x;
                double dy = states[jobj].y - states[iobj].y;

                double r = std::sqrt(dx * dx + dy * dy + SMALL);

                dudt += G_const * states[jobj].mass * dx / std::pow(r, 3);
                dvdt += G_const * states[jobj].mass * dy / std::pow(r, 3);

            }

            ydots.emplace_back(OrbitState(states[iobj].mass, dxdt, dydt, dudt, dvdt));
        }

        return ydots;
    }

    std::vector<OrbitState>
    single_step_yoshida(const std::vector<OrbitState>& states, const double dt) {

        // we can just work on the new state directly, so start by simply
        // copying the old state to the new state
        std::vector<OrbitState> states_new;
        for (const auto& s : states) {
            states_new.emplace_back(OrbitState(s.mass, s.x, s.y, s.u, s.v));
        }
        
        // drift with the current velocity
        for (auto& s : states_new) {
            s.x += c[0] * s.u * dt;
            s.y += c[0] * s.v * dt;
        }

        for (int i = 0; i < 3; ++i) {
            // evaluate the force
            auto ydots = rhs(states_new);

            // kick followed by drift
            for (int n = 0; n < n_objects; ++n) {
                states_new[n].u += d[i] * ydots[n].u * dt;
                states_new[n].v += d[i] * ydots[n].v * dt;

                states_new[n].x += c[i+1] * states_new[n].u * dt;
                states_new[n].y += c[i+1] * states_new[n].v * dt;
            }

        }

        return states_new;
    }

};


int main() {

    std::vector<double> m_planets{0.0054, 0.0074, 0.0071};
    std::vector<double> a_planets{71.6, 41.4, 16.3};

    SolarSystem s(1.5, m_planets, a_planets, 1000);

    for (const auto& h : s.history) {
        for (const auto& o : h) {
            std::cout << o << std::endl;
        }
    }

    std::cout << "initial hill separation = " << s.compute_hill_separation(s.history[0]) << std::endl;

}
