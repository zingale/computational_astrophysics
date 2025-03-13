// this driver is used to create a sample radial velocity
// curve for an FFT example

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ranges>
#include <vector>
#include <print>

#include "solar_system.H"

int main() {

    constexpr double M_star{0.3};

    std::vector<double> m_planets{1.e-3, 4.e-5, 3.e-5, 1.e-5}; //, 0.5e-5};
    std::vector<double> a_planets{2, 1, 0.5, 0.2}; // , 0.1};

    SolarSystem s(M_star, m_planets, a_planets, 1000);

    for (const auto& states : s.history) {
        for (const auto& o : states) {
            std::cout << o << std::endl;
        }
    }

    std::cout << "initial hill separation = "
              << s.compute_hill_separation(s.history[0]) << std::endl;

    double max_P{-1.0};
    for (const auto [n, a] : a_planets | std::views::enumerate) {
        double P = std::sqrt(std::pow(a, 3) / M_star);
        max_P = std::max(max_P, P);
        std::println("planet {} period = {:7.3f} years", n, P);
    }

    // observe once a week
    s.integrate(1.0/52.0, 25 * max_P);

    std::ofstream of("radial_velocity.dat");

    for (std::size_t n = 0; n < s.history.size(); ++n) {
        const auto& states = s.history[n];

        // only output the radial velocity of the star
        // (which we'll take as just the x-velocity)
        of << std::setw(14) << s.times[n];
        of << std::setw(14) << states[0].u << std::endl;
    }

    std::ofstream of2("planet_positions.dat");

    for (std::size_t n = 0; n < s.history.size(); ++n) {
        const auto& states = s.history[n];

        // only output the radial velocity of the star
        // (which we'll take as just the x-velocity)
        of2 << std::setw(14) << s.times[n];
        for (const auto [n, p] : states | std::views::enumerate) {
            if (n == 0) {
                // skip the star
                continue;
            }

            of2 << std::setw(14) << p.x << " ";
            of2 << std::setw(14) << p.y << " ";
        }
        of2 << std::endl;
    }

}
