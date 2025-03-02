// this driver is used to create a sample radial velocity
// curve for an FFT example

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "solar_system.H"

int main() {

    constexpr double M_star{0.3};

    std::vector<double> m_planets{1.e-3, 4.e-5, 3.e-5, 1.e-5, 0.5e-5};
    std::vector<double> a_planets{2, 1, 0.5, 0.2, 0.1};

    SolarSystem s(M_star, m_planets, a_planets, 1000);

    for (const auto& states : s.history) {
        for (const auto& o : states) {
            std::cout << o << std::endl;
        }
    }

    std::cout << "initial hill separation = "
              << s.compute_hill_separation(s.history[0]) << std::endl;

    double max_P = std::sqrt(std::pow(a_planets[0], 3) / M_star);
    std::cout << "maximum period = " << max_P << std::endl;

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
}
