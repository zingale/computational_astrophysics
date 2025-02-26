#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "solar_system.H"

int main() {

    std::vector<double> m_planets{0.0054, 0.0074, 0.0071};
    std::vector<double> a_planets{71.6, 41.4, 16.3};

    SolarSystem s(1.5, m_planets, a_planets, 1000);

    for (const auto& states : s.history) {
        for (const auto& o : states) {
            std::cout << o << std::endl;
        }
    }

    std::cout << "initial hill separation = "
              << s.compute_hill_separation(s.history[0]) << std::endl;

    s.integrate(0.05, 50000);

    std::ofstream of("planetary_stability.dat");

    for (std::size_t n = 0; n < s.history.size(); ++n) {
        const auto& states = s.history[n];

        of << std::setw(14) << s.times[n];
        for (const auto& object : states) {
            of << object;
        }
        of << std::endl;
    }

}
