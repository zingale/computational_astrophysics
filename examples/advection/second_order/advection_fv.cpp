#include <iostream>
#include <functional>

#include <fvgrid.H>
#include <initial_conditions.H>

std::vector<double> flux_update(const FVGrid& g, const double u, const std::vector<double>& a) {
    // compute -div{F} for linear advection

    // first compute the slope

    std::vector<double> da(a.size(), 0.0);

    for (int i = g.ilo-1; i <= g.ihi+1; ++i) {
        // left slope
        double dl = a[i] - a[i-1];

        // right slope
        double dr = a[i+1] - a[i];

        // minmod
        if (dl * dr > 0.0) {
            da[i] = std::abs(dl) < std::abs(dr) ? dl : dr;
        } else {
            da[i] = 0.0;
        }
    }

    // now we do the Riemann problem and the state construction all at once,
    // since upwinding always takes the left state

    std::vector<double> aint(a.size(), 0.0);

    for (int i = g.ilo; i <= g.ihi+1; ++i) {
        aint[i] = a[i-1] + 0.5 * da[i-1];
    }

    // finally compute the flux difference

    std::vector<double> flux_diff(a.size(), 0.0);

    for (int i = g.ilo; i <= g.ihi; ++i) {
        flux_diff[i] = u * (aint[i] - aint[i+1]) / g.dx;
    }

    return flux_diff;

}

FVGrid advection_mol(const int nx, const double u, const double C,
                     const int num_periods, std::function<void(FVGrid&)> init_cond) {

    FVGrid g(nx, 2, 0.0, 1.0);

    // time information

    double dt = C * g.dx / u;
    double t{0.0};

    double tmax = num_periods * (g.xmax - g.xmin) / u;

    init_cond(g);

    // evolution loop

    std::vector<double> a_new(g.nq, 0.0);

    while (t < tmax) {

        if (t + dt > tmax) {
            dt = tmax - t;
        }

        // second-order RK

        g.fill_BCs(g.a);

        auto k1 = flux_update(g, u, g.a);

        std::vector<double> atmp(g.nq, 0.0);
        for (int i = g.ilo; i <= g.ihi; ++i) {
            atmp[i] = g.a[i] + 0.5 * dt * k1[i];
        }

        g.fill_BCs(atmp);

        auto k2 = flux_update(g, u, atmp);

        for (int i = g.ilo; i <= g.ihi; ++i) {
            g.a[i] += dt * k2[i];
        }

        t += dt;

    }

    return g;
}


int main() {

    int nx{128};
    double u{1.0};
    double C{0.5};

    int num_periods{1};

    auto g = advection_mol(nx, u, C, num_periods, tophat);

    // output

    for (int i = g.ilo; i <= g.ihi; ++i) {
        std::cout << g.x[i] << "   " << g.a[i] << std::endl;
    }

}
