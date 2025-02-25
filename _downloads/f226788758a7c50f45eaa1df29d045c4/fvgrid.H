#ifndef FDGRID_H
#define FDGRID_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>

class FVGrid {
    // a finite-volume class with nx points and ng ghost points

    public:

    int nx;
    int ng;
    int xmin;
    int xmax;

    int nq;

    double dx;

    // indices of the valid domain (i.e. no ghost cells)
    int ilo;
    int ihi;

    // storage for the coordinate
    std::vector<double> x;

    // storage for the solution
    std::vector<double> a;

    FVGrid(const int _nx, const int _ng, const double _xmin, const double _xmax)
        : nx(_nx), ng(_ng), xmin(_xmin), xmax(_xmax), nq(nx + 2*ng),
          ilo(ng), ihi(ng+nx-1),
          x(nq, 0.0), a(nq, 0.0)
    {
        assert(nx > 0);
        assert(ng >= 0);
        assert(xmax > xmin);

        dx = (xmax - xmin) / static_cast<double>(nx);

        for (int i = 0; i < nq; ++i) {
            x[i] = xmin + static_cast<double>(i - ng + 0.5) * dx;
        }

    }

    void fill_BCs(std::vector<double>& adummy) {

        // periodic BCs

        // left edge of domain
        for (int i = 0; i < ng; ++i) {
            adummy[ilo-1-i] = adummy[ihi-i];
        }

        // right edge of domain
        for (int i = 0; i < ng; ++i) {
            adummy[ihi+1+i] = adummy[ilo+i];
        }
    }

    double norm(const std::vector<double>& e) {
        assert(static_cast<int>(e.size()) == nx + 2*ng);

        double l2{0.0};
        for (int i = ilo; i <= ihi; ++i) {
            l2 += e[i] * e[i];
        }
        l2 = std::sqrt(dx * l2);

        return l2;
    }

};
#endif
