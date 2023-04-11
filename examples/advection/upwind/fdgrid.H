#ifndef FDGRID_H
#define FDGRID_H

#include <iostream>
#include <cassert>
#include <vector>

class Grid {
    // a finite-difference grid class with nx points and ng ghost points

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

    Grid(const int _nx, const int _ng, const double _xmin, const double _xmax)
        : nx(_nx), ng(_ng), xmin(_xmin), xmax(_xmax), nq(nx + 2*ng),
          ilo(ng), ihi(ng+nx-1),
          x(nq, 0.0), a(nq, 0.0)
    {
        assert(nx > 0);
        assert(ng >= 0);
        assert(xmax > xmin);

        dx = (xmax - xmin) / static_cast<double>(nx);

        for (int i = 0; i < nq; ++i) {
            x[i] = xmin + static_cast<double>(i - ng) * dx;
        }

    }

    void fill_BCs() {

        // periodic BCs

        // left edge of domain
        for (int i = 0; i < ng; ++i) {
            a[ilo-1-i] = a[ihi-i];
        }

        // right edge of domain
        for (int i = 0; i < ng; ++i) {
            a[ihi+1+i] = a[ilo+i];
        }
    }

};
#endif
