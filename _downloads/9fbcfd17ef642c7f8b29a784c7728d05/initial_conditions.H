#ifndef INITIAL_CONDITIONS_H
#define INITIAL_CONDITIONS_H

#include <fvgrid.H>

inline
void tophat(FVGrid& g) {
    // tophat initial conditions

    for (int i = g.ilo; i <= g.ihi; ++i) {
        if (g.x[i] > 1.0 / 3.0 && g.x[i] <= 2.0 / 3.0) {
            g.a[i] = 1.0;
        } else {
            g.a[i] = 0.0;
        }
    }

}
#endif
