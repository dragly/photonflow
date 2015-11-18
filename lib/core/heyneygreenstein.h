#ifndef HEYNEYGREENSTEIN_H
#define HEYNEYGREENSTEIN_H

#include "../core/randomnumbergenerator.h"

namespace Distribution {
inline double heyneyGreenstein(double g, RNG &rng)
{
    double g2 = g*g;
    double eta = rng.RandomFloat();
    double cosTheta = 0.0;
    if(g > 1e-2) {
        double k = ((1 - g2) / (1 - g + 2.0 * g * eta));
        double k2 = k*k;
        cosTheta = 1.0 / (2.0 * g) * (1.0 + g2 - k2);
    } else {
        cosTheta = 2*eta - 1;
    }
    return cosTheta;
}
}

#endif // HEYNEYGREENSTEIN_H

