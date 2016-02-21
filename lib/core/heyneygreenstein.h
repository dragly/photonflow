#ifndef HEYNEYGREENSTEIN_H
#define HEYNEYGREENSTEIN_H

#include "../core/randomnumbergenerator.h"
#include "../core/common.h"

namespace Distribution {
inline double heyneyGreenstein(double g, RNG &rng)
{
    double eps = 1e-16;
    double epsi = (1 - 2*eps); // used to avoid 0 and 1
    double g2 = g*g;
    double eta = eps + epsi * rng.randomFloat();
    double cosTheta = 0.0;
    if(g > 1e-2) {
        double k = ((1 - g2) / (1 - g + 2.0 * g * eta));
        double k2 = k*k;
        cosTheta = 1.0 / (2.0 * g) * (1.0 + g2 - k2);
    } else {
        cosTheta = 2*eta - 1;
    }
    photonflowAssert(cosTheta >= -1.0 && cosTheta <= 1.0);
    return cosTheta;
}
}

#endif // HEYNEYGREENSTEIN_H

