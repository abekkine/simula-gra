#ifndef UNIVERSE_H_
#define UNIVERSE_H_

#include "display.h"

namespace universe {

    enum DistributionType {
        D_ROTATING_SPHERICAL = 1001,
        D_RANDOM_SPHERICAL,
        D_RANDOM_CUBICAL,
        D_UNIFORM_CUBICAL,
        D_UNIFORM_CUBICAL_W_NOISE
    };

    extern int NumBodies;
    extern double SpreadRadius;
    extern int MassDistribution;

    void init();
    void update();
    void render();
}

#endif // UNIVERSE_H_

