


#include "data/ions/particle_initializers/fluid_particle_initializer.h"


void maxwellianVelocity(std::array<double, 3> V, std::array<double, 3> Vth,
                        std::mt19937_64 generator, std::array<double, 3>& partVelocity)
{
    std::normal_distribution<> maxwellX(V[0], Vth[0]);
    std::normal_distribution<> maxwellY(V[1], Vth[1]);
    std::normal_distribution<> maxwellZ(V[2], Vth[2]);

    partVelocity[0] = maxwellX(generator);
    partVelocity[1] = maxwellY(generator);
    partVelocity[2] = maxwellZ(generator);
}

