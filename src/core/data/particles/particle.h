#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_H


#include <array>


struct Particle
{
    double weight{0};
    double charge{0};

    double Ex = 0, Ey = 0, Ez = 0; //! electric field seen by the particle
    double Bx = 0, By = 0, Bz = 0; //! magnetic field seen by the particle

    std::array<double, 3> v    = {{0, 0, 0}}; //! velocity
    std::array<float, 3> delta = {{0, 0, 0}}; //! normalized position within a cell (in [0,1[)
    std::array<int, 3> iCell   = {{0, 0, 0}}; //! cell index in the three directions


    Particle()                       = default;
    Particle(Particle const& source) = default;
    Particle(Particle&& source)      = default;
    Particle& operator=(Particle const& source) = default;
    Particle& operator=(Particle&& source) = default;

    Particle(double weight, double charge, std::array<int, 3> iCell, std::array<float, 3> delta,
             std::array<double, 3> v)
        : weight{weight}
        , charge{charge}
        , v{v}
        , delta{delta}
        , iCell{iCell}
    {
    }
};



#endif
