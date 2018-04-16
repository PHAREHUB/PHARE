#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H

#include "particle.h"

#include <cstddef>
#include <vector>


namespace PHARE
{
template<std::size_t dim>
using ParticleArray = std::vector<Particle<dim>>;
}


#endif
