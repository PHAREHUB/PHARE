#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_H


#include <cstddef>
#include <vector>

#include "particle.h"

namespace PHARE
{
namespace core
{
    // TODO make a real particleArray class that has copy-deleted Ctor
    template<std::size_t dim>
    using ParticleArray = std::vector<Particle<dim>>;
} // namespace core
} // namespace PHARE


#endif
