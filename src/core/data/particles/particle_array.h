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

    template<std::size_t dim>
    void empty(ParticleArray<dim>& array)
    {
        array.erase(std::begin(array), std::end(array));
    }

    template<std::size_t dim>
    void swap(ParticleArray<dim>& array1, ParticleArray<dim>& array2)
    {
        std::swap(array1, array2);
    }

} // namespace core
} // namespace PHARE


#endif
