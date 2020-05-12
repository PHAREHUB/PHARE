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
    template<std::size_t dim, typename Float = double>
    using ParticleArray = std::vector<Particle<dim, Float>>;

    template<std::size_t dim, typename Float = double>
    void empty(ParticleArray<dim, Float>& array)
    {
        array.erase(std::begin(array), std::end(array));
    }

    template<std::size_t dim, typename Float = double>
    void swap(ParticleArray<dim, Float>& array1, ParticleArray<dim, Float>& array2)
    {
        std::swap(array1, array2);
    }

} // namespace core
} // namespace PHARE


#endif
