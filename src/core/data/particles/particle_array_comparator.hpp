#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_COMPARATOR
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_COMPARATOR

#include "core/data/particles/comparing/particles_comparing.hpp"

namespace PHARE::core
{


template<typename PS0, typename PS1>
auto static compare_particles(PS0 const& ps0, PS1 const& ps1, double const atol = 1e-15)
{
    using Comparator
        = ParticlesComparator<PS0::layout_mode, PS0::alloc_mode, PS1::layout_mode, PS1::alloc_mode>;

    return Comparator{}(ps0, ps1, atol);
}




} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_COMPARATOR */
