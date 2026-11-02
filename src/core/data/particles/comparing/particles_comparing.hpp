#ifndef PHARE_CORE_DATA_PARTICLES_COMPARING_PARTICLES_COMPARING
#define PHARE_CORE_DATA_PARTICLES_COMPARING_PARTICLES_COMPARING

#include "core/data/particles/particle_array_def.hpp"

#include "core/data/particles/comparing/detail/def_comparing.hpp"
#include "core/data/particles/comparing/detail/aos_comparing.hpp"
#include "core/data/particles/comparing/detail/soa_comparing.hpp"

namespace PHARE::core
{


template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
template<typename PS0, typename PS1>
EqualityReport
ParticlesComparator<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
    PS0 const& ps0, PS1 const& ps1, double const atol)
{
    return particles_equals(ps0, ps1, atol);
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_COMPARING_PARTICLES_COMPARING */
