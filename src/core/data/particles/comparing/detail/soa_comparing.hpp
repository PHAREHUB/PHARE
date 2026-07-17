
#ifndef PHARE_CORE_DATA_PARTICLES_COMPARING_DETAIL_SOA_COMPARING
#define PHARE_CORE_DATA_PARTICLES_COMPARING_DETAIL_SOA_COMPARING

#include "core/utilities/memory.hpp"
#include "core/utilities/equality.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/comparing/detail/def_comparing.hpp"

namespace PHARE::core
{
using enum LayoutMode;
using enum AllocatorMode;


template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<SoA, CPU, SoA, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                   double const atol)
{
    return EqualityReport{true};
}

template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<SoA, CPU, AoS, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                   double const atol)
{
    return EqualityReport{};
}

template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<SoAVX, CPU, AoS, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                     double const atol)
{
    return EqualityReport{};
}
template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<AoS, CPU, SoAVX, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                     double const atol)
{
    return ParticlesComparator<SoAVX, CPU, AoS, CPU>{}(ps1, ps0, atol);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_COMPARING_DETAIL_SOA_COMPARING */
