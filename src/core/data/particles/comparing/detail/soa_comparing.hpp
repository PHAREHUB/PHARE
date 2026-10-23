
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


// compares via per-index accessors only (iCell(i)/v(i)/delta(i)) - AoS and SoA both
// provide them, regardless of storage layout, and SoA has no flat iterator
template<typename PS0, typename PS1>
EqualityReport index_based_particles_equals(PS0 const& ps0, PS1 const& ps1, double const atol)
{
    if (ps0.size() != ps1.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ps0.size()) + " vs "
                                         + std::to_string(ps1.size())};

    for (std::size_t i = 0; i < ps0.size(); ++i)
    {
        std::string const idx = std::to_string(i);
        if (ps0.iCell(i) != ps1.iCell(i))
            return EqualityReport{false, "icell mismatch at index: " + idx, i};

        if (!float_equals(ps0.v(i), ps1.v(i), atol))
            return EqualityReport{false, "v mismatch at index: " + idx, i};

        if (!float_equals(ps0.delta(i), ps1.delta(i), atol))
            return EqualityReport{false, "delta mismatch at index: " + idx, i};
    }

    return EqualityReport{true};
}


template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<SoA, CPU, SoA, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                   double const atol)
{
    return index_based_particles_equals(ps0, ps1, atol);
}

template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<SoA, CPU, AoS, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                   double const atol)
{
    return index_based_particles_equals(ps0, ps1, atol);
}

template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<AoS, CPU, SoA, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                   double const atol)
{
    return index_based_particles_equals(ps0, ps1, atol);
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
