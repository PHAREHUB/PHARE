#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_COMPARATOR
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_COMPARATOR

#include "core/utilities/equality.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <string>
#include <cstddef>

namespace PHARE::core
{

// --- public API ---

// hides the per layout/alloc impl behind operator() - see implementation section below
template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesComparator;


template<typename PS0, typename PS1>
auto static compare_particles(PS0 const& ps0, PS1 const& ps1, double const atol = 1e-15)
{
    using Comparator
        = ParticlesComparator<PS0::layout_mode, PS0::alloc_mode, PS1::layout_mode, PS1::alloc_mode>;

    return Comparator{}(ps0, ps1, atol);
}


// --- implementations ---

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesComparator
{
    static_assert(all_are<LayoutMode>(src_layout_mde, dst_layout_mde));
    static_assert(all_are<AllocatorMode>(src_alloc_mde, dst_alloc_mde));

    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;
    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<typename PS0, typename PS1>
    EqualityReport operator()(PS0 const& ps0, PS1 const& ps1, double const atol = 1e-15);
};


// shared by any layout combination whose particle arrays expose per-index accessors
// (iCell(i)/v(i)/delta(i)) - AoS and SoA both do, regardless of storage layout
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


template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
template<typename PS0, typename PS1>
EqualityReport
ParticlesComparator<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
    PS0 const& ps0, PS1 const& ps1, double const atol)
{
    return particles_equals(ps0, ps1, atol);
}


using enum LayoutMode;
using enum AllocatorMode;

template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<AoS, CPU, AoS, CPU>::operator()(PS0 const& ps0, PS1 const& ps1,
                                                                   double const atol)
{
    return index_based_particles_equals(ps0, ps1, atol);
}

template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<AoSMapped, CPU, AoS, CPU>::operator()(PS0 const& ps0,
                                                                         PS1 const& ps1,
                                                                         double const atol)
{
    return ParticlesComparator<AoS, CPU, AoS, CPU>{}(ps0, ps1, atol);
}


template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<AoSPC, CPU, AoSPC, CPU>::operator()(PS0 const& ps0,
                                                                       PS1 const& ps1,
                                                                       double const atol)
{
    if (ps0.size() != ps1.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ps0.size()) + " vs "
                                         + std::to_string(ps1.size())};

    ParticlesComparator<AoS, CPU, AoS, CPU> comparator;
    for (auto const& cell : ps0.local_box())
        if (auto eq = comparator(ps0(cell), ps1(cell), atol); !eq)
            return eq;

    return EqualityReport{true};
}


template<>
template<typename PS0, typename PS1>
EqualityReport ParticlesComparator<AoSPCTS, CPU, AoSPCTS, CPU>::operator()(PS0 const& ps0,
                                                                           PS1 const& ps1,
                                                                           double const atol)
{
    if (ps0.size() != ps1.size())
        return EqualityReport{false, "different sizes: " + std::to_string(ps0.size()) + " vs "
                                         + std::to_string(ps1.size())};

    if (ps0().size() != ps1().size())
        return EqualityReport{false, "different tile counts: " + std::to_string(ps0().size())
                                         + " vs " + std::to_string(ps1().size())};

    ParticlesComparator<AoSPC, CPU, AoSPC, CPU> comparator;
    for (std::size_t ti = 0; ti < ps0().size(); ++ti)
        if (auto eq = comparator(ps0()[ti](), ps1()[ti](), atol); !eq)
            return eq;

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


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_COMPARATOR */
