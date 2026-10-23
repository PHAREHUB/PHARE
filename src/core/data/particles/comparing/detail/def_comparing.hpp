
#ifndef PHARE_CORE_DATA_PARTICLES_COMPARING_DETAIL_DEF_COMPARING
#define PHARE_CORE_DATA_PARTICLES_COMPARING_DETAIL_DEF_COMPARING

#include "core/utilities/equality.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE::core
{


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


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_COMPARING_DETAIL_DEF_COMPARING */
