
#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_DEF_REFINING
#define PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_DEF_REFINING

#include "core/data/particles/particle_array_def.hpp"


namespace PHARE::amr
{

template<auto layout_mde, auto alloc_mde>
struct ParticlesRefiner
{
    static_assert(core::all_are<core::LayoutMode>(layout_mde));
    static_assert(core::all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    template<auto type, typename Src, typename Dst, typename Box_t, typename Refiner,
             typename Transformer>
    void operator()(Src const&, Dst&, Box_t const&, Refiner, Transformer);
};


} // namespace PHARE::amr


#endif /* PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_DEF_REFINING */
