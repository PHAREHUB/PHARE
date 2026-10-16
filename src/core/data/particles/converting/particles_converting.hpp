#ifndef PHARE_CORE_DATA_PARTICLES_CONVERTING_PARTICLES_CONVERTING
#define PHARE_CORE_DATA_PARTICLES_CONVERTING_PARTICLES_CONVERTING

#include "core/data/particles/particle_array_def.hpp"

#include "core/data/particles/converting/detail/def_converting.hpp"
#include "core/data/particles/converting/detail/aos_converting.hpp"
#include "core/data/particles/converting/detail/soa_converting.hpp"

namespace PHARE::core
{

// generic fallthrough
template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
    Src const& src, GridLayout const& layout)
{
    auto dst = make_particles<Dst>(layout);


    // not optimized but *should work*
    dst.reserve(dst.size() + src.size());
    std::copy(src.begin(), src.end(), std::back_inserter(dst));

    return dst;
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_CONVERTING_PARTICLES_CONVERTING */
