// IWYU pragma: private, include "core/data/particles/particle_array_appender.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_PARTICLES_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_PARTICLES_APPENDING

// #include "core/utilities/memory.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"
#include "core/data/particles/appending/detail/aos_appending.hpp"
#include "core/data/particles/appending/detail/soa_appending.hpp"

namespace PHARE::core
{

// generic fallthrough
template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
    Src const& src, Dst& dst)
{
    // not optimized but *should work*
    // if it doesn't, make a new impl
    dst.reserve(dst.size() + src.size());
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}




} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_PARTICLES_APPENDING */
