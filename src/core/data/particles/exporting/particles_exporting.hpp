#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_PARTICLES_EXPORTING
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_PARTICLES_EXPORTING

#include "core/utilities/memory.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"
#include "core/data/particles/exporting/detail/aos_exporting.hpp"
// #include "core/data/particles/exporting/detail/soa_exporting.hpp"

namespace PHARE::core
{


// // generic fallthrough
// template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
// template<auto type, typename Src, typename Dst>
// void ParticlesAppender<src_layout_mde, src_alloc_mde, dst_layout_mde, dst_alloc_mde>::operator()(
//     Src const& src, Dst& dst)
// {
//     // not optimized but *should work*
//     dst.reserve(dst.size() + src.size());
//     std::copy(src.begin(), src.end(), std::back_inserter(dst));
// }



} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_PARTICLES_EXPORTING */
