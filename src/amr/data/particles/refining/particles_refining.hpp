#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_PARTICLES_REFINING
#define PHARE_AMR_DATA_PARTICLES_REFINING_PARTICLES_REFINING

#include "amr/data/particles/refining/detail/def_refining.hpp"
#include "amr/data/particles/refining/detail/aos_refining.hpp"
// #include "core/data/particles/refining/detail/soa_refining.hpp"

namespace PHARE::amr
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



} // namespace PHARE::amr

#endif /* PHARE_AMR_DATA_PARTICLES_REFINING_PARTICLES_REFINING */
