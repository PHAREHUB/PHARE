#ifndef PHARE_CORE_DATA_PARTICLES_SERIALIZING_PARTICLES_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_SERIALIZING_PARTICLES_SERIALIZING

#include "core/data/particles/serializing/detail/def_serializing.hpp"
#include "core/data/particles/serializing/detail/aos_serializing.hpp"
#include "core/data/particles/serializing/detail/soa_serializing.hpp"

namespace PHARE::core
{



// template<auto src_layout_mde, auto src_alloc_mde>
// template<typename Dst, typename Src>
// void ParticlesDeserializer<src_layout_mde, src_alloc_mde>::operator()(std::string const&
// file_name,
//                                                                       Dst& dst)
// {
//     ParticlesDeserializer<Src::layout_mode, CPU>{}.template readN<Src>(
//         file_name, [&](auto const& particles) {
//             for (auto const& p : particles)
//                 dst.push_back(p);
//         });
// }

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_SERIALIZING_PARTICLES_SERIALIZING */
