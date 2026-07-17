// IWYU pragma: private, include "core/data/particles/converting/particles_converting.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_AOS_CONVERTING
#define PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_AOS_CONVERTING

#include "core/utilities/memory.hpp"
#include "core/data/particles/converting/detail/def_converting.hpp"

namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;

template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoS, CPU, AoS, CPU>::operator()(Src const& src, GridLayout const& layout)
{
    return src;
}

template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSTS, CPU, AoS, CPU>::operator()(Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
        std::copy(tile().begin(), tile().end(), std::back_inserter(out));

    return out;
}
template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSCMTS, CPU, AoS, CPU>::operator()(Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());
    for (auto const& tile : src())
        std::copy(tile().begin(), tile().end(), std::back_inserter(out));
    return out;
}

template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSPCTS, CPU, AoS, CPU>::operator()(Src const& src,
                                                            GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());
    for (auto const& tile : src())
    {
        // full per-tile box: particles that left the patch domain live in tile ghost
        // cells; tile ghost cells inside the patch stay empty (owner tiles hold those)
        auto const& cps = tile();
        for (auto const& bix : cps.local_box())
            std::copy(cps(bix).begin(), cps(bix).end(), std::back_inserter(out));
    }
    return out;
}

template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<AoSTS, GPU_UNIFIED, AoS, CPU>::operator()( //
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.resize(src.size());
    std::size_t offset = 0;

    for (auto const& tile : src())
    {
        auto& particles = tile();
        if (particles.size() == 0)
            continue;
        gpu::copy(out.data() + offset, particles.data(), particles.size());
        offset += particles.size();
    }
    return out;
}




} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_AOS_CONVERTING */
