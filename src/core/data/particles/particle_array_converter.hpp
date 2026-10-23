#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CONVERTER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CONVERTER

#include "core/data/particles/particle_array_detail.hpp"
#include "core/data/particles/particle_array_sorter.hpp"

namespace PHARE::core
{

// --- public API ---

// hides the per layout/alloc impl behind operator() - see implementation section below
template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesConverter;


template<typename Dst, typename Src, typename GridLayout>
auto static convert_particles_from(Src const& src, GridLayout const& layout)
{
    using Converter
        = ParticlesConverter<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>;

    return Converter{}.template operator()<Dst>(src, layout);
}

template<typename Dst, typename Src, typename GridLayout>
auto static convert_particles(Src const& src, GridLayout const& layout)
{
    if constexpr (std::is_same_v<Dst, Src>)
        return src;
    else
        return convert_particles_from<Dst>(src, layout);
}

template<typename Dst, typename Src, typename GridLayout>
auto static convert_particles_and_sort(Src const& src, GridLayout const& layout)
{
    auto out = convert_particles<Dst>(src, layout);
    sort_particles(out, grow(layout.AMRBox(), GridLayout::options.particle_ghost_width));
    return out;
}

template<auto layout_mode, auto O, typename GridLayout>
auto convert_to(ParticleArray<O> const& src, GridLayout const& layout)
{
    using Parts = ParticleArray<O.with_layout(layout_mode).with_storage(StorageMode::VECTOR)>;

    return convert_particles<Parts>(src, layout);
}


// --- implementations ---

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesConverter
{
    static_assert(all_are<LayoutMode>(src_layout_mde, dst_layout_mde));
    static_assert(all_are<AllocatorMode>(src_alloc_mde, dst_alloc_mde));

    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;

    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<typename Dst, typename Src, typename GridLayout>
    Dst operator()(Src const& src, GridLayout const& layout);
};


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
Dst ParticlesConverter<AoSPCTS, CPU, AoS, CPU>::operator()(Src const& src, GridLayout const& layout)
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
Dst ParticlesConverter<AoSCMTS, CPU, AoS, CPU>::operator()(Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());
    for (auto const& tile : src()) // one flat AoSMapped array per tile
        std::copy(tile().begin(), tile().end(), std::back_inserter(out));
    return out;
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_CONVERTER */
