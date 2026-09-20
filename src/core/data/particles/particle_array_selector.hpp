#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP

#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

// --- public API ---

// hides the per layout/alloc impl behind select()/count() - see implementation section below
template<auto layout_mde, auto alloc_mde>
struct ParticlesSelector;

template<ParticleType particle_type = ParticleType::Domain, typename Src, typename Dst,
         typename box_t>
void select_particles(Src const& src, Dst& dst, box_t const& box)
{
    using Selector = ParticlesSelector<Src::layout_mode, Src::alloc_mode>;
    Selector::template select<particle_type>(src, dst, box);
}

template<ParticleType particle_type = ParticleType::Domain, typename Src, typename Dst,
         typename box_t, typename Transformer>
void select_particles(Src const& src, Dst& dst, box_t const& box, Transformer&& transformer)
{
    using Selector = ParticlesSelector<Src::layout_mode, Src::alloc_mode>;
    Selector::template select<particle_type>(src, dst, box, transformer);
}

template<ParticleType particle_type = ParticleType::Domain, typename Src, typename box_t>
std::size_t count_particles(Src const& src, box_t const& box)
{
    using Selector = ParticlesSelector<Src::layout_mode, Src::alloc_mode>;
    return Selector::template count<particle_type>(src, box);
}

// --- implementations ---

template<auto layout_mde, auto alloc_mde>
struct ParticlesSelector
{
    static_assert(all_are<LayoutMode>(layout_mde));
    static_assert(all_are<AllocatorMode>(alloc_mde));

    auto constexpr static layout_mode = layout_mde;
    auto constexpr static alloc_mode  = alloc_mde;

    // particle_type says whether src physically lives in the domain or in a
    // (tiled) ghost layer -- tiled layouts need this to know whether to
    // intersect against a tile's domain box or its (clamped) ghost box
    template<ParticleType particle_type = ParticleType::Domain, typename SrcParticles,
             typename DstParticles, typename box_t>
    static void select(SrcParticles const&, DstParticles&, box_t const&);

    template<ParticleType particle_type = ParticleType::Domain, typename SrcParticles,
             typename DstParticles, typename box_t, typename Shift>
    static void select(SrcParticles const&, DstParticles&, box_t const&, Shift&&);

    template<ParticleType particle_type = ParticleType::Domain, typename SrcParticles,
             typename box_t>
    static std::size_t count(SrcParticles const&, box_t const&);
};


using enum LayoutMode;
using enum AllocatorMode;

// AoS, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    for (auto const& p : src)
        if (isIn(p, box))
            dst.emplace_back(p);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    for (auto const& p : src)
        if (isIn(p, box))
            dst.emplace_back(shift_particle(p, shifter));
}

// AoSPCTS, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPCTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    // ghost_box(), not the tile's bare domain box, for ghost-layer selections: a
    // tile's per-cell buckets extend into its own (clamped, uniquely-owned) ghost
    // cells, which is where patchGhost/levelGhost particles actually live -- see
    // TileSetParticles::nbr_particles_in
    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        auto& src_pc = src_tile();
        auto const tile_box
            = particle_type == ParticleType::Domain ? src_pc.box() : src_pc.ghost_box();
        if (auto const overlap = tile_box * box)
        {
            auto const lcl_src_box = src_pc.local_box(*overlap);
            for (auto it = lcl_src_box.begin(); it != lcl_src_box.end(); ++it)
                for (auto const& p : src_pc(*it))
                    dst.push_back(p);
        }
    }

    dst.template on_appended<particle_type>();
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSPCTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        auto& src_pc = src_tile();
        auto const tile_box
            = particle_type == ParticleType::Domain ? src_pc.box() : src_pc.ghost_box();
        if (auto const overlap = tile_box * box)
        {
            auto const lcl_src_box = src_pc.local_box(*overlap);
            for (auto it = lcl_src_box.begin(); it != lcl_src_box.end(); ++it)
                for (auto const& p : src_pc(*it))
                    dst.push_back(shift_particle(p, shifter));
        }
    }

    dst.template on_appended<particle_type>();
}

// AoSPC, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto const& p : sv)
            dv.emplace_back(p);
    }
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSPC, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box - fn); // BEWARE
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto p : sv)
        {
            p.iCell() = (Point{p.iCell()} + fn).toArray();
            dv.emplace_back(p);
        }
    }
}

// AoSMapped, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSMapped, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    src.export_particles(box, dst);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSMapped, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    auto const offseter = [&](auto const& particle) { return shift_particle(particle, shifter); };

    src.export_particles(box, dst, offseter);
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSMapped, CPU>::count(SrcParticles const& src, box_t const& box)
{
    return src.nbr_particles_in(box);
}

// AoS, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    for (auto const& p : src)
        if (isIn(p, box))
            ++n;
    return n;
}

// AoSPCTS, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPCTS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    for (auto& tile : const_cast<SrcParticles&>(src)())
    {
        auto& pc            = tile();
        auto const tile_box = particle_type == ParticleType::Domain ? pc.box() : pc.ghost_box();
        if (auto const overlap = tile_box * box)
        {
            auto const lcl_box = pc.local_box(*overlap);
            for (auto it = lcl_box.begin(); it != lcl_box.end(); ++it)
                n += pc(*it).size();
        }
    }
    return n;
}

// AoSPC, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPC, CPU>::count(SrcParticles const& src, box_t const& box)
{
    auto const overlap = src.ghost_box() * box; // box may extend past src's own range
    if (not overlap)
        return 0;
    auto const lcl_box = src.local_box(*overlap);
    std::size_t n      = 0;
    for (auto it = lcl_box.begin(); it != lcl_box.end(); ++it)
        n += src(*it).size();
    return n;
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP */
