#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP

#include "core/utilities/box/box_spans.hpp"
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

// AoSCMTS, CPU
// Domain particles live only in their own tile, ghost-layer particles (level ghosts) only in
// their clamp-owner tile (TileSet::tag_cells_), so ghost cells are taken from that tile.
// fn(tile, box) gets each tile's owned sub-box of box, to use with the tile's cellmap
template<ParticleType particle_type, typename SrcParticles, typename Fn>
void on_owned_cell_mapped_boxes(SrcParticles const& src, auto const& box, Fn&& fn)
{
    using Box_t       = Box<int, SrcParticles::dimension>;
    auto const& tiles = src();
    for (auto const& tile : tiles)
    {
        if constexpr (particle_type == ParticleType::Domain)
        {
            if (auto const overlap = *tile * box)
                fn(tile, *overlap);
        }
        else if (auto const overlap = grow(*tile, src.ghost_cells()) * box)
            src.amr_lcl_idx(*overlap, [&](auto const& bix, auto const& lix) {
                if (tiles.at(lix) == &tile)
                    fn(tile, Box_t{bix, bix});
            });
    }
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSCMTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    on_owned_cell_mapped_boxes<particle_type>(
        src, box, [&](auto const& tile, auto const& b) { tile().export_particles(b, dst); });
    dst.template on_appended<particle_type>();
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSCMTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    auto const offseter = [&](auto const& particle) { return shift_particle(particle, shifter); };
    on_owned_cell_mapped_boxes<particle_type>(src, box, [&](auto const& tile, auto const& b) {
        tile().export_particles(b, dst, offseter);
    });
    dst.template on_appended<particle_type>();
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

// AoSCMTS, CPU
template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSCMTS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    on_owned_cell_mapped_boxes<particle_type>(
        src, box, [&](auto const& tile, auto const& b) { n += tile().nbr_particles_in(b); });
    return n;
}

} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_SELECTOR_HPP */
