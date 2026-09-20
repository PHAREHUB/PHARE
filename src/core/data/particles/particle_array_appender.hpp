#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_APPENDER
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_APPENDER

#include "core/data/particles/particle_array_def.hpp"

#include <iterator>
#include <cstddef>

namespace PHARE::core
{

// --- public API ---

// hides the per layout/alloc impl behind operator() - see implementation section below
template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesAppender;


// Deals with one box at a time - no ParticleType, no domain/ghost distinction, so callers
// passing overlapping boxes across multiple calls will over-reserve the overlap. Resolves
// cells through the container's own local_cell/local_tile_cell (not a box*tile guess), so
// ghost-layer cells clamp-owned by a border tile (see TileSet::tag_cells_) land correctly.
template<typename Dst, typename Box_t>
void reserve(Dst& dst, Box_t const& box, std::size_t const& ppc)
{
    using enum LayoutMode;

    if constexpr (Dst::layout_mode == AoSPC)
    {
        // PerCellVector exposes indexed per-cell access directly, ghost cells included
        // (its NdArray spans the whole ghost_box, no tiling to worry about).
        for (auto const& bix : box)
            dst(dst.local_cell(*bix)).reserve(ppc);
    }
    else if constexpr (Dst::layout_mode == AoSPCTS)
    {
        // tiles are per-cell (PerCellVector) here, not flat - no tile-level reserve exists,
        // reserve each cell in the overlap directly instead. Approximate: box*tile misses
        // ghost cells clamp-owned outside a tile's own box, but reserve is just a hint.
        for (auto& tile : dst())
            if (auto const overlap = box * tile)
                for (auto const& bix : tile().local_box(*overlap))
                    tile()(bix).reserve(ppc);
    }
    else if constexpr (Dst::layout_mode == AoSCMTS)
    {
        // one flat vector per tile - same box*tile approximation as above
        for (auto& tile : dst())
            if (auto const overlap = box * tile)
                tile().reserve(tile().size() + ppc * overlap->size());
    }
    else
    {
        dst.reserve(ppc * box.size());
    }
}


// reserves box, then for every cell in box emplaces ppc particles built by fn(icell, id) -
// id is a running counter across the whole box, not reset per cell. emplace_back already
// resolves the owning cell/tile per particle (including clamp-owned ghost cells for tiled
// layouts, via the same particles_.at() lookup the span side uses) - no faster path exists
// on the vector side for tiled layouts, which expose that lookup only through their span.
template<auto type, typename Dst, typename Box_t>
void add_particles(Dst& dst, Box_t const& box, std::size_t const& ppc, auto&& fn)
{
    reserve(dst, box, ppc);

    std::size_t id = dst.size();
    for (auto const& bix : box)
        for (std::size_t i = 0; i < ppc; ++i)
            dst.emplace_back(fn(*bix, id++));

    dst.template on_appended<type>();
}


template<auto type, typename Src, typename Dst>
void append_particles(Src const& src, Dst& dst)
{
    using Appending
        = ParticlesAppender<Src::layout_mode, Src::alloc_mode, Dst::layout_mode, Dst::alloc_mode>;

    PHARE_DEBUG_DO(int const old_size = dst.size();)

    Appending{0, src.size()}.template operator()<type>(src, dst);

    // tiled level ghosts are duplicated per reaching tile, see append_to_tiles
    if constexpr (!(type == ParticleType::LevelGhost and is_tiled(Dst::layout_mode)))
        assert(dst.size() == old_size + src.size());
}


// --- implementations ---

template<auto src_layout_mde, auto src_alloc_mde, auto dst_layout_mde, auto dst_alloc_mde>
struct ParticlesAppender
{
    static_assert(all_are<LayoutMode>(src_layout_mde, dst_layout_mde));
    static_assert(all_are<AllocatorMode>(src_alloc_mde, dst_alloc_mde));

    auto constexpr static src_layout_mode = src_layout_mde;
    auto constexpr static src_alloc_mode  = src_alloc_mde;

    auto constexpr static dst_layout_mode = dst_layout_mde;
    auto constexpr static dst_alloc_mode  = dst_alloc_mde;

    template<auto type, typename Src, typename Dst>
    void operator()(Src const& src, Dst& dst);

    std::size_t const start;
    std::size_t const end;
};

using LM = LayoutMode;
using AM = AllocatorMode;


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    // reserve_flat's upfront estimate is average-case, not a hard bound - reserving here on
    // every chunk would force an exact-size reallocation instead of amortized push_back
    // growth for any shortfall.
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoSMapped, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}




// tiled dst: level ghosts are duplicated into every tile whose grown box reaches their cell
// (as the refiner does), everything else goes to the owning tile only
template<auto type, typename Src, typename Dst>
void append_to_tiles(Src const& src, Dst& dst)
{
    if constexpr (type == ParticleType::LevelGhost)
    {
        for (auto& tile : dst())
        {
            auto const tile_ghost_box = [&]() {
                if constexpr (Dst::layout_mode == LM::AoSCMTS)
                    return tile().box();
                else
                    return tile().ghost_box();
            }();
            for (auto const& p : src)
                if (isIn(p.iCell(), tile_ghost_box))
                    tile().emplace_back(p);
        }
    }
    else
        for (auto const& p : src)
            dst.emplace_back(p);

    dst.template on_appended<type>();
}

template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPCTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSPCTS, CPU>::operator()");
    append_to_tiles<type>(src, dst);
}




template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSCMTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSCMTS, CPU>::operator()");
    append_to_tiles<type>(src, dst);
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoSPC, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template on_appended<type>();
}


} // namespace PHARE::core

#endif /* PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_APPENDER */
