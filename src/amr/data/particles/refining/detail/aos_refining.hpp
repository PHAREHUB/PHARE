// IWYU pragma: private, include "amr/data/particles/refining/particles_refining.hpp"

#ifndef PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER
#define PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER


#include "core/data/particles/particle_array_def.hpp"
#include "core/def.hpp"

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_selector.hpp"


#include "amr/utilities/box/amr_box.hpp"
#include "amr/data/particles/refining/detail/def_refining.hpp"

#include <SAMRAI/hier/BoxContainer.h>


namespace PHARE::amr
{

using enum core::LayoutMode;
using enum AllocatorMode;


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    reserve_flat<type, options>(src, dst, dst_boxes, dst_amr_box);
    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box = phare_box_from<Src::dimension>(samrai_box);
        stream_split_particles<type, options>(src, dst, dst_box, fn0, fn1);
    }
}

template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSMapped, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;

    auto const& dstBox   = dst.box();
    auto const ghost_nbr = (dst.ghost_box().shape()[0] - dst.box().shape()[0]) / 2;

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const& splitBox    = dst_box;
        auto const coarseDstBox = coarsen_box(dst_box);

        for (auto const& src_tile : src())
        {
            if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
                continue;

            if constexpr (type == core::ParticleType::Domain)
            {
                for (auto& dst_tile : dst())
                    if (auto const growbox = grow(dst_tile, ghost_nbr); splitBox * growbox)
                        if (auto const overlap = dst_box * dst_tile)
                            reserve_flat<type, options>(src_tile(), dst_tile(), *overlap,
                                                        dst_amr_box, 1.);

                // DOUBLE LOOP == LESS COSTLY LONG TERM
                for (auto& dst_tile : dst())
                    if (auto const growbox = grow(dst_tile, ghost_nbr); splitBox * growbox)
                        if (auto const overlap = dst_box * dst_tile)
                            stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap,
                                                                  fn0, fn1);
            }
            else if constexpr (type == core::ParticleType::Ghost)
            {
                for (auto& dst_tile : dst())
                {
                    // Only grow in patch-boundary directions to avoid double-counting ghost
                    // cells near tile boundaries in 2D/3D (adjacent tiles share grow overlap).
                    auto excl_box = grow(dst_tile, ghost_nbr);
                    for (std::size_t d = 0; d < Src::dimension; ++d)
                    {
                        if (dst_tile.lower[d] != dstBox.lower[d])
                            excl_box.lower[d] = dst_tile.lower[d];
                        if (dst_tile.upper[d] != dstBox.upper[d])
                            excl_box.upper[d] = dst_tile.upper[d];
                    }
                    if (*(dst.box() * excl_box) != excl_box)
                        if (auto const overlap = dst_box * excl_box)
                        {
                            reserve_flat<type, options>(src_tile(), dst_tile(), *overlap,
                                                        dst_amr_box, 1.);
                            stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap,
                                                                  fn0, fn1);
                        }
                }
            }
            else
            {
                assert(false);
            }
        }
    }
}



template<auto type, typename Dst, typename Tile, typename Box_t>
std::optional<Box_t> refine_tile_overlap(Dst const& dst, Tile const& dst_tile, Box_t const& dst_box,
                                         std::size_t const ghost_nbr)
{
    // dst_box can extend past this array's own domain (periodic wrap, a restricted/
    // temporary ParticlesData SAMRAI builds internally, or the level ghost layer).
    // Particles are never duplicated across tiles, so each bound's margin is owned
    // by whichever tile touches the domain edge it's adjacent to (same clamp rule
    // tag_cells_ uses, so TileSet::at(ghost_cell) finds them); a tile that doesn't
    // own a bound just clips to its own edge there instead of the whole region
    // being rejected. Readers needing a neighbour's level ghosts gather via at().
    auto const growbox = grow(dst_tile, ghost_nbr);
    if (not(dst_box * growbox))
        return std::nullopt;

    auto const& domain_box = dst.box();
    auto const& tile_box   = *dst_tile;
    auto region            = dst_box;
    for (std::size_t d = 0; d < Dst::dimension; ++d)
    {
        if (dst_box.lower[d] < domain_box.lower[d])
        {
            bool const owns_edge = domain_box.lower[d] >= tile_box.lower[d]
                                   and domain_box.lower[d] <= tile_box.upper[d];
            region.lower[d]
                = owns_edge ? dst_box.lower[d] : std::max(dst_box.lower[d], tile_box.lower[d]);
        }
        else
            region.lower[d] = std::max(dst_box.lower[d], tile_box.lower[d]);

        if (dst_box.upper[d] > domain_box.upper[d])
        {
            bool const owns_edge = domain_box.upper[d] >= tile_box.lower[d]
                                   and domain_box.upper[d] <= tile_box.upper[d];
            region.upper[d]
                = owns_edge ? dst_box.upper[d] : std::min(dst_box.upper[d], tile_box.upper[d]);
        }
        else
            region.upper[d] = std::min(dst_box.upper[d], tile_box.upper[d]);

        if (region.lower[d] > region.upper[d])
            return std::nullopt;
    }
    return region;
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSPCTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    // per-tile destination is per-cell storage, with no single flat buffer for reserve()
    // to size - see the shrink_to_fit_in comment below.
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    auto const ghost_nbr                      = options.particle_ghost_width;

    // Per-cell reserve: each dst cell's own coarse-parent cell count (not a box-wide average)
    // - uneven patterns (3D's 6-particle Pink) still spread unevenly across a coarse cell's
    // children, but that's a per-cell rounding slop, not the systematic box-total error a
    // uniform box-average would introduce.
    static_assert(type == core::ParticleType::Domain || type == core::ParticleType::Ghost);

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const coarseDstBox = coarsen_box(dst_box);

        auto const dst_tile_overlap = [&](auto& dst_tile) {
            return refine_tile_overlap<type>(dst, dst_tile, dst_box, ghost_nbr);
        };

        // reserve pass: one .reserve() call per dst cell, via reserve_flat_count against the
        // whole (tiled) src - count_particles<AoSPCTS> sums each tile's own domain box only
        // (never ghost cells), so cells sitting in more than one neighbor tile's ghost halo
        // don't get summed twice the way checking every src_tile's ghost_box used to.
        for (auto& dst_tile : dst())
            if (auto const overlap = dst_tile_overlap(dst_tile))
            {
                auto& dst_pc = dst_tile();
                for (auto const& bix : *overlap)
                {
                    core::Box<int, Src::dimension> const cell_box{*bix, *bix};
                    dst_pc(dst_pc.local_cell(*bix))
                        .reserve(reserve_flat_count<type, options>(
                            src, cell_box, dst_amr_box,
                            per_cell_density_factor<options, type, Src>()));
                }
            }

        for (auto const& src_tile : src())
        {
            if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
                continue;

            for (auto& dst_tile : dst())
                if (auto const overlap = dst_tile_overlap(dst_tile))
                    stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap, fn0,
                                                          fn1);
        }
    }

    dst.template on_appended<type>();
}


template<>
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSCMTS, CPU>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    // one flat cell-mapped array per tile: same per-tile ownership as AoSPCTS, but each dst
    // tile has a single flat buffer, so reserve per tile overlap instead of per cell
    auto&& [dst, src, dst_boxes, dst_amr_box] = args;
    auto const ghost_nbr                      = options.particle_ghost_width;

    static_assert(type == core::ParticleType::Domain || type == core::ParticleType::Ghost);

    for (auto const& samrai_box : dst_boxes)
    {
        auto const dst_box      = phare_box_from<Src::dimension>(samrai_box);
        auto const coarseDstBox = coarsen_box(dst_box);

        for (auto& dst_tile : dst())
            if (auto const overlap = refine_tile_overlap<type>(dst, dst_tile, dst_box, ghost_nbr))
                reserve_flat<type, options>(src, dst_tile(), *overlap, dst_amr_box,
                                            per_cell_density_factor<options, type, Src>());

        for (auto const& src_tile : src())
        {
            if (not(coarseDstBox * grow(src_tile, ghost_nbr)))
                continue;

            for (auto& dst_tile : dst())
                if (auto const overlap
                    = refine_tile_overlap<type>(dst, dst_tile, dst_box, ghost_nbr))
                    stream_split_particles<type, options>(src_tile(), dst_tile(), *overlap, fn0,
                                                          fn1);
        }
    }

    dst.template on_appended<type>();
}



template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0, auto fn1)
{
    ParticlesRefiner<AoS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSMapped, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                          auto fn1)
{
    ParticlesRefiner<AoSMapped, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSTS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                      auto fn1)
{
    ParticlesRefiner<AoSTS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}


template<> // slow
template<auto options, auto type, typename Src, typename Dst>
void ParticlesRefiner<AoSPCTS, GPU_UNIFIED>::operator()(RefineArgs<Src, Dst>& args, auto fn0,
                                                        auto fn1)
{
    ParticlesRefiner<AoSPCTS, CPU>{}.template operator()<options, type>(args, fn0, fn1);
}




} // namespace PHARE::amr


#endif /* PHARE_AMR_DATA_PARTICLES_REFINING_DETAIL_AOS_REFINER */
