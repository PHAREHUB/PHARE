
#ifndef PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER
#define PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER


#include "core/def.hpp"
#include "core/def/phare_config.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_appender.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/exporting/detail/def_exporting.hpp"


#include <stdexcept>


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


template<>
template<typename Src, typename Dst, typename Box_t>
void ParticlesExporter<AoSTS, CPU>::move_particles(Src& src, Dst& dst, Box_t const& box,
                                                   std::size_t const growby)
{
    throw std::runtime_error("do not use");

    static_assert(Src::layout_mode == Dst::layout_mode);
    assert(src().size() == dst().size());

    auto const old_size = src.size();
    for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
    {
        auto& src_tile     = src()[tidx];
        auto const src_box = grow(src_tile, growby);
        if (!src_tile().size() or !(box * src_box))
            continue;

        auto& dst_tile        = dst()[tidx];
        auto const not_in_box = partition_particles_not_in(src_tile(), box);
        auto const start      = not_in_box.size();
        auto const end        = src_tile().size();
        auto const size       = end - start;

        PHARE_DEBUG_DO({
            for (std::size_t i = 0; i < start; ++i)
            {
                assert(not isIn(src_tile()[i], box));
                assert(not isIn(src_tile()[i], box));
            }
            for (std::size_t i = start; i < end; ++i)
            {
                assert(isIn(src_tile()[i], box));
            }
            for (std::size_t i = 0; i < dst_tile().size(); ++i)
            {
                assert(isIn(dst_tile()[i], dst_tile));
            }
        })

        if (!size)
            continue;

        append_particles<ParticleType::Domain>(src_tile().view(start, size), dst_tile());

        PHARE_DEBUG_DO({
            for (std::size_t i = 0; i < start; ++i)
            {
                assert(not isIn(src_tile()[i], box));
            }
            for (std::size_t i = start; i < end; ++i)
            {
                assert(isIn(src_tile()[i], box));
            }
            for (std::size_t i = 0; i < dst_tile().size(); ++i)
            {
                assert(isIn(dst_tile()[i], dst_tile));
            }
        })


        if (old_size)
        {
            auto const& p = src_tile()[0];
            assert(start);
        }

        src_tile().resize(start);
    }

    src.sync();
    dst.sync();

    auto const new_size = src.size();
    if (old_size)
    {
        assert(new_size);
    }
}


template<> // slow
template<typename Src, typename Dst, typename Box_t>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_particles(Src& src, Dst& dst, Box_t const& box,
                                                           std::size_t const growby)
{
    ParticlesExporter<AoSTS, CPU>{}.move_particles(src, dst, box, growby);
}


template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::move_in_domain(Dst& dst, Src& src, Box<int, dim> const& domain_box)
{
    for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
    {
        auto& src_tile = src()[tidx];
        auto& dst_tile = dst()[tidx];

        auto const tile_ghost_box = grow(src_tile, 1);
        assert(src_tile == dst_tile);

        if (tile_ghost_box * domain_box == tile_ghost_box) // not border tile
            continue;

        auto const end = src_tile().size();

        for (std::size_t i = end; i-- > 0;)
        {
            auto const& p = src_tile()[i];
            if (isIn(p, domain_box))
            {
                dst.push_back(p);
                src_tile().vector().erase(src_tile().vector().begin() + i);
            }
        }
    }

    src.sync();
    dst.sync();
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                        Box<int, dim> const& domain_box,
                                                        Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);

    for (std::size_t tidx = 0; tidx < src().size(); ++tidx)
    {
        auto& src_tile = src()[tidx];
        auto& dst_tile = dst()[tidx];

        auto const tile_ghost_box = grow(src_tile, 1);
        assert(src_tile == dst_tile);

        if (tile_ghost_box * domain_box == tile_ghost_box) // not border tile
            continue;

        auto const end = src_tile().size();
        for (std::size_t i = end; i-- > 0;)
        {
            auto const& p = src_tile()[i];
            for (auto const& ghost_layer : ghost_layer_boxes)
            {
                if (isIn(p, ghost_layer))
                {
                    dst_tile().push_back(p);
                    src_tile().vector().erase(src_tile().vector().begin() + i);
                    break;
                }
            }
        }
    }

    src.sync();
    dst.sync();
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                        Box<int, dim> const& domain_box,
                                                        Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


template<> // slow
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_in_domain(Dst& dst, Src& src,
                                                           Box<int, dim> const& domain_box)
{
    ParticlesExporter<AoSTS, CPU>{}.move_in_domain(dst, src, domain_box);
}

template<> // slow
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_in_ghost_layer(Dst& dst, Src& src,
                                                                Box<int, dim> const& domain_box,
                                                                Box<int, dim> const& ghost_box)
{
    ParticlesExporter<AoSTS, CPU>{}.move_in_ghost_layer(dst, src, domain_box, ghost_box);
}

template<> // slow
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::move_in_ghost_layer(Dst& dst, Src& src,
                                                                Box<int, dim> const& domain_box,
                                                                Boxes const& ghost_boxes)
{
    ParticlesExporter<AoSTS, CPU>{}.move_in_ghost_layer(dst, src, domain_box, ghost_boxes);
}


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    src.erase(box);
    src.sortMapping();
}
template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    for (auto const& box : boxes)
        src.erase(box);
}

template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSMapped, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    throw std::runtime_error("todo");
}
template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSMapped, GPU_UNIFIED>::delete_particles_not_in(Src& src,
                                                                        Boxes const& boxes)
{
    throw std::runtime_error("todo");
}

template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    auto const ghost_nbr = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& tile : src())
    {
        if (tile().size() == 0)
            continue;

        if (auto const tile_gbox = grow(*tile, ghost_nbr); tile_gbox * box)
        {
            auto const range = partition_particles(tile(), box);
            auto const resiz = range.size();

            assert(resiz <= tile().size());
            tile().resize(resiz);
        }
    }

    src.sync();
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto const ghost_nbr = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& tile : src())
    {
        if (tile().size() == 0)
            continue;

        if (auto const tile_gbox = grow(*tile, ghost_nbr); tile_gbox * src.box() != tile_gbox)
        {                                                           // only patch border tiles
            auto const ranges = partition_particles(tile(), boxes); // what
            auto const resiz  = sum_from(ranges, [](auto const& r) { return r.size(); });

            assert(resiz <= tile().size());
            tile().resize(resiz);
        }
    }

    src.sync();
}



template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    for (auto& tile : src())
    {
        if (tile().size() == 0)
            continue;

        auto const isect = box * *tile;

        if (!isect)
        {
            tile().clear();
        }
        else if (*isect != *tile)
        {
            auto const range = makeRange(tile(), 0, tile().size());
            auto const kept  = tile().partition(
                std::move(range), [&](auto const& cell) { return isIn(cell, *isect); });
            tile().resize(kept.size());
            tile().remap();
        }
    }

    src.sync();
}


template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSCMTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    using box_t = Box<int, Src::dimension>;

    for (auto& tile : src())
    {
        if (tile().size() == 0)
            continue;

        std::vector<box_t> static isects;
        isects.clear();
        for (auto const& box : boxes)
            if (auto const isect = box * *tile)
                isects.push_back(*isect);

        if (isects.empty())
        {
            tile().clear();
        }
        else
        {
            auto const range = makeRange(tile(), 0, tile().size());
            auto const kept  = tile().partition(
                std::move(range), [&](auto const& cell) { return isIn(cell, isects); });
            tile().resize(kept.size());
            tile().remap();
        }
    }

    src.sync();
}


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles_not_in(Src& src,
                                                                    Box<int, dim> const& box)
{
    ParticlesExporter<AoSTS, CPU>{}.delete_particles_not_in(src, box);
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSTS, GPU_UNIFIED>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    ParticlesExporter<AoSTS, CPU>{}.delete_particles_not_in(src, boxes);
}


template<>
template<typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, CPU>::delete_particles_not_in(Src& src, Box<int, dim> const& box)
{
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for (auto const& bix : cps.local_box(gb))
        {
            auto& cell_parts = cps(bix);
            if (cell_parts.size() == 0)
                continue;
            std::array<int, dim> gcell;
            for (std::size_t d = 0; d < dim; ++d)
                gcell[d] = gb.lower[d] + static_cast<int>(bix[d]);
            if (!isIn(gcell, box))
                cell_parts.clear();
        }
        cps.template sync<2, ParticleType::Ghost>();
    }
    src.sync();
}

template<>
template<typename Src, typename Boxes>
void ParticlesExporter<AoSPCTS, CPU>::delete_particles_not_in(Src& src, Boxes const& boxes)
{
    auto constexpr dim = Src::dimension;
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for (auto const& bix : cps.local_box(gb))
        {
            auto& cell_parts = cps(bix);
            if (cell_parts.size() == 0)
                continue;
            std::array<int, dim> gcell;
            for (std::size_t d = 0; d < dim; ++d)
                gcell[d] = gb.lower[d] + static_cast<int>(bix[d]);
            bool in_any = false;
            for (auto const& b : boxes)
                if (isIn(gcell, b)) { in_any = true; break; }
            if (!in_any)
                cell_parts.clear();
        }
        cps.template sync<2, ParticleType::Ghost>();
    }
    src.sync();
}


// AoSPCTS domain <-> patch-ghost/level-ghost particle exchange: per-cell aware, mirrors
// delete_particles_not_in<AoSPCTS, CPU> traversal. `dst.push_back()` resolves the
// receiving tile+cell from the particle's own iCell(), so src/dst tiles need not align.
template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, CPU>::move_in_domain(Dst& dst, Src& src,
                                                     Box<int, dim> const& domain_box)
{
    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for (auto const& bix : cps.local_box(gb))
        {
            auto& cell_parts = cps(bix);
            for (std::size_t i = cell_parts.size(); i-- > 0;)
            {
                if (isIn(cell_parts.iCell(i), domain_box))
                {
                    dst.push_back(cell_parts[i]);
                    cell_parts.assign(cell_parts.size() - 1, i);
                    cell_parts.pop_back();
                }
            }
        }
        cps.template sync<2, ParticleType::Ghost>();
    }
    src.sync();
    dst.sync();
}

template<>
template<typename Dst, typename Src, std::size_t dim>
void ParticlesExporter<AoSPCTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Box<int, dim> const& ghost_box)
{
    auto const ghost_layer_boxes = ghost_box.remove(domain_box);

    for (auto& tile : src())
    {
        auto& cps      = tile();
        auto const& gb = cps.ghost_box();
        for (auto const& bix : cps.local_box(gb))
        {
            auto& cell_parts = cps(bix);
            for (std::size_t i = cell_parts.size(); i-- > 0;)
            {
                auto const& icell = cell_parts.iCell(i);
                for (auto const& ghost_layer : ghost_layer_boxes)
                {
                    if (isIn(icell, ghost_layer))
                    {
                        dst.push_back(cell_parts[i]);
                        cell_parts.assign(cell_parts.size() - 1, i);
                        cell_parts.pop_back();
                        break;
                    }
                }
            }
        }
        cps.template sync<2, ParticleType::Ghost>();
    }
    src.sync();
    dst.sync();
}

template<>
template<typename Dst, typename Src, std::size_t dim, typename Boxes>
void ParticlesExporter<AoSPCTS, CPU>::move_in_ghost_layer(Dst& dst, Src& src,
                                                          Box<int, dim> const& domain_box,
                                                          Boxes const& ghost_boxes)
{
    for (auto const& gb : ghost_boxes)
        this->move_in_ghost_layer(dst, src, domain_box, gb);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_EXPORTING_DETAIL_AOS_EXPORTER */
