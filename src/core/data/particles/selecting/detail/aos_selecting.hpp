
#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING


// #include "core/utilities/memory.hpp"
#include "core/logger.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_partitioner.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"


#include <stdexcept>

namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;


//

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

//

// AoSTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSTS, CPU>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    // Domain particles live in a tile's own (domain) box; patch/level ghost particles
    // live only in a tile's own (clamped, uniquely-owned) ghost cells -- see
    // TileSetParticles::nbr_particles_in for the same grow()-based reasoning.
    auto const tile_box = [&](auto& tile) {
        if constexpr (particle_type == ParticleType::Domain)
            return *tile;
        else
        {
            auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;
            return grow(tile, ghost_cells);
        }
    };

    for (auto& src_tile : const_cast<SrcParticles&>(src)())
        if (auto const src_box_overlap = tile_box(src_tile) * box)
            for (auto const& p : src_tile())
                if (isIn(p, *src_box_overlap))
                    dst.push_back(p);

    if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
        dst.template sync<2>();
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoSTS, CPU>::select( // box is unshifted global AMR indexing
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shifter)
{
    using enum LayoutMode;
    static_assert(any_in(DstParticles::layout_mode, AoS, AoSTS), "Otherwise not implemented!");

    auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;

    for (auto& src_tile : const_cast<SrcParticles&>(src)())
    {
        if (auto const src_overlap = grow(src_tile, ghost_cells) * box)
        {
            if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
            {
                for (auto const& p : src_tile())
                    if (isIn(p, *src_overlap))
                        dst.push_back(shift_particle(p, shifter));
            }
            else
            {
                auto const copy = [&](auto& s, auto& d, auto const& b) {
                    auto const in_dst = partition_particles(s, b);
                    if (in_dst.size() == 0)
                        return;

                    auto const old_dst_size = d.size();
                    d.resize(d.size() + in_dst.size());

                    mem::copy<CPU>(d.data() + old_dst_size, s.data(), in_dst.size());
                    if (any(shifter, [](auto const e) { return e != 0; }))
                        for (std::size_t i = old_dst_size; i < d.size(); i++)
                            d[i].iCell() = *((shifter) + d[i].iCell());
                };
                copy(src_tile(), dst, *src_overlap);
            }
        }
    }

    if constexpr (DstParticles::layout_mode == LayoutMode::AoSTS)
        dst.template sync<2>();
}


//

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

    if constexpr (any_in(DstParticles::layout_mode, AoSTS, AoSPCTS))
        dst.template sync<2>();
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

    if constexpr (any_in(DstParticles::layout_mode, AoSTS, AoSPCTS))
        dst.template sync<2>();
}


//

// AoS, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<AoS, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

//

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

//

// AoSPC, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select( //
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
void ParticlesSelector<AoSPC, GPU_UNIFIED>::select( //
    SrcParticles const& src, DstParticles& dst, box_t const& box, Shift&& shift)
{ // unoptimized
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box - shift);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());
    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
    {
        auto& sv = src(*src_it);
        auto& dv = dst(*dst_it);
        dv.reserve(dv.size() + sv.size());
        for (auto const& pit : sv)
        {
            auto p    = pit.copy();
            p.iCell() = (Point{p.iCell()} + shift).toArray();
            dv.emplace_back(p);
        }
    }
}

//

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

//

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

// AoS, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoS, GPU_UNIFIED>::count(SrcParticles const& src, box_t const& box)
{
    throw std::runtime_error("finish this");
    return 0;
}

// AoSTS, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSTS, CPU>::count(SrcParticles const& src, box_t const& box)
{
    std::size_t n = 0;
    for (auto& tile : const_cast<SrcParticles&>(src)())
    {
        auto const box_for_overlap = [&]() {
            if constexpr (particle_type == ParticleType::Domain)
                return *tile;
            else
            {
                auto const ghost_cells = (src.ghost_box().shape()[0] - src.box().shape()[0]) / 2;
                return grow(tile, ghost_cells);
            }
        }();
        if (auto const overlap = box_for_overlap * box)
            for (auto const& p : tile())
                if (isIn(p, *overlap))
                    ++n;
    }
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
        auto& pc = tile();
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
    auto const lcl_box = src.local_box(box);
    std::size_t n      = 0;
    for (auto it = lcl_box.begin(); it != lcl_box.end(); ++it)
        n += src(*it).size();
    return n;
}

// AoSPC, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<AoSPC, GPU_UNIFIED>::count(SrcParticles const& src, box_t const& box)
{
    auto const lcl_box = src.local_box(box);
    std::size_t n      = 0;
    for (auto it = lcl_box.begin(); it != lcl_box.end(); ++it)
        n += src(*it).size();
    return n;
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_AOS_SELECTING */
