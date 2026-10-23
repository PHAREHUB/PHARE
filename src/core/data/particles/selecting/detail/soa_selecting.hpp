#ifndef PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_SOA_SELECTING
#define PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_SOA_SELECTING

#include "core/utilities/memory.hpp"
#include "core/data/particles/selecting/detail/def_selecting.hpp"

#include <stdexcept>


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;

// SoA, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<SoA, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                         box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<SoA, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                         box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<SoA, CPU>::count(SrcParticles const& src, box_t const& box)
{
    throw std::runtime_error("finish this");
    return 0;
}


//

// SoA, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<SoA, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                 box_t const& box)
{
    throw std::runtime_error("finish this");
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<SoA, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                 box_t const& box, Shift&& fn)
{
    throw std::runtime_error("finish this");
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<SoA, GPU_UNIFIED>::count(SrcParticles const& src, box_t const& box)
{
    throw std::runtime_error("finish this");
    return 0;
}


//

// SoAPC, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<SoAPC, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                           box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());

    auto copy = [](auto&... args) {
        auto const& [src_cell, dst_cell] = std::forward_as_tuple(args...);
        auto dst_tuple                   = dst_cell.as_tuple();
        auto src_tuple                   = src_cell.as_tuple();
        if (dst_cell.capacity() < dst_cell.size() + src_cell.size())
            dst_cell.reserve(dst_cell.size() + src_cell.size());
        for_N<std::tuple_size_v<decltype(src_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<alloc_mode>(a.data() + dst_cell.size(), b.data(), src_cell.size());
        });
        dst_cell.resize(dst_cell.size() + src_cell.size());
    };

    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
        copy(src(*src_it), dst(*dst_it));
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<SoAPC, CPU>::select(SrcParticles const& src, DstParticles& dst,
                                           box_t const& box, Shift&& fn)
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
        for (auto const& pit : sv)
        {
            auto p    = pit.copy();
            p.iCell() = (Point{p.iCell()} + fn).toArray();
            dv.emplace_back(p);
        }
    }
}


//

// SoAPC, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t>
void ParticlesSelector<SoAPC, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                   box_t const& box)
{
    auto const lcl_src_box = src.local_box(box);
    auto const lcl_dst_box = dst.local_box(box);
    assert(lcl_src_box.shape() == lcl_dst_box.shape());

    auto copy = [](auto&... args) {
        auto const& [src_cell, dst_cell] = std::forward_as_tuple(args...);
        auto dst_tuple                   = dst_cell.as_tuple();
        auto src_tuple                   = src_cell.as_tuple();
        if (dst_cell.capacity() < dst_cell.size() + src_cell.size())
            dst_cell.reserve(dst_cell.size() + src_cell.size());
        for_N<std::tuple_size_v<decltype(src_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<alloc_mode>(a.data() + dst_cell.size(), b.data(), src_cell.size());
        });
        dst_cell.resize(dst_cell.size() + src_cell.size());
    };

    auto src_it = lcl_src_box.begin();
    auto dst_it = lcl_dst_box.begin();
    for (; src_it != lcl_src_box.end(); ++src_it, ++dst_it)
        copy(src(*src_it), dst(*dst_it));
}

template<>
template<ParticleType particle_type, typename SrcParticles, typename DstParticles, typename box_t,
         typename Shift>
void ParticlesSelector<SoAPC, GPU_UNIFIED>::select(SrcParticles const& src, DstParticles& dst,
                                                   box_t const& box, Shift&& shift)
{
    // unoptimized
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

// SoAPC, CPU

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<SoAPC, CPU>::count(SrcParticles const& src, box_t const& box)
{
    auto const lcl_box = src.local_box(box);
    std::size_t n      = 0;
    for (auto it = lcl_box.begin(); it != lcl_box.end(); ++it)
        n += src(*it).size();
    return n;
}

// SoAPC, GPU_UNIFIED

template<>
template<ParticleType particle_type, typename SrcParticles, typename box_t>
std::size_t ParticlesSelector<SoAPC, GPU_UNIFIED>::count(SrcParticles const& src, box_t const& box)
{
    auto const lcl_box = src.local_box(box);
    std::size_t n      = 0;
    for (auto it = lcl_box.begin(); it != lcl_box.end(); ++it)
        n += src(*it).size();
    return n;
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SELECTING_DETAIL_SOA_SELECTING */
