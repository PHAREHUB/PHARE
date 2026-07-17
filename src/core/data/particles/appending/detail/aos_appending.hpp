// IWYU pragma: private, include "core/data/particles/appending/particles_appending.hpp"
#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_AOS_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_AOS_APPENDING


#include "core/utilities/memory.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"


#include <iterator>

namespace PHARE::core
{

using LM = LayoutMode;
using AM = AllocatorMode;


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    dst.reserve(dst.size() + src.size());
    std::copy(src.begin(), src.end(), std::back_inserter(dst));
}




template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSTS, CPU>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template sync<2, type>();
}

template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSCMTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSCMTS, CPU>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template sync<2, type>();
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPCTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSPCTS, CPU>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template sync<2, type>();
}




template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoS, AM::CPU, LM::AoSPC, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    for (auto const& p : src)
        dst.emplace_back(p);

    dst.template sync<2, type>();
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSPC, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        auto& dst_arr        = dst(lix);
        auto const curr_size = dst_arr.size();
        Dst::resize(dst_arr, curr_size + size);
        gpu::copy(dst_arr.data() + curr_size, src.data() + src_start, size);
        src_start += size;
    };

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::AoSTS, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    PHARE_LOG_SCOPE(3, "ParticlesAppender<AoSMapped, CPU, AoSTS, GPU_UNIFIED>::operator()");

    for (auto const& p : src)
        dst.emplace_back(p);

    // auto const overlap = src.box() * dst.ghost_box();
    // if (!overlap)
    //     return;

    // std::int32_t src_start = 0;
    // auto finish            = [&](auto const lix, auto const size) {
    //     auto& particles = dst(lix);
    //     mem::copy<GPU_UNIFIED>(particles.data() + particles.size(), src.data() + src_start,
    //     size); particles.resize(particles.size() + size); src_start += size;
    // };

    // for (auto& tile : dst())
    //     tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
    //                        return src.nbr_particles_in(bix);
    //                    }));

    // dst.reset_views();

    // for (auto const bix : *overlap)
    //     if (auto const size = src.nbr_particles_in(bix))
    //         finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}



} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_AOS_APPENDING */
