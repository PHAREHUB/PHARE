// IWYU pragma: private, include "core/data/particles/appending/particles_appending.hpp"
#ifndef PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOA_APPENDING
#define PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOA_APPENDING

#include "core/def/phare_config.hpp"
#include "core/utilities/memory.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/arrays/particle_array_soavx.hpp"

// #include "core/data/particles/particle_array.hpp"
// #include "core/data/particles/particle_packer.hpp"
// #include "core/data/particles/particle_array_converter.hpp"
#include "core/data/particles/appending/detail/def_appending.hpp"


namespace PHARE::core
{

using LM = LayoutMode;
using AM = AllocatorMode;


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::SoA, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    dst.reserve(dst.size() + src.size());
    for (auto const& p : src)
        dst.push_back(p);
}


// AoSMapped-CPU to SoAPC-GPU_UNIFIED
template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::SoAPC, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    static constexpr std::uint8_t N = 128;

    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    SoAVXArray<Src::dimension, N> tmp;

    auto copy = [](auto&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        auto dst_tuple                     = particles.as_tuple();
        auto tmp_tuple                     = tmp.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(tmp_tuple);
            mem::copy<dst_alloc_mode>(a.data() + particles.size(), b.data(), S);
        });
        particles.resize(particles.size() + S);
    };

    std::int32_t src_start = 0;
    auto do_copy           = [&](auto&&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        for (std::size_t pidx = 0; pidx < S; ++pidx)
            tmp.assign(src[pidx + src_start], pidx);
        copy(particles, tmp, i, S);
        src_start += S;
    };


    std::size_t tmp_start = 0;
    auto const tmp_tuple  = tmp.as_tuple();
    auto finish           = [&](auto const lix, auto const size) {
        auto& particles = dst(lix);
        particles.reserve(particles.size() + size);
        auto const full_count = size / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
            do_copy(particles, tmp, i, N);
        if (auto const remaining = size % N)
            do_copy(particles, tmp, i, remaining);

        particles.resize(particles.size() + size);
    };

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}



template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::SoATS, AM::GPU_UNIFIED>::operator()( //
    Src const& src, Dst& dst)
{
    static constexpr std::uint8_t N = 128;

    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    auto copy = [](auto&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        auto dst_tuple                     = particles.as_tuple();
        auto tmp_tuple                     = tmp.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(tmp_tuple);
            mem::copy<dst_alloc_mode>(a.data() + particles.size(), b.data(), S);
        });
        particles.resize(particles.size() + S);
    };

    std::int32_t src_start = 0;
    auto do_copy           = [&](auto&&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        for (std::size_t pidx = 0; pidx < S; ++pidx)
            tmp.assign(src[pidx + src_start], pidx);
        copy(particles, tmp, i, S);
        src_start += S;
    };

    SoAArray<Src::dimension, N> tmp;

    auto finish = [&](auto const lix, auto const size) {
        assert(dst().at(lix));
        auto& tile            = *dst().at(lix);
        auto& particles       = tile();
        auto const full_count = size / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
            do_copy(particles, tmp, i, N);
        if (auto const remaining = size % N)
            do_copy(particles, tmp, i, remaining);
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    dst.reset_views();

    for (auto const bix : *overlap)
        if (auto const size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}

template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::SoATS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    std::int32_t src_start = 0;
    auto finish            = [&](auto const lix, auto const size) {
        dst(lix).append(src, src_start, size);
        src_start += size;
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}



template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::SoAVX, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    dst.reserve(dst.size() + src.size());
    for (auto const& p : src)
        dst.push_back(p);
}


template<>
template<auto type, typename Src, typename Dst>
void ParticlesAppender<LM::AoSMapped, AM::CPU, LM::SoAVXTS, AM::CPU>::operator()( //
    Src const& src, Dst& dst)
{
    static constexpr std::uint8_t N = 128;

    auto const overlap = src.box() * dst.ghost_box();
    if (!overlap)
        return;

    auto copy = [](auto&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        auto dst_tuple                     = particles.as_tuple();
        auto tmp_tuple                     = tmp.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(dst_tuple);
            auto& b = std::get<vi>(tmp_tuple);
            mem::copy<dst_alloc_mode>(a.data() + particles.size(), b.data(), S);
        });
        particles.resize(particles.size() + S);
    };

    std::int32_t src_start = 0;

    auto do_copy = [&](auto&&... args) {
        auto const& [particles, tmp, i, S] = std::forward_as_tuple(args...);
        for (std::size_t pidx = 0; pidx < S; ++pidx)
            tmp.assign(src[pidx + src_start], pidx);
        copy(particles, tmp, i, S);
        src_start += S;
    };

    SoAVXArray<Src::dimension, N> tmp;

    auto finish = [&](auto const lix, auto const size) {
        auto& particles       = dst(lix);
        auto const full_count = size / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
            do_copy(particles, tmp, i, N);
        if (auto const remaining = size % N)
            do_copy(particles, tmp, i, remaining);
    };

    for (auto& tile : dst())
        tile().reserve(tile().size() + sum_from(*tile, [&](auto const bix) {
                           return src.nbr_particles_in(bix);
                       }));

    for (auto const bix : *overlap)
        if (auto size = src.nbr_particles_in(bix))
            finish(dst.local_cell(bix), size);

    dst.template sync<2, type>();
}



} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_APPENDING_DETAIL_SOA_APPENDING */
