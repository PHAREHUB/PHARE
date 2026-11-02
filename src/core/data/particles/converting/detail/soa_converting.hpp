// IWYU pragma: private, include "core/data/particles/converting/particles_converting.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_SOA_CONVERTING
#define PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_SOA_CONVERTING

#include "core/utilities/memory.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"
#include "core/data/particles/converting/detail/def_converting.hpp"

namespace PHARE::core
{
using enum LayoutMode;
using enum AllocatorMode;

template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<SoA, CPU, AoS, CPU>::operator()( //
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());
    for (std::size_t i = 0; i < src.size(); ++i)
        out.emplace_back(src.copy(i));
    assert(src.size() == out.size());
    return out;
}

template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<SoAVX, CPU, AoS, CPU>::operator()( //
    Src const& src, GridLayout const& layout)
{
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());
    for (std::size_t i = 0; i < src.size(); ++i)
        out.emplace_back(src.copy(i));
    assert(src.size() == out.size());
    return out;
}


template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<SoATS, CPU, AoS, CPU>::operator()( //
    Src const& src, GridLayout const& layout)
{
    static constexpr std::uint8_t N = 128;

    auto copy = [](auto&... args) {
        auto const& [tmp, particles, i, S] = std::forward_as_tuple(args...);
        auto tmp_tuple                     = tmp.as_tuple();
        auto src_tuple                     = particles.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(tmp_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<CPU>(a.data(), b.data() + N * i, S);
        });
    };

    SoAArray<GridLayout::dimension, N> tmp;
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
    {
        auto const& particles = tile();
        auto const full_count = particles.size() / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
        {
            copy(tmp, particles, i, N);
            out.append(tmp, 0, N);
        }
        if (auto const remaining = particles.size() % N)
        {
            copy(tmp, particles, i, remaining);
            out.append(tmp, 0, remaining);
        }
    }
    assert(src.size() == out.size());
    return out;
}


template<>
template<typename Dst, typename Src, typename GridLayout>
Dst ParticlesConverter<SoATS, GPU_UNIFIED, AoS, CPU>::operator()( //
    Src const& src, GridLayout const& layout)
{
    static constexpr std::uint8_t N = 128;

    auto copy = [](auto&... args) {
        auto const& [tmp, particles, i, S] = std::forward_as_tuple(args...);
        auto tmp_tuple                     = tmp.as_tuple();
        auto src_tuple                     = particles.as_tuple();
        for_N<std::tuple_size_v<decltype(tmp_tuple)>>([&](auto vi) {
            auto& a = std::get<vi>(tmp_tuple);
            auto& b = std::get<vi>(src_tuple);
            mem::copy<GPU_UNIFIED>(a.data(), b.data() + N * i, S);
        });
    };

    SoAArray<GridLayout::dimension, N> tmp;
    auto out = make_particles<Dst>(layout);
    out.reserve(src.size());

    for (auto const& tile : src())
    {
        auto const& particles = tile();
        auto const full_count = particles.size() / N;
        std::size_t i         = 0;
        for (; i < full_count; ++i)
        {
            copy(tmp, particles, i, N);
            out.append(tmp, 0, N);
        }
        if (auto const remaining = particles.size() % N)
        {
            copy(tmp, particles, i, remaining);
            out.append(tmp, 0, remaining);
        }
    }
    assert(src.size() == out.size());
    return out;
}

} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_CONVERTING_DETAIL_SOA_CONVERTING */
