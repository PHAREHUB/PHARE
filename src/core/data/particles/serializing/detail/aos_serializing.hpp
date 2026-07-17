// IWYU pragma: private, include "core/data/particles/serializing/particles_serializing.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_AOS_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_AOS_SERIALIZING

#include "core/utilities/span.hpp"
#include "core/data/particles/serializing/detail/def_serializing.hpp"


namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;

// AoS
template<>
template<typename Src>
void ParticlesSerializer<AoS, CPU>::operator()(std::string const& file_name, Src const& src)
{
    using Particle_t = typename Src::value_type;

    open_file_to_write(file_name).write(src.vector().data(), src.size());
}


template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<AoS, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    static_assert(Dst::alloc_mode == Src::alloc_mode);
    static_assert(layout_mode == Src::layout_mode || any_in(Src::layout_mode, AoS, AoSMapped));

    using Particle_t = typename Dst::value_type;

    auto file            = open_file_from_start(file_name);
    auto const curr_size = dst.size();
    dst.resize(curr_size + (file.size() / sizeof(Particle_t)));

    file.read<Particle_t>(dst.vector().data() + curr_size, dst.size() - curr_size);
}


template<>
template<typename Src, std::uint16_t N, typename Fn>
void ParticlesDeserializer<AoS, CPU>::readN(std::string const& file_name, Fn fn)
{
    using Particle_t  = typename Src::value_type;
    using Particles_t = std::array<Particle_t, N>;

    auto file           = open_file_from_start(file_name);
    auto const new_size = file.size() / sizeof(Particle_t);
    Particles_t particles;

    std::size_t i = 0;
    for (; i < new_size; i += N)
        fn(make_span(file.read(particles.data(), N), N));

    auto const last = new_size % N;
    for (; i < new_size; ++i)
        fn(make_span(file.read(particles.data(), last), last));
}


// AoSMapped
template<>
template<typename Src>
void ParticlesSerializer<AoSMapped, CPU>::operator()(std::string const& file_name, Src const& src)
{
    ParticlesSerializer<AoS, CPU>{}(file_name, src);
}


template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<AoSMapped, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    ParticlesDeserializer<AoS, CPU>{}(file_name, dst);
    dst.remap();
}



template<>
template<typename Src, std::uint16_t N, typename Fn>
void ParticlesDeserializer<AoSMapped, CPU>::readN(std::string const& file_name, Fn fn)
{
    ParticlesDeserializer<AoS, CPU>{}.template readN<Src, N>(file_name, fn);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_AOS_SERIALIZING */
