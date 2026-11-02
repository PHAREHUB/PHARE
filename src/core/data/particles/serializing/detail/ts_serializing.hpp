// IWYU pragma: private, include "core/data/particles/serializing/particles_serializing.hpp"

#ifndef PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_TS_SERIALIZING
#define PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_TS_SERIALIZING

#include "core/data/particles/particle_array.hpp"
#include "core/data/particles/serializing/detail/def_serializing.hpp"

namespace PHARE::core
{

using enum LayoutMode;
using enum AllocatorMode;

// AoSTS / AoSPCTS - particles are spread across tiles (and, for AoSPCTS, per-cell
// buckets within each tile), so there is no single contiguous buffer to write/read in
// one shot as the flat AoS layouts do. Write/read one particle at a time via
// per_particle/emplace_back rather than staging through an intermediate flat buffer.
// value_type is the same Particle_t regardless of layout_mode (see
// ParticleArray::value_type), so the file format is byte-compatible with the flat
// AoS/AoSMapped serializers.
template<auto layout_mode>
struct TiledParticlesSerializing
{
    template<typename Src>
    static void write(std::string const& file_name, Src const& src, auto& self)
    {
        auto file = self.open_file_to_write(file_name);
        per_particle(src, [&](auto const& p) { file.write(&p, 1); });
    }

    template<typename Dst, typename Src>
    static void read(std::string const& file_name, Dst& dst, auto& self)
    {
        static_assert(Dst::alloc_mode == Src::alloc_mode);
        using Particle_t = typename Dst::value_type;

        auto file             = self.open_file_from_start(file_name);
        auto const nbr_parts  = file.size() / sizeof(Particle_t);
        Particle_t p;
        for (std::size_t i = 0; i < nbr_parts; ++i)
            dst.emplace_back(*file.read(&p, 1));

        dst.template on_appended<ParticleType::Domain>();
    }
};


template<>
template<typename Src>
void ParticlesSerializer<AoSTS, CPU>::operator()(std::string const& file_name, Src const& src)
{
    TiledParticlesSerializing<AoSTS>::write(file_name, src, *this);
}

template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<AoSTS, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    TiledParticlesSerializing<AoSTS>::template read<Dst, Src>(file_name, dst, *this);
}

// file format is a flat Particle_t array, same as AoS/AoSMapped - reuse their chunked
// reader for cross-format deserialization (e.g. into SoA/SoAVX) instead of duplicating it.
template<>
template<typename Src, std::uint16_t N, typename Fn>
void ParticlesDeserializer<AoSTS, CPU>::readN(std::string const& file_name, Fn fn)
{
    ParticlesDeserializer<AoS, CPU>{}.template readN<Src, N>(file_name, fn);
}


template<>
template<typename Src>
void ParticlesSerializer<AoSPCTS, CPU>::operator()(std::string const& file_name, Src const& src)
{
    TiledParticlesSerializing<AoSPCTS>::write(file_name, src, *this);
}

template<>
template<typename Dst, typename Src>
void ParticlesDeserializer<AoSPCTS, CPU>::operator()(std::string const& file_name, Dst& dst)
{
    TiledParticlesSerializing<AoSPCTS>::template read<Dst, Src>(file_name, dst, *this);
}

template<>
template<typename Src, std::uint16_t N, typename Fn>
void ParticlesDeserializer<AoSPCTS, CPU>::readN(std::string const& file_name, Fn fn)
{
    ParticlesDeserializer<AoS, CPU>{}.template readN<Src, N>(file_name, fn);
}


} // namespace PHARE::core


#endif /* PHARE_CORE_DATA_PARTICLES_SERIALIZING_DETAIL_TS_SERIALIZING */
