#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP

#include "core/data/particles/particle_array_def.hpp"
#include "core/data/particles/particle_array_type_options.hpp"

// Mutators
#include "core/data/particles/arrays/particle_array_cm_ts.hpp"

// Impls
#include "core/data/particles/arrays/particle_array_aos.hpp"
#include "core/data/particles/arrays/particle_array_soa.hpp"

#include <cassert>

namespace PHARE::core
{

template<auto opts>
class ParticleArray;

template<auto opts, auto o, auto layout_mode, auto storage_mode>
struct ParticleArrayLayoutResolver;

template<auto opts, auto o = ParticleArrayTypeOptions_t<opts>::FROM(opts)>
class ParticleArrayResolver
{
    using resolver_t = ParticleArrayLayoutResolver<opts, o, opts.layout_mode, opts.storage_mode>;

public:
    struct strings;
    using value_type              = resolver_t::value_type;
    auto static constexpr type_id = strings::type_id;
};

template<auto opts, auto o>
struct ParticleArrayResolver<opts, o>::strings
{
    std::string_view static constexpr alloc_mode   = enum_name(opts.alloc_mode);
    std::string_view static constexpr layout_mode  = enum_name(opts.layout_mode);
    std::string_view static constexpr storage_mode = enum_name(opts.storage_mode);
    std::string_view static constexpr cma          = ",";
    std::string_view static constexpr _dim         = to_string_view_v<std::size_t, opts.dim>;

    auto static constexpr type_id
        = join_string_views_v<_dim, cma, layout_mode, cma, alloc_mode, cma, storage_mode, cma>;
};


template<auto opts, auto o = ParticleArrayTypeOptions_t<opts>::FROM(opts)>
using ResolvedParticleArray_t = ParticleArrayResolver<opts, o>::value_type;


template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoS, StorageMode::VECTOR>
{
    using value_type = AoSParticles<AoSVector, o>;
};
template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoS, StorageMode::ARRAY>
{
    using value_type = AoSParticles<AoSArray, o>;
};
template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoS, StorageMode::SPAN>
{
    using value_type = AoSParticles<AoSSpan, o>;
};

template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoSMapped, StorageMode::VECTOR>
{
    using value_type = AoSMappedParticles<AoSMappedVector, o>;
};

template<auto opts, auto o> // NOT DEFINED CAUSE NOT SENSICAL!
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoSMapped, StorageMode::ARRAY>;

template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoSMapped, StorageMode::SPAN>
{
    using value_type = AoSMappedParticles<AoSMappedSpan, o>;
};

template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoSCMTS, StorageMode::VECTOR>
{
    using Inner      = ParticleArray<opts.with_layout(LayoutMode::AoSMapped)>;
    using value_type = MappedTileSetParticles<MappedTileSetVector<Inner>>;
};

template<auto opts, auto o> // NOT DEFINED CAUSE NOT SENSICAL!
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoSCMTS, StorageMode::ARRAY>;

template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::AoSCMTS, StorageMode::SPAN>
{
    using Inner      = ParticleArray<opts.with_layout(LayoutMode::AoSMapped)>;
    using value_type = MappedTileSetParticles<MappedTileSetSpan<Inner>>;
};

// internal only - see LayoutMode::SoA
template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::SoA, StorageMode::VECTOR>
{
    using value_type = SoAParticles<SoAVector<opts.dim, opts.alloc_mode>>;
};
template<auto opts, auto o> // NOT DEFINED CAUSE NOT SENSICAL!
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::SoA, StorageMode::ARRAY>;
template<auto opts, auto o>
struct ParticleArrayLayoutResolver<opts, o, LayoutMode::SoA, StorageMode::SPAN>
{
    using value_type = SoAParticles<SoASpan<opts.dim, opts.alloc_mode, opts._const_>>;
};

template<typename Particles_t>
void per_particle_iterate(Particles_t& particles, auto const&& fn)
{
    for (auto& p : particles)
        fn(p);
}

template<typename Particles_t>
    requires(Particles_t::layout_mode == LayoutMode::AoS)
void per_particle(Particles_t& particles, auto const&& fn)
{
    per_particle_iterate(particles, std::move(fn));
}

template<typename Particles_t>
    requires(Particles_t::layout_mode == LayoutMode::AoSMapped)
void per_particle(Particles_t& particles, auto const&& fn)
{
    per_particle_iterate(particles, std::move(fn));
}

// one flat (AoSMapped) particle array per tile
template<typename Particles_t>
    requires(Particles_t::layout_mode == LayoutMode::AoSCMTS)
void per_particle(Particles_t& particles, auto const&& fn)
{
    for (auto& tile : particles())
        per_particle_iterate(tile(), std::move(fn));
}

} // namespace PHARE::core

#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_DETAIL_HPP*/
