#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TYPE_OPTIONS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TYPE_OPTIONS_HPP

#include "core/data/particles/particle_array_def.hpp"

namespace PHARE::core
{

// VECTOR/SPAN: same shape for every layout - only AoS actually forwards this into its
// resolved type (see ParticleArrayLayoutResolver in particle_array_detail.hpp); other
// layouts just need `o` to exist for ParticleArrayResolver's default template argument.
template<auto opts, auto layout_mode, auto storage_mode>
struct ParticleArrayTypeOptions
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};

// ARRAY: fixed capacity, no allocator
template<auto opts, auto layout_mode>
struct ParticleArrayTypeOptions<opts, layout_mode, StorageMode::ARRAY>
{
    std::size_t dim;
    std::size_t size;

    ParticleArrayTypeOptions static constexpr FROM(auto o, auto&&... args)
    {
        return {o.dim, args...};
    }
};

// per-cell layouts are never ARRAY-backed
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPC, StorageMode::ARRAY>; // nonsense
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPCTS, StorageMode::ARRAY>; // nonsense


template<auto opts>
using ParticleArrayTypeOptions_t
    = ParticleArrayTypeOptions<opts, opts.layout_mode, opts.storage_mode>;


} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TYPE_OPTIONS_HPP*/
