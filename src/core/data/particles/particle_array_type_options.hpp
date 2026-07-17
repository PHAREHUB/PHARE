#ifndef PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TYPE_OPTIONS_HPP
#define PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TYPE_OPTIONS_HPP


#include "core/data/particles/particle_array_def.hpp"



namespace PHARE::core
{


template<auto opts, auto layout_mode, auto storage_mode>
struct ParticleArrayTypeOptions;


template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoS, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoS, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};

template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoS, StorageMode::ARRAY>
{
    std::size_t dim;
    std::size_t size;

    ParticleArrayTypeOptions static constexpr FROM(auto o, auto&&... args)
    {
        return {o.dim, args...};
    }
};


template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSMapped, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSMapped, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSMapped, StorageMode::ARRAY>
{
    std::size_t dim;
    std::size_t size;

    ParticleArrayTypeOptions static constexpr FROM(auto o, auto&&... args)
    {
        return {o.dim, args...};
    }
};


template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPC, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;
    bool thread_safe = false;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPC, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;
    bool thread_safe = false;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPC, StorageMode::ARRAY>;


template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPCTS, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions<opts, LayoutMode::AoSPC, StorageMode::VECTOR> per_cell_opts{
        dim, alloc_mode, true};


    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPCTS, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions<opts, LayoutMode::AoSPC, StorageMode::SPAN> per_cell_opts{
        dim, alloc_mode, true};

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSPCTS, StorageMode::ARRAY>;


template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSTS, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSTS, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSTS, StorageMode::ARRAY>; // nonsense



template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSCMTS, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSCMTS, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::AoSCMTS, StorageMode::ARRAY>
{
    std::size_t dim;
    std::size_t size;

    ParticleArrayTypeOptions static constexpr FROM(auto o, auto&&... args)
    {
        return {o.dim, args...};
    }
};



template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::SoA, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::SoA, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};

template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::SoA, StorageMode::ARRAY>
{
    std::size_t dim;
    std::size_t size;

    ParticleArrayTypeOptions static constexpr FROM(auto o, auto&&... args)
    {
        return {o.dim, args...};
    }
};


template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::SoAVX, StorageMode::VECTOR>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};
template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::SoAVX, StorageMode::SPAN>
{
    std::size_t dim;
    AllocatorMode alloc_mode;

    ParticleArrayTypeOptions static constexpr FROM(auto o) { return {o.dim, o.alloc_mode}; }
};

template<auto opts>
struct ParticleArrayTypeOptions<opts, LayoutMode::SoAVX, StorageMode::ARRAY>
{
    std::size_t dim;
    std::size_t size;

    ParticleArrayTypeOptions static constexpr FROM(auto o, auto&&... args)
    {
        return {o.dim, args...};
    }
};



template<auto opts>
using ParticleArrayTypeOptions_t
    = ParticleArrayTypeOptions<opts, opts.layout_mode, opts.storage_mode>;




} // namespace PHARE::core


#endif /*PHARE_CORE_DATA_PARTICLES_PARTICLE_ARRAY_TYPE_OPTIONS_HPP*/
