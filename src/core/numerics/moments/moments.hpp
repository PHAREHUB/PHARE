#ifndef MOMENTS_HPP
#define MOMENTS_HPP


#include "core/data/particles/particle_array_def.hpp"
#include "core/data/field/field_tiles.hpp"

#include <type_traits>

#include <stdexcept>


namespace PHARE
{
namespace core
{
    // whole fields, or only tile tile_idx of tiled fields
    template<typename Ions>
    void resetMoments(Ions& ions, auto const&... tile_idx)
    {
        for (auto& pop : ions)
        {
            maybe_tile_at(pop.particleDensity(), tile_idx...).zero();
            maybe_tile_at(pop.flux(), tile_idx...).zero();
        }
    }


    struct DomainDeposit
    {
    };

    struct LevelGhostDeposit
    {
    };


    // whole patch, or only into tile tile_idx of tiled moments
    template<typename DepositTag>
    void depositParticles(auto& ions, auto& layout, auto interpolate, DepositTag,
                          auto const&... tile_idx)
    {
        static_assert(std::is_same_v<DepositTag, DomainDeposit>
                      or std::is_same_v<DepositTag, LevelGhostDeposit>);

        auto constexpr type = std::is_same_v<DepositTag, DomainDeposit> ? ParticleType::Domain
                                                                        : ParticleType::LevelGhost;
        for (auto& pop : ions)
        {
            auto& particleDensity = pop.particleDensity();
            auto& flux            = pop.flux();
            auto& particles       = [&]() -> auto& {
                if constexpr (type == ParticleType::Domain)
                    return pop.domainParticles();
                else
                    return pop.levelGhostParticlesOld();
            }();

            if constexpr (sizeof...(tile_idx) > 0)
                interpolate.template particleToMesh_tile<type>(tile_idx..., particles,
                                                               particleDensity, flux);
            else
                interpolate.template particleToMesh<type>(particles, layout, particleDensity, flux);
        }
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_HPP
