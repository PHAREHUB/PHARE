#ifndef MOMENTS_HPP
#define MOMENTS_HPP


#include "core/numerics/interpolator/interpolator.hpp"

#include <ranges>


namespace PHARE
{
namespace core
{
    template<typename Ions>
    void resetMoments(Ions& ions)
    {
        for (auto& pop : ions)
        {
            pop.particleDensity().zero();
            pop.flux().zero();
        }
    }


    // deposits those of `particles` that are in cells owned by `ownership`
    // (a LevelGhostOwnership, see level_ghost_ownership.hpp)
    void depositOwned(auto& interpolate, auto& particles, auto& density, auto& flux,
                      auto const& layout, auto const& ownership, double const coef = 1.)
    {
        auto owned = particles | std::views::filter([&](auto const& particle) {
                         return ownership.owns(particle.iCell);
                     });
        interpolate(owned, density, flux, layout, coef);
    }


    struct DomainDeposit
    {
    };

    template<typename Ownership>
    struct LevelGhostDeposit
    {
        Ownership const& ownership;
    };


    template<typename DepositTag>
    void depositParticles(auto& ions, auto& layout, auto interpolate, DepositTag const& tag)
    {
        for (auto& pop : ions)
        {
            auto& particleDensity = pop.particleDensity();
            auto& flux            = pop.flux();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            {
                auto& partArray = pop.domainParticles();
                interpolate(partArray, particleDensity, flux, layout);
            }
            else // LevelGhostDeposit
                depositOwned(interpolate, pop.levelGhostParticlesOld(), particleDensity, flux,
                             layout, tag.ownership);
        }
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_HPP
