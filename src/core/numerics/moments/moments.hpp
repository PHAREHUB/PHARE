#ifndef MOMENTS_HPP
#define MOMENTS_HPP


#include "core/numerics/interpolator/interpolator.hpp"

#include <ranges>


namespace PHARE
{
namespace core
{
    // zeroes the moments the time loop deposits: population densities and E-node fluxes (the
    // moment-node flux is deposited for outputs only, see HybridLazyMoments)
    template<typename Ions>
    void resetMoments(Ions& ions)
    {
        for (auto& pop : ions)
        {
            pop.particleDensity().zero();
            pop.fluxE().zero();
        }
    }


    // those of `particles` that are in cells owned by `ownership`
    // (a LevelGhostOwnership, see level_ghost_ownership.hpp)
    auto owned(auto& particles, auto const& ownership)
    {
        return particles | std::views::filter([o = &ownership](auto const& particle) {
                   return o->owns(particle.iCell);
               });
    }

    // the time loop's deposit of those of `particles` owned by `ownership`
    void depositOwned(auto& interpolate, auto& particles, auto& density, auto& fluxE,
                      auto const& layout, auto const& ownership, double const coef = 1.)
    {
        interpolate.depositDensityAndFluxE(owned(particles, ownership), density, fluxE, layout,
                                           coef);
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
            auto& fluxE           = pop.fluxE();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
                interpolate.depositDensityAndFluxE(pop.domainParticles(), particleDensity, fluxE,
                                                   layout);
            else // LevelGhostDeposit
                depositOwned(interpolate, pop.levelGhostParticlesOld(), particleDensity, fluxE,
                             layout, tag.ownership);
        }
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_HPP
