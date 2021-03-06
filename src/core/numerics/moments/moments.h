#ifndef MOMENTS_H
#define MOMENTS_H

#include <iterator>

#include "core/numerics/interpolator/interpolator.h"


namespace PHARE
{
namespace core
{
    template<typename Ions>
    void resetMoments(Ions& ions)
    {
        for (auto& pop : ions)
        {
            pop.density().zero();
            pop.flux().zero();
        }
    }


    struct DomainDeposit
    {
    };

    struct PatchGhostDeposit
    {
    };
    struct LevelGhostDeposit
    {
    };


    template<typename Ions, typename GridLayout, typename DepositTag>
    void depositParticles(Ions& ions, GridLayout& layout,
                          Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate,
                          DepositTag)
    {
        for (auto& pop : ions)
        {
            auto& density = pop.density();
            auto& flux    = pop.flux();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            {
                auto& partArray = pop.domainParticles();
                interpolate(std::begin(partArray), std::end(partArray), density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, PatchGhostDeposit>)
            {
                auto& partArray = pop.patchGhostParticles();
                interpolate(std::begin(partArray), std::end(partArray), density, flux, layout);
            }
            else if constexpr (std::is_same_v<DepositTag, LevelGhostDeposit>)
            {
                auto& partArray = pop.levelGhostParticlesOld();
                interpolate(std::begin(partArray), std::end(partArray), density, flux, layout);
            }
        }
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_H
