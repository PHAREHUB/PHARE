#ifndef MOMENTS_HPP
#define MOMENTS_HPP


#include "core/numerics/interpolator/interpolator.hpp"

#include <stdexcept>


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
            pop.chargeDensity().zero();
            pop.flux().zero();
        }
    }


    struct DomainDeposit
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
            auto& particleDensity = pop.particleDensity();
            auto& chargeDensity   = pop.chargeDensity();
            auto& flux            = pop.flux();

            if constexpr (std::is_same_v<DepositTag, DomainDeposit>)
            {
                auto& partArray = pop.domainParticles();
                interpolate(partArray, particleDensity, chargeDensity, flux, layout);
            }
            // else if constexpr (std::is_same_v<DepositTag, LevelGhostDeposit>)
            // {
            //     auto& partArray = pop.levelGhostParticlesOld();
            //     interpolate(partArray, particleDensity, chargeDensity, flux, layout);
            // }
            else
                throw std::runtime_error("unknown deposit tag");
        }
    }

} // namespace core
} // namespace PHARE



#endif // MOMENTS_HPP
