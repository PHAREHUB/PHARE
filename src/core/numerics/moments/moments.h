#ifndef MOMENTS_H
#define MOMENTS_H

#include <iterator>

#include "core/numerics/interpolator/interpolator.h"


namespace PHARE
{
namespace core
{
    template<typename Ions, typename GridLayout>
    void
    computeIonMoments(Ions& ions, GridLayout& layout,
                      Interpolator<GridLayout::dimension, GridLayout::interp_order> interpolate)
    {
        for (auto& pop : ions)
        {
            auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
            auto& ghosts                 = pop.patchGhostParticles();
            auto& domain                 = pop.domainParticles();

            auto& density = pop.density();
            auto& flux    = pop.flux();

            // need to reset flux and density before interpolating particles
            flux.zero();
            density.zero();

            interpolate(std::begin(domain), std::end(domain), density, flux, layout);
            interpolate(std::begin(ghosts), std::end(ghosts), density, flux, layout);
            interpolate(std::begin(levelGhostParticlesOld), std::end(levelGhostParticlesOld),
                        density, flux, layout);
        }

        ions.computeDensity();
        ions.computeBulkVelocity();
    }


} // namespace core
} // namespace PHARE



#endif // MOMENTS_H
