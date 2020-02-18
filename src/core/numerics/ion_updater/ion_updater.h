#ifndef PHARE_ION_UPDATER_H
#define PHARE_ION_UPDATER_H


#include "core/utilities/box/box.h"
#include "core/utilities/particle_selector/particle_selector.h"

#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/boundary_condition/boundary_condition.h"
#include "core/numerics/moments/moments.h"

#include "core/data/ions/ions.h"

#include "initializer/data_provider.h"



#include <memory>

// TODO alpha coef for interpolating new and old levelGhost should be given somehow...


namespace PHARE::core
{
enum class UpdaterMode { moments_only = 1, particles_and_moments = 2 };

template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField          = typename Ions::vecfield_type;
    using ParticleArray     = typename Ions::particle_array_type;
    using ParticleSelector  = typename PHARE::core::ParticleSelector<Box>;
    using PartIterator      = typename ParticleArray::iterator;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher            = PHARE::core::Pusher<dimension, PartIterator, Electromag, Interpolator,
                                       ParticleSelector, BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, PartIterator, Electromag, Interpolator,
                                                 ParticleSelector, BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict& dict)
        : pusher_{makePusher(dict["name"].template to<std::string>())}
    {
    }

    void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                           UpdaterMode = UpdaterMode::particles_and_moments);


    void updateIons(Ions& ions, GridLayout const& layout);


private:
    void updateMomentsOnly_(Ions& ions, Electromag const& em, GridLayout const& layout);

    void updateAll_(Ions& ions, Electromag const& em, GridLayout const& layout);
};




template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 GridLayout const& layout,
                                                                 double dt, UpdaterMode mode)
{
    resetMoments(ions);
    pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

    if (mode == UpdaterMode::moments_only)
    {
        updateMomentsOnly_(ions, em, layout);
    }
    else
    {
        updateAll_(ions, em, layout);
    }
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions, GridLayout const& layout)
{
    setNansOnGhosts(ions, layout);
    ions.computeDensity();
    ions.computeBulkVelocity();
}



template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_(Ions& ions, Electromag const& em,
                                                                  GridLayout const& layout)
{
    auto inDomainSelector = ParticleSelector{layout.AMRBox()};

    for (auto& pop : ions)
    {
        ParticleArray tmpDomain;
        ParticleArray tmpPatchGhost;
        ParticleArray tmpLevelGhost;

        auto pushAndAccumulate = [&](auto& inputArray, auto& outputArray) {
            outputArray.resize(inputArray.size());

            auto inRange  = makeRange(std::begin(inputArray), std::end(inputArray));
            auto outRange = makeRange(std::begin(outputArray), std::end(outputArray));

            auto newEnd = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_,
                                        inDomainSelector, layout);

            interpolator_(std::begin(outputArray), newEnd, pop.density(), pop.flux(), layout);
        };

        pushAndAccumulate(pop.domainParticles(), tmpDomain);
        pushAndAccumulate(pop.patchGhostParticles(), tmpPatchGhost);
        pushAndAccumulate(pop.levelGhostParticles(), tmpLevelGhost);
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_
   evolves moments and particles from time n to n+1
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateAll_(Ions& ions, Electromag const& em,
                                                          GridLayout const& layout)
{
    auto inDomainSelector = ParticleSelector{layout.AMRBox()};

    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeRange(std::begin(domainParticles), std::end(domainParticles));

        auto firstOutside = pusher_->move(domainPartRange, domainPartRange, em, pop.mass(),
                                          interpolator_, inDomainSelector, layout);

        domainParticles.erase(firstOutside, std::end(domainParticles));

        auto pushAndCopyInDomain = [&](auto& particleArray) {
            auto range  = makeRange(std::begin(particleArray), std::end(particleArray));
            auto newEnd = pusher_->move(range, range, em, pop.mass(), interpolator_,
                                        inDomainSelector, layout);
            std::copy(std::begin(particleArray), newEnd, std::back_inserter(domainParticles));
        };


        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());


        interpolator_(std::begin(domainParticles), std::end(domainParticles), pop.density(),
                      pop.flux(), layout);
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_H
