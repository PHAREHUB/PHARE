#ifndef PHARE_ION_UPDATER_H
#define PHARE_ION_UPDATER_H


#include "core/utilities/box/box.h"
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
    using PartIterator      = typename ParticleArray::iterator;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher            = PHARE::core::Pusher<dimension, PartIterator, Electromag, Interpolator,
                                       BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, PartIterator, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
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
    fixMomentGhosts(ions, layout);
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
    auto domainBox = layout.AMRBox();

    auto inDomainBox = [&domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, domainBox);
    };

    auto constexpr partGhostWidth = GridLayout::ghostWidthForParticles();
    auto ghostBox{domainBox};
    ghostBox.grow(partGhostWidth);

    auto ghostSelector = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };


    for (auto& pop : ions)
    {
        ParticleArray tmpDomain;
        ParticleArray tmpPatchGhost;
        ParticleArray tmpLevelGhost;

        // first push all domain particles
        // push them while still inDomainBox
        // accumulate those inDomainBox

        tmpDomain.resize(pop.domainParticles().size());

        auto inRange
            = makeRange(std::begin(pop.domainParticles()), std::end(pop.domainParticles()));
        auto outRange = makeRange(std::begin(tmpDomain), std::end(tmpDomain));

        auto newEnd
            = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_, inDomainBox, layout);

        interpolator_(std::begin(tmpDomain), newEnd, pop.density(), pop.flux(), layout);


        // then push patch and level ghost particles
        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // some will leave the ghost area
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto& inputArray, auto& outputArray) {
            outputArray.resize(inputArray.size());

            auto inRange  = makeRange(std::begin(inputArray), std::end(inputArray));
            auto outRange = makeRange(std::begin(outputArray), std::end(outputArray));

            auto firstGhostOut = pusher_->move(inRange, outRange, em, pop.mass(), interpolator_,
                                               ghostSelector, layout);

            auto endInDomain = std::partition(firstGhostOut, std::end(outputArray), inDomainBox);

            interpolator_(firstGhostOut, endInDomain, pop.density(), pop.flux(), layout);
        };

        // std::cout << "nbr levelghost : " << pop.levelGhostParticles().size() << "\n";
        pushAndAccumulateGhosts(pop.patchGhostParticles(), tmpPatchGhost);
        pushAndAccumulateGhosts(pop.levelGhostParticles(), tmpLevelGhost);
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
    auto constexpr partGhostWidth = GridLayout::ghostWidthForParticles();
    auto domainBox                = layout.AMRBox();
    auto ghostBox{domainBox};
    ghostBox.grow(partGhostWidth);

    auto ghostSelector = [&ghostBox, &domainBox](auto const& part) {
        auto cell = cellAsPoint(part);
        return core::isIn(cell, ghostBox) and !core::isIn(cell, domainBox);
    };

    auto inDomainSelector
        = [&domainBox](auto const& part) { return core::isIn(cellAsPoint(part), domainBox); };


    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();

        auto domainPartRange = makeRange(std::begin(domainParticles), std::end(domainParticles));

        auto firstOutside = pusher_->move(domainPartRange, domainPartRange, em, pop.mass(),
                                          interpolator_, inDomainSelector, layout);

        domainParticles.erase(firstOutside, std::end(domainParticles));


        auto pushAndCopyInDomain = [&](auto& particleArray) {
            auto range = makeRange(std::begin(particleArray), std::end(particleArray));
            auto firstOutGhostBox
                = pusher_->move(range, range, em, pop.mass(), interpolator_, ghostSelector, layout);

            std::copy_if(firstOutGhostBox, std::end(particleArray),
                         std::back_inserter(domainParticles), inDomainSelector);

            particleArray.erase(firstOutGhostBox, std::end(particleArray));
        };


        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());

        interpolator_(std::begin(domainParticles), std::end(domainParticles), pop.density(),
                      pop.flux(), layout);
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_H
