#ifndef PHARE_CORE_LOL_ION_UPDATER_H
#define PHARE_CORE_LOL_ION_UPDATER_H

#include "core/utilities/box/box.h"
#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/pusher/pusher.h"


#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/boundary_condition/boundary_condition.h"
#include "core/numerics/moments/moments.h"

#include "core/data/ions/ions.h"

#include "initializer/data_provider.h"


#include "core/numerics/ion_updater/ion_updater.h" // for enum class UpdaterMode

#include <memory>


// TODO alpha coef for interpolating new and old levelGhost should be given somehow...


namespace PHARE::core
{
template<typename Ions, typename Electromag, typename GridLayout>
// temporary name
class LOL_IonUpdater
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
    using Pusher = PHARE::core::KirovPusher<dimension, PartIterator, Electromag, Interpolator,
                                            BoundaryCondition, GridLayout>;

    static std::size_t get_threads(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict["pusher"].contains("threads"))
            return dict["pusher"]["threads"].template to<std::size_t>();
        return 1;
    }

public:
    LOL_IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_threads_{get_threads(dict)}
        , pusher_{std::make_unique<Pusher>(pusher_threads_)}
    //, pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
        assert(pusher_threads_ > 0);
    }

    void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                           UpdaterMode = UpdaterMode::particles_and_moments);


    void updateIons(Ions& ions, GridLayout const& layout);


private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, PartIterator, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;


    void updateMomentsOnly_(Ions& ions, Electromag const& em, GridLayout const& layout);

    void updateAll_(Ions& ions, Electromag const& em, GridLayout const& layout);

    std::size_t pusher_threads_ = 1;
    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;
};




template<typename Ions, typename Electromag, typename GridLayout>
void LOL_IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions,
                                                                     Electromag const& em,
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
void LOL_IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions, GridLayout const& layout)
{
    fixMomentGhosts(ions, layout);
    ions.computeDensity();
    ions.computeBulkVelocity();
}



template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief LOL_IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void LOL_IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_(Ions& ions,
                                                                      Electromag const& em,
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
        // first push all domain particles
        // push them while still inDomainBox
        // accumulate those inDomainBox

        {
            auto copy  = pop.domainParticles();
            auto range = makeRange(copy);

            auto newEnd = pusher_->move(range, em, pop.mass(), interpolator_, inDomainBox, layout);

            interpolator_(std::begin(copy), newEnd, pop.density(), pop.flux(), layout);
        }


        // then push patch and level ghost particles
        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // some will leave the ghost area
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto const& ghost_particles) {
            auto copy  = ghost_particles;
            auto range = makeRange(copy);

            auto firstGhostOut
                = pusher_->move(range, em, pop.mass(), interpolator_, ghostSelector, layout);

            auto endInDomain = std::partition(firstGhostOut, std::end(copy), inDomainBox);

            interpolator_(firstGhostOut, endInDomain, pop.density(), pop.flux(), layout);
        };

        pushAndAccumulateGhosts(pop.patchGhostParticles());
        pushAndAccumulateGhosts(pop.levelGhostParticles());
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief LOL_IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_
   evolves moments and particles from time n to n+1
 */
void LOL_IonUpdater<Ions, Electromag, GridLayout>::updateAll_(Ions& ions, Electromag const& em,
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


    // push domain particles, erase from array those leaving domain
    // push patch and level ghost particles that are in ghost area (==ghost box without domain)
    // copy patch and ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in domain are to be interpolated on mesh.


    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();

        auto domainPartRange = makeRange(domainParticles);

        auto firstOutside = pusher_->move(domainPartRange, em, pop.mass(), interpolator_,
                                          inDomainSelector, layout);

        domainParticles.erase(firstOutside, std::end(domainParticles));


        auto pushAndCopyInDomain = [&](auto& particleArray) {
            auto range = makeRange(particleArray);

            auto firstOutGhostBox
                = pusher_->move(range, em, pop.mass(), interpolator_, ghostSelector, layout);

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
