#ifndef PHARE_ION_UPDATER_H
#define PHARE_ION_UPDATER_H


#include "core/utilities/box/box.h"
#include "core/utilities/particle_selector/particle_selector.h"

#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/boundary_condition/boundary_condition.h"
#include "core/numerics/moments/moments.h"

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

    template<typename GhostFiller>
    void update(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                GhostFiller&& fillGhosts, UpdaterMode = UpdaterMode::particles_and_moments);


private:
    template<typename GhostFiller>
    void updateMomentsOnly_(Ions& ions, Electromag const& em, GridLayout const& layout,
                            GhostFiller&& fillGhosts);

    template<typename GhostFiller>
    void updateAll_(Ions& ions, Electromag const& em, GridLayout const& layout,
                    GhostFiller&& fillGhosts);


    void setNaNsOnGhosts_(Ions& ions, GridLayout const& layout)
    {
        auto ix0 = layout.physicalStartIndex(QtyCentering::primal, Direction::X);
        auto ix1 = layout.physicalEndIndex(QtyCentering::primal, Direction::X);
        auto ix2 = layout.ghostEndIndex(QtyCentering::primal, Direction::X);


        for (auto& pop : ions)
        {
            for (auto ix = 0u; ix < ix0; ++ix) // leftGhostNodes
            {
                auto& density = pop.density();
                auto& flux    = pop.flux();

                auto& fx = flux.getComponent(Component::X);
                auto& fy = flux.getComponent(Component::Y);
                auto& fz = flux.getComponent(Component::Z);

                density(ix) = NAN;
                fx(ix)      = NAN;
                fy(ix)      = NAN;
                fz(ix)      = NAN;
            }

            for (auto ix = ix1 + 1; ix <= ix2; ++ix)
            {
                auto& density = pop.density();
                auto& flux    = pop.flux();

                auto& fx = flux.getComponent(Component::X);
                auto& fy = flux.getComponent(Component::Y);
                auto& fz = flux.getComponent(Component::Z);

                density(ix) = NAN;
                fx(ix)      = NAN;
                fy(ix)      = NAN;
                fz(ix)      = NAN;
            }
        }
    }
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename GhostFiller>
void IonUpdater<Ions, Electromag, GridLayout>::update(Ions& ions, Electromag const& em,
                                                      GridLayout const& layout, double dt,
                                                      GhostFiller&& fillGhosts, UpdaterMode mode)
{
    resetMoments(ions);
    pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

    if (mode == UpdaterMode::moments_only)
    {
        updateMomentsOnly_(ions, em, layout, std::forward<GhostFiller>(fillGhosts));
    }
    else
    {
        updateAll_(ions, em, layout, std::forward<GhostFiller>(fillGhosts));
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename GhostFiller>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_(Ions& ions, Electromag const& em,
                                                                  GridLayout const& layout,
                                                                  GhostFiller&& fillGhosts)
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

    fillGhosts();
    setNaNsOnGhosts_(ions, layout);
    ions.computeDensity();
    ions.computeBulkVelocity();
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename GhostFiller>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_
   evolves moments and particles from time n to n+1
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateAll_(Ions& ions, Electromag const& em,
                                                          GridLayout const& layout,
                                                          GhostFiller&& fillGhosts)
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

    fillGhosts();
    setNaNsOnGhosts_(ions, layout);
    ions.computeDensity();
    ions.computeBulkVelocity();
}



} // namespace PHARE::core


#endif // ION_UPDATER_H
