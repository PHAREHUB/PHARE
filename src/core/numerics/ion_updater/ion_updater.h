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


namespace PHARE
{
namespace core
{
    enum class UpdaterMode { moments_only = 1, particles_and_moments = 2 };

    template<typename Ions, typename Electromag, typename GridLayout>
    class IonUpdater
    {
    private:
        static constexpr auto dimension    = GridLayout::dimension;
        static constexpr auto interp_order = GridLayout::interp_order;

        using Box               = PHARE::core::Box<int, dimension>;
        using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
        using VecField          = typename Ions::vecfield_type;
        using ParticleArray     = typename Ions::particle_array_type;
        using ParticleSelector  = typename PHARE::core::ParticleSelector<Box>;
        using PartIterator      = typename ParticleArray::iterator;
        using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
        using Pusher = PHARE::core::Pusher<dimension, PartIterator, Electromag, Interpolator,
                                           ParticleSelector, BoundaryCondition, GridLayout>;

        constexpr static auto makePusher
            = PHARE::core::PusherFactory::makePusher<dimension, PartIterator, Electromag,
                                                     Interpolator, ParticleSelector,
                                                     BoundaryCondition, GridLayout>;

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
    };




    template<typename Ions, typename Electromag, typename GridLayout>
    template<typename GhostFiller>
    void IonUpdater<Ions, Electromag, GridLayout>::update(Ions& ions, Electromag const& em,
                                                          GridLayout const& layout, double dt,
                                                          GhostFiller&& fillGhosts,
                                                          UpdaterMode mode)
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
    void IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_(Ions& ions,
                                                                      Electromag const& em,
                                                                      GridLayout const& layout,
                                                                      GhostFiller&& fillGhosts)
    {
        //
    }



    template<typename Ions, typename Electromag, typename GridLayout>
    template<typename GhostFiller>
    void IonUpdater<Ions, Electromag, GridLayout>::updateAll_(Ions& ions, Electromag const& em,
                                                              GridLayout const& layout,
                                                              GhostFiller&& fillGhosts)
    {
        auto inDomainSelector = ParticleSelector{layout.AMRBox()};

        for (auto& pop : ions)
        {
            auto& domainParticles = pop.domainParticles();
            auto domainPartRange
                = makeRange(std::begin(domainParticles), std::end(domainParticles));

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

            fillGhosts();

            interpolator_(std::begin(pop.patchGhostParticles()),
                          std::end(pop.patchGhostParticles()), pop.density(), pop.flux(), layout);


            double alpha = 0.5;
            interpolator_(std::begin(pop.levelGhostParticlesNew()),
                          std::end(pop.levelGhostParticlesNew()), pop.density(), pop.flux(), layout,
                          /*coef = */ alpha);


            interpolator_(std::begin(pop.levelGhostParticlesOld()),
                          std::end(pop.levelGhostParticlesOld()), pop.density(), pop.flux(), layout,
                          /*coef = */ (1. - alpha));
        }
    }



} // namespace core
} // namespace PHARE


#endif // ION_UPDATER_H
