#ifndef PHARE_ION_UPDATER_H
#define PHARE_ION_UPDATER_H


#include "core/utilities/box/box.h"
#include "core/utilities/particle_selector/particle_selector.h"

#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/pusher/pusher.h"
#include "core/numerics/pusher/pusher_factory.h"
#include "core/numerics/boundary_condition/boundary_condition.h"

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

    public:
        IonUpdater(PHARE::initializer::PHAREDict& dict)
            : pusher_{makePusher(dict["name"].template to<std::string>())}
        {
        }

        void update(Ions& ions, Electromag const& em, GridLayout const& layout,
                    UpdaterMode = UpdaterMode::particles_and_moments);


    private:
        void updateMomentsOnly_(Ions& ions, Electromag const& em, GridLayout const& layout);

        void updateAll_(Ions& ions, Electromag const& em, GridLayout const& layout);
    };



    template<typename Ions, typename Electromag, typename GridLayout>
    void IonUpdater<Ions, Electromag, GridLayout>::update(Ions& ions, Electromag const& em,
                                                          GridLayout const& layout,
                                                          UpdaterMode mode)
    {
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
    void IonUpdater<Ions, Electromag, GridLayout>::updateMomentsOnly_(Ions& ions,
                                                                      Electromag const& em,
                                                                      GridLayout const& layout)
    {
        //
    }



    template<typename Ions, typename Electromag, typename GridLayout>
    void IonUpdater<Ions, Electromag, GridLayout>::updateAll_(Ions& ions, Electromag const& em,
                                                              GridLayout const& layout)
    {
        for (auto& pop : ions)
        {
            auto& domainParticles = pop.domainParticles();
            auto domainPartRange
                = makeRange(std::begin(domainParticles), std::end(domainParticles));

            // auto firstOutside = pusher_
        }
    }



} // namespace core
} // namespace PHARE


#endif // ION_UPDATER_H
