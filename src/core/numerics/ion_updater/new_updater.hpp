#ifndef PHARE_NEW_UPDATER_HPP
#define PHARE_NEW_UPDATER_HPP


#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/pusher/fusher.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/data/ions/ions.hpp"

#include "initializer/data_provider.hpp"

#include "core/logger.hpp"

#include "ion_updater.hpp"

#include <cstddef>
#include <memory>


namespace PHARE::core::other
{


template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box           = PHARE::core::Box<int, dimension>;
    using Interpolator  = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField      = typename Ions::vecfield_type;
    using ParticleArray = typename Ions::particle_array_type;
    using Particle_t    = typename ParticleArray::Particle_t;
    using PartIterator  = typename ParticleArray::iterator;

    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::other::Pusher<dimension, ParticleArray, Electromag, Interpolator,
                                              BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher = PHARE::core::other::PusherFactory::makePusher<
        dimension, ParticleArray, Electromag, Interpolator, BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
    }

    void updatePopulations(Ions& ions, Electromag const& em, GridLayout const& layout, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions);


    void reset()
    {
        // clear memory
        tmp_particles_ = std::move(ParticleArray{Box{}});
    }


private:
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, GridLayout const& layout);

    void updateAndDepositAll_(Ions& ions, Electromag const& em, GridLayout const& layout);


    // dealloced on regridding/load balancing coarsest
    ParticleArray tmp_particles_{Box{}}; //{std::make_unique<ParticleArray>(Box{})};
};




template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 GridLayout const& layout,
                                                                 double dt, UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    resetMoments(ions);
    pusher_->setMeshAndTimeStep(layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
    {
        updateAndDepositDomain_(ions, em, layout);
    }
    else
    {
        updateAndDepositAll_(ions, em, layout);
    }
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions)
{
    ions.computeDensity();
    ions.computeBulkVelocity();
}



template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
                                                                       Electromag const& em,
                                                                       GridLayout const& layout)
{
    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();

    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositDomain_");

    auto const domainBox = layout.AMRBox();
    auto const ghostBox  = grow(domainBox, partGhostWidth);

    auto const inDomainBox = [&](auto& particle) { return isIn(particle, domainBox); };
    auto const inGhostBox  = [&](auto& particle) { return isIn(particle, ghostBox); };

    auto const domainChanger
        = [&](auto& particles, auto& particle, auto const& oldCell, auto const idx) {
              if (!inDomainBox(particle))
                  return particles.swap_last_reduce_by_one(oldCell, idx);
              else
                  particles.change_icell(particle, oldCell, idx);
              return false;
          };

    auto const ghostChanger
        = [&](auto& particles, auto& particle, auto const& oldCell, auto const idx) {
              if (!inGhostBox(particle))
                  return particles.swap_last_reduce_by_one(oldCell, idx);
              else
                  particles.change_icell(particle, oldCell, idx);
              return false;
          };

    for (auto& pop : ions)
    {
        auto& domain = pop.domainParticles();

        pusher_->move(domain, em, pop.mass(), interpolator_, layout, domainChanger);
        interpolator_(domain, pop.density(), pop.flux(), layout);

        auto pushAndAccumulateGhosts = [&](auto& particles, bool copyInDomain = false) {
            pusher_->move(particles, em, pop.mass(), interpolator_, layout, ghostChanger, true);

            if (particles.size() == 0)
                return;
            auto enteredInDomain = particles.partition(inDomainBox);
            interpolator_(enteredInDomain, pop.density(), pop.flux(), layout);

            if (copyInDomain)
            {
                domain.reserve(domain.size() + enteredInDomain.size());
                std::copy(enteredInDomain.begin(), enteredInDomain.end(),
                          std::back_inserter(domain));
            }
        };

        pushAndAccumulateGhosts(pop.patchGhostParticles(), true);
        pushAndAccumulateGhosts(tmp_particles_.replace_from(pop.levelGhostParticles()));
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositAll_(Ions& ions,
                                                                    Electromag const& em,
                                                                    GridLayout const& layout)
{
    auto constexpr partGhostWidth = GridLayout::nbrParticleGhosts();

    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositAll_");

    auto const domainBox = layout.AMRBox();
    auto const ghostBox  = grow(domainBox, partGhostWidth);

    auto const inDomainBox = [&](auto& particle) { return isIn(particle, domainBox); };
    auto const inGhostBox  = [&](auto& particle) { return isIn(particle, ghostBox); };

    auto const domainChanger
        = [&](auto& particles, auto& particle, auto const& oldCell, auto const idx) {
              if (!inDomainBox(particle))
                  return particles.swap_last_reduce_by_one(oldCell, idx);
              else
                  particles.change_icell(particle, oldCell, idx);
              return false;
          };

    auto const ghostChanger
        = [&](auto& particles, auto& particle, auto const& oldCell, auto const idx) {
              if (!inGhostBox(particle))
                  return particles.swap_last_reduce_by_one(oldCell, idx);
              else
                  particles.change_icell(particle, oldCell, idx);
              return false;
          };

    auto const inGhostLayer
        = [&](auto& cell) { return isIn(cell, ghostBox) and !isIn(cell, domainBox); };

    for (auto& pop : ions)
    {
        auto& domain = pop.domainParticles();

        pusher_->move(domain, em, pop.mass(), interpolator_, layout, domainChanger);

        auto pushAndCopyInDomain = [&](auto& particles) {
            pusher_->move(particles, em, pop.mass(), interpolator_, layout, ghostChanger, true);

            if (particles.size() == 0)
                return;
            particles.export_particles(domain,
                                       [&](auto const& cell) { return isIn(cell, domainBox); });

            auto inGhostLayerRange = particles.partition(inGhostLayer);
            particles.erase(makeRange(particles, inGhostLayerRange.iend(), particles.size()));
        };

        pushAndCopyInDomain(pop.patchGhostParticles());
        pushAndCopyInDomain(pop.levelGhostParticles());

        interpolator_(makeIndexRange(domain), pop.density(), pop.flux(), layout);
    }
}



} // namespace PHARE::core::other

#endif // PHARE_NEW_UPDATER_HPP
