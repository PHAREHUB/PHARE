#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"

#include "initializer/data_provider.hpp"


#include <memory>


namespace PHARE::core
{
enum class UpdaterMode { domain_only = 1, all = 2 };



template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box           = PHARE::core::Box<int, dimension>;
    using Interpolator  = PHARE::core::Interpolator<dimension, interp_order>;
    using ParticleArray = Ions::particle_array_type;
    using Particle_t    = ParticleArray::Particle_t;
    using PartIterator  = ParticleArray::iterator;

    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher            = BorisPusher<dimension, ParticleArray, Electromag, Interpolator,
                                          BoundaryCondition, GridLayout>;

private:
    Pusher pusher_;
    Interpolator interpolator_;

public:
    IonUpdater() = default;
    IonUpdater(auto const& /*dict*/) {}

    template<typename Boxing_t>
    void updatePopulations(Ions& ions, Electromag const& em, Boxing_t const& boxing, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions);


    void reset() { /* noop */ }


private:
    template<typename Boxing_t>
    void updateCopy(Ions& ions, Electromag const& em, Boxing_t const& boxing);

    template<typename Boxing_t>
    void updateInplace(Ions& ions, Electromag const& em, Boxing_t const& boxing);
};



template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 Boxing_t const& boxing, double dt,
                                                                 UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    resetMoments(ions);
    pusher_.setMeshAndTimeStep(boxing.layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
    {
        updateCopy(ions, em, boxing);
    }
    else
    {
        updateInplace(ions, em, boxing);
    }
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions)
{
    ions.computeChargeDensity();
    ions.computeBulkVelocity();
}

// this is to detach how we partition particles from the updater directly
template<typename IonUpdater_t, typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Particle_t                     = IonUpdater_t::Particle_t;
    using Box_t                          = IonUpdater_t::Box;
    using ParticleArray_t                = IonUpdater_t::ParticleArray;
    using ParticleRange                  = IndexRange<ParticleArray_t>;
    using Selector_t                     = std::function<ParticleRange(ParticleRange&)>;

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    bool isInGhostBox(Particle_t& p) const { return isIn(p, ghostBox); };
    bool isInDomainBox(Particle_t& p) const { return isIn(p, domainBox); };
    bool isInNonLevelGhostBox(auto const& icell) const { return isIn(icell, nonLevelGhostBox); };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inDomainBox = [domainBox = domainBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return core::isIn(cell, domainBox); });
    };

    Selector_t const inGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return isIn(cell, ghostBox); });
    };

    Selector_t const inNonLevelGhostBox
        = [nonLevelGhostBox = nonLevelGhostBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(Point{cell}, nonLevelGhostBox);
              });
          };

    Selector_t const inGhostLayer
        = [ghostBox = ghostBox, domainBox = domainBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(cell, ghostBox) and !isIn(cell, domainBox);
              });
          };

    Selector_t const outsideGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return !isIn(cell, ghostBox); });
    };
};


template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateCopy(Ions& ions, Electromag const& em,
                                                          Boxing_t const& boxing)
{
    bool constexpr copy_particle = true;

    PHARE_LOG_SCOPE(3, "IonUpdater::updateCopy");

    for (auto& pop : ions)
        pusher_.template move<copy_particle>(pop, em, interpolator_, boxing);
}


template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateInplace(Ions& ions, Electromag const& em,
                                                             Boxing_t const& boxing)
{
    bool constexpr copy_particle = false;

    PHARE_LOG_SCOPE(3, "IonUpdater::updateInplace");

    for (auto& pop : ions)
        pusher_.template move<copy_particle>(pop, em, interpolator_, boxing);
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
