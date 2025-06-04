#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/logger.hpp"
// #include "core/data/ions/ions.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
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
    using This = IonUpdater<Ions, Electromag, GridLayout>;

public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField          = typename Ions::vecfield_type;
    using ParticleArray     = typename Ions::particle_array_type;
    using Particle_t        = typename ParticleArray::Particle_t;
    using PartIterator      = typename ParticleArray::iterator;
    using ParticleRange     = IndexRange<ParticleArray>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::Pusher<dimension, ParticleRange, Electromag, Interpolator,
                                       BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, ParticleRange, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
    {
    }

    template<typename Boxing_t>
    void updatePopulations(Ions& ions, Electromag const& em, Boxing_t const& boxing, double dt,
                           UpdaterMode = UpdaterMode::all);


    void updateIons(Ions& ions);


    void reset()
    {
        // clear memory
        tmp_particles_ = std::move(ParticleArray{Box{}});
    }


private:
    template<typename Boxing_t>
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, Boxing_t const& boxing);

    template<typename Boxing_t>
    void updateAndDepositAll_(Ions& ions, Electromag const& em, Boxing_t const& boxing);


    // dealloced on regridding/load balancing coarsest
    ParticleArray tmp_particles_{Box{}}; //{std::make_unique<ParticleArray>(Box{})};
};




template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(Ions& ions, Electromag const& em,
                                                                 Boxing_t const& boxing, double dt,
                                                                 UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    resetMoments(ions);
    pusher_->setMeshAndTimeStep(boxing.layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
    {
        updateAndDepositDomain_(ions, em, boxing);
    }
    else
    {
        updateAndDepositAll_(ions, em, boxing);
    }
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::updateIons(Ions& ions)
{
    ions.computeDensity();
    ions.computeBulkVelocity();
}

template<typename IonUpdater_t, typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = IonUpdater_t::Box;
    using Selector_t                     = IonUpdater_t::Pusher::ParticleSelector;

    GridLayout_t const layout;
    Box_t const nonLevelGhostBox;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

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
              return particleRange.array().partition(
                  particleRange, [&](auto const& cell) { return isIn(cell, nonLevelGhostBox); });
          };

    Selector_t const inGhostLayer
        = [ghostBox = ghostBox, domainBox = domainBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(cell, ghostBox) and !isIn(cell, domainBox);
              });
          };
};

/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
                                                                       Electromag const& em,
                                                                       Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositDomain_");

    auto const& layout = boxing.layout;

    for (auto& pop : ions)
    {
        auto& domain = (tmp_particles_ = pop.domainParticles()); // make local copy

        // first push all domain particles twice
        // accumulate those inNonLevelGhostBox
        auto outRange = makeIndexRange(domain);
        auto allowed = outRange = pusher_->move(outRange, outRange, em, pop.mass(), interpolator_,
                                                layout, boxing.noop, boxing.inNonLevelGhostBox);

        interpolator_(allowed, pop.density(), pop.flux(), layout);

        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy


            auto outRange = makeIndexRange(tmp_particles_);

            auto enteredInDomain = pusher_->move(outRange, outRange, em, pop.mass(), interpolator_,
                                                 layout, boxing.inGhostBox, boxing.inDomainBox);


            interpolator_(enteredInDomain, pop.density(), pop.flux(), layout);
        };

        // !TODO REVISE!
        // After this function is done domain particles overlaping ghost layers of neighbor patches
        // are sent to these neighbor's patchghost particle array.
        // After being pushed, some patch ghost particles may enter the domain. These need to be
        // copied into the domain array so they are transfered to the neighbor patch
        // ghost array and contribute to moments there too.
        // On the contrary level ghost particles entering the domain here do not need to be copied
        // since they contribute to nodes that are not shared with neighbor patches an since
        // level border nodes will receive contributions from levelghost old and new particles

        if (pop.levelGhostParticles().size())
            pushAndAccumulateGhosts(pop.levelGhostParticles());
    }
}


/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositAll_(Ions& ions,
                                                                    Electromag const& em,
                                                                    Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositAll_");

    auto const& layout = boxing.layout;

    // push domain particles, erase from array those leaving domain
    // push level ghost particles that are in ghost area (==ghost box without domain)
    // copy ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in non level ghost box are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeIndexRange(domainParticles);

        auto inDomain = pusher_->move(domainPartRange, domainPartRange, em, pop.mass(),
                                      interpolator_, layout, boxing.noop, boxing.inDomainBox);

        auto now_ghosts = makeRange(domainParticles, inDomain.iend(), domainParticles.size());
        auto const not_level_ghosts = boxing.inNonLevelGhostBox(now_ghosts);

        // copy out new patch ghosts
        auto& patchGhost = pop.patchGhostParticles();
        patchGhost.reserve(patchGhost.size() + not_level_ghosts.size());
        std::copy(not_level_ghosts.begin(), not_level_ghosts.end(), std::back_inserter(patchGhost));

        domainParticles.erase(now_ghosts); // drop all ghosts

        if (pop.levelGhostParticles().size())
        {
            auto particleRange = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange
                = pusher_->move(particleRange, particleRange, em, pop.mass(), interpolator_, layout,
                                boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(
                domainParticles, [&](auto const& cell) { return isIn(cell, boxing.domainBox); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        interpolator_(makeIndexRange(domainParticles), pop.density(), pop.flux(), layout);
        interpolator_(makeIndexRange(patchGhost), pop.density(), pop.flux(), layout);
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
