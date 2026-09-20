#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP

#include "core/def/phare_config.hpp"

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/interpolator/interpolating.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"

#include "core/numerics/pusher/boris/detail/multi_boris.hpp"
#include "core/data/particles/particle_array_exporter.hpp"

#include "ion_updater_def.hpp" // = UpdaterMode

#include <string>
#include <optional>
#include <unordered_map>


namespace PHARE::core
{

/**
 * @brief IonUpdater moves/deposits a single patch's ions/electromag, given directly.
 */
template<typename ParticleArray_t, typename GridLayout>
class IonUpdater
{
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;
    using ParticleRange                = IndexRange<ParticleArray_t>;
    using Selector_t                   = std::function<ParticleRange(ParticleRange&)>;
    using Interpolator                 = PHARE::core::Interpolator<dimension, interp_order>;
    using BoundaryCondition            = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher                       = PHARE::core::BorisPusher<dimension>;

public:
    using Boxing_t = UpdaterCellMapSelectionBoxing<Selector_t, GridLayout>;

    void updatePopulations(auto& ions, auto const& em, Boxing_t const& boxing, double dt,
                           UpdaterMode mode = UpdaterMode::all);

    void reset() { tmp_particles_ = ParticleArray_t{}; }

private:
    void updateAndDepositDomain_(auto& ions, auto const& em, Boxing_t const& boxing);
    void updateAndDepositAll_(auto& ions, auto const& em, Boxing_t const& boxing);

    Pusher pusher_;
    Interpolator interpolator_;
    std::optional<ParticleArray_t> tmp_particles_;
};


template<typename ParticleArray_t, typename GridLayout>
void IonUpdater<ParticleArray_t, GridLayout>::updatePopulations(auto& ions, auto const& em,
                                                                Boxing_t const& boxing, double dt,
                                                                UpdaterMode mode)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    resetMoments(ions);
    pusher_.setMeshAndTimeStep(boxing.layout.meshSize(), dt);

    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(ions, em, boxing);
    else
        updateAndDepositAll_(ions, em, boxing);
}


/**
 * @brief IonUpdater::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
template<typename ParticleArray_t, typename GridLayout>
void IonUpdater<ParticleArray_t, GridLayout>::updateAndDepositDomain_(auto& ions, auto const& em,
                                                                      Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositDomain_");

    auto const& layout = boxing.layout;

    for (auto& pop : ions)
    {
        tmp_particles_ = pop.domainParticles(); // make local copy
        auto& domain   = *tmp_particles_;

        // first push all domain particles twice
        // accumulate those inNonLevelGhostBox
        auto outRange = makeIndexRange(domain);
        auto allowed = outRange = pusher_.move(outRange, outRange, em, pop.mass(), interpolator_,
                                               layout, boxing.noop, boxing.inNonLevelGhostBox);

        interpolator_(allowed, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy

            auto outRange = makeIndexRange(*tmp_particles_);

            auto enteredInDomain = pusher_.move(outRange, outRange, em, pop.mass(), interpolator_,
                                                layout, boxing.inGhostBox, boxing.inDomainBox);

            interpolator_(enteredInDomain, pop.particleDensity(), pop.chargeDensity(), pop.flux(),
                          layout);
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
 * @brief IonUpdater::updateAndDepositAll_
   evolves moments and particles from time n to n+1
 */
template<typename ParticleArray_t, typename GridLayout>
void IonUpdater<ParticleArray_t, GridLayout>::updateAndDepositAll_(auto& ions, auto const& em,
                                                                   Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(1, "IonUpdater::updateAndDepositAll_");

    auto const& layout = boxing.layout;

    // push domain particles, erase from array those leaving domain
    // push level ghost particles that are in ghost area (==ghost box without domain)
    // copy ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in non level ghost box are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeIndexRange(domainParticles);

        auto inDomain = pusher_.move(domainPartRange, domainPartRange, em, pop.mass(),
                                     interpolator_, layout, boxing.noop, boxing.inDomainBox);

        auto now_ghosts = makeRange(domainParticles, inDomain.iend(), domainParticles.size());
        auto const not_level_ghosts = boxing.inNonLevelGhostBox(now_ghosts);

        // copy out new patch ghosts
        auto& patchGhost = pop.patchGhostParticles();
        patchGhost.reserve(patchGhost.size() + not_level_ghosts.size());
        std::copy(not_level_ghosts.begin(), not_level_ghosts.end(), std::back_inserter(patchGhost));

        PHARE_DEBUG_DO({
            auto const outsideGhostBox = boxing.outsideGhostBox(now_ghosts);
            for (auto const& particle : outsideGhostBox)
            {
                PHARE_LOG_LINE_SS(particle);
                auto const nearbyBox = grow(Box(particle.iCell(), particle.iCell()), 3);
                for (auto const& xyz : em.E)
                    if (auto const overlap = nearbyBox * layout.AMRGhostBoxFor(xyz))
                        for (auto const [bix, lix] : layout.amr_lcl_idx(*overlap))
                        {
                            PHARE_LOG_LINE_SS(xyz.name() << " at:" << bix << ":" << xyz(lix));
                        }
            }
            if (outsideGhostBox.size())
                throw core::DictionaryException{}("ID", "Updater::outsideGhostBox");
        })

        domainParticles.erase(now_ghosts); // drop all ghosts

        if (pop.levelGhostParticles().size())
        {
            auto particleRange = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange
                = pusher_.move(particleRange, particleRange, em, pop.mass(), interpolator_, layout,
                               boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(
                domainParticles, [&](auto const& cell) { return isIn(cell, boxing.domainBox); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        interpolator_( //
            domainParticles, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);
        interpolator_( //
            patchGhost, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);
    }
}


/**
 * @brief ParallelIonUpdater moves/deposits ions across every patch of a level in parallel,
 * given a ModelAccessor rather than direct ions/electromag references.
 */
template<typename ParticleArray_t, typename GridLayout>
class ParallelIonUpdater
{
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;
    using Interpolator_t               = Interpolator<dimension, interp_order, /*atomic=*/false>;
    using Interpolating_t              = Interpolating<dimension, interp_order, /*atomic=*/false>;

public:
    using Boxing_t = UpdaterSelectionBoxing<GridLayout>;

    auto constexpr static use_main_thread = MultiBorisOptions{}.use_main_thread;

    void updatePopulations(auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings,
                           double const& dt, UpdaterMode mode = UpdaterMode::all);

    void reset() {}

private:
    void updateAndDepositDomain_(auto& accessor,
                                 std::unordered_map<std::string, Boxing_t> const& boxings);
    void updateAndDepositAll_(auto& accessor,
                              std::unordered_map<std::string, Boxing_t> const& boxings);

    double dt_ = 0;
};


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updatePopulations(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings, double const& dt,
    UpdaterMode mode)
{
    PHARE_LOG_SCOPE(2, "IonUpdater::updatePopulations");

    for (std::size_t i = 0; i < accessor.size(); ++i)
    {
        auto view      = accessor[i];
        auto [ions, _] = view.args;
        resetMoments(ions);
    }
    dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(accessor, boxings);
    else
        updateAndDepositAll_(accessor, boxings);
}


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updateAndDepositDomain_(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdater::updateAndDepositDomain_");

    if (accessor.size() == 0)
        return;

    using Accessor_t = std::remove_reference_t<decltype(accessor)>;
    MultiBoris<Accessor_t, Interpolator_t>{dt_, accessor}.template move<MultiBorisMode::COPY>(
        boxings);
}


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updateAndDepositAll_(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings)
{
    PHARE_LOG_SCOPE(1, "IonUpdater::updateAndDepositAll_");

    if (accessor.size() == 0)
        return;

    using Accessor_t = std::remove_reference_t<decltype(accessor)>;
    MultiBoris<Accessor_t, Interpolator_t>{dt_, accessor}.move(boxings);

    auto post_move_sync = [&](auto const i) mutable {
        auto view                 = accessor[i];
        auto [ions, _]            = view.args;
        auto const patch_id       = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);

        auto const per_pop = [&](auto& pop) {
            auto& domain = pop.domainParticles();
            delete_particles_not_in(domain, patch_boxings.nonLevelGhostBox);
            move_in_ghost_layer(pop.patchGhostParticles(), domain, patch_boxings.domainBox,
                                patch_boxings.nonLevelGhostBox);
            move_in_domain(domain, pop.levelGhostParticles(), patch_boxings.domainBox);
            delete_particles_not_in(pop.levelGhostParticles(), patch_boxings.ghostBox);
            delete_particles_not_in(domain, patch_boxings.domainBox);
        };

        for (auto& pop : ions)
            per_pop(pop);
    };

    auto deposit = [&](auto const i) mutable {
        auto view            = accessor[i];
        auto [ions, _]       = view.args;
        auto const patch_id  = view.patchID();
        auto const& boxing_i = boxings.at(patch_id);
        Interpolating_t interp;
        for (auto& pop : ions)
        {
            interp.particleToMesh(pop.domainParticles(), boxing_i.layout, pop.particleDensity(),
                                  pop.chargeDensity(), pop.flux());
            interp.particleToMesh(pop.patchGhostParticles(), boxing_i.layout, pop.particleDensity(),
                                  pop.chargeDensity(), pop.flux());
        }
    };

    if constexpr (use_main_thread)
    {
        for (std::size_t i = 0; i < accessor.size(); ++i)
            post_move_sync(i);
        for (std::size_t i = 0; i < accessor.size(); ++i)
            deposit(i);
    }
    else
    {
        auto& tp = ThreadPool::INSTANCE();
        for (std::size_t i = 0; i < accessor.size(); ++i)
            tp.async([&post_move_sync, i] { post_move_sync(i); });
        tp.sync();
        for (std::size_t i = 0; i < accessor.size(); ++i)
            tp.async([&deposit, i] { deposit(i); });
        tp.sync();
    }
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
