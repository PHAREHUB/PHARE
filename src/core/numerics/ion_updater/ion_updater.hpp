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

// after a REF move: domain particles leaving the owned domain become patch ghosts or are
// deleted, level ghosts entering it become domain particles
template<typename ParticleArray_t>
void post_move_sync(ParticleArray_t& domain, ParticleArray_t& patchGhost,
                    ParticleArray_t& levelGhost, auto const& boxing)
{
    auto const& nonOwnedBoxes = boxing.nonOwnedBoxes;

    delete_particles_not_in(domain, boxing.nonLevelGhostBox);
    move_in_ghost_layer(patchGhost, domain, boxing.domainBox, boxing.nonLevelGhostBox);
    if (!nonOwnedBoxes.empty())
    {
        for (auto const& box : nonOwnedBoxes)
            move_in_domain(patchGhost, domain, box);
        delete_particles_in(levelGhost, nonOwnedBoxes);
    }
    move_in_domain(domain, levelGhost, boxing.domainBox);
    delete_particles_not_in(levelGhost, boxing.ghostBox);
    delete_particles_not_in(domain, boxing.domainBox);
}

// AoSCMTS: same result in one pass per tile. A particle's fate only depends on its cell, and
// domain particles were already moved to the tile owning their cell (on_moved), so only the
// cells of a tile's ghost layer outside the patch domain, and its non-owned domain cells, are
// visited, through the tile's cellmap. Level ghosts only live in patch border tiles
template<typename ParticleArray_t>
    requires(ParticleArray_t::layout_mode == LayoutMode::AoSCMTS)
void post_move_sync(ParticleArray_t& domain, ParticleArray_t& patchGhost,
                    ParticleArray_t& levelGhost, auto const& boxing)
{
    auto const& domainBox     = boxing.domainBox;
    auto const& nonOwnedBoxes = boxing.nonOwnedBoxes;
    auto const ghost_cells    = domain.ghost_cells();

    auto const owned = [&](auto const& cell) {
        return isIn(cell, domainBox) and (nonOwnedBoxes.empty() or !isIn(cell, nonOwnedBoxes));
    };
    auto const touches_non_owned = [&](auto const& box) {
        for (auto const& nob : nonOwnedBoxes)
            if (box * nob)
                return true;
        return false;
    };

    std::vector<std::size_t> erased;

    for (std::size_t ti = 0; ti < domain().size(); ++ti)
    {
        auto& tile = domain()[ti];
        auto& ps   = tile();
        if (ps.size() == 0)
            continue;

        auto const reach = grow(*tile, ghost_cells); // the tile's cellmap box
        if (domainBox * reach == reach and !touches_non_owned(*tile))
            continue; // interior tile, every particle is owned

        erased.clear();
        for (auto const& cell : reach)
        {
            if (owned(cell))
                continue;
            // non-owned domain cells and the patch ghost layer -> patch ghost, else deleted
            bool const to_patch_ghost
                = isIn(cell, domainBox) or isIn(cell, boxing.nonLevelGhostBox);
            for (auto const idx : ps.map(cell))
            {
                if (to_patch_ghost)
                    patchGhost.emplace_back(ps[idx]);
                erased.push_back(idx);
            }
        }
        domain.erase_in_tile(ti, erased);
    }

    for (std::size_t ti = 0; ti < levelGhost().size(); ++ti)
    {
        auto& ps = levelGhost()[ti]();
        erased.clear();
        for (std::size_t idx = 0; idx < ps.size(); ++idx)
        {
            // owned -> domain, still in the ghost layer -> kept, else (non-owned, gone) deleted
            auto const& cell = ps[idx].iCell();
            if (owned(cell))
                domain.emplace_back(ps[idx]);
            else if (!isIn(cell, domainBox) and isIn(cell, boxing.ghostBox))
                continue;
            erased.push_back(idx);
        }
        levelGhost.erase_in_tile(ti, erased);
    }

    domain.on_appended();
    patchGhost.on_appended();
    levelGhost.on_appended();
}

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
        auto allowed  = outRange
            = pusher_.move(outRange, outRange, em, pop.charge(), pop.mass(), interpolator_, layout,
                           boxing.noop, boxing.inNonLevelGhostBox);

        interpolator_(allowed, pop.particleDensity(), pop.flux(), layout);

        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inOwnedDomain

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy

            auto outRange = makeIndexRange(*tmp_particles_);

            auto enteredInDomain
                = pusher_.move(outRange, outRange, em, pop.charge(), pop.mass(), interpolator_,
                               layout, boxing.inGhostBox, boxing.inOwnedDomain);

            interpolator_(enteredInDomain, pop.particleDensity(), pop.flux(), layout);
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

        auto inDomain
            = pusher_.move(domainPartRange, domainPartRange, em, pop.charge(), pop.mass(),
                           interpolator_, layout, boxing.noop, boxing.inOwnedDomain);

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
                = pusher_.move(particleRange, particleRange, em, pop.charge(), pop.mass(),
                               interpolator_, layout, boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(domainParticles,
                                           [&](auto const& cell) { return boxing.isOwned(cell); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        interpolator_( //
            domainParticles, pop.particleDensity(), pop.flux(), layout);
        interpolator_( //
            patchGhost, pop.particleDensity(), pop.flux(), layout);
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

protected:
    // per patch steps, shared with mkn_xyz::ParallelIonUpdater (different dispatch)
    static void reset_moments(auto& accessor);
    static void post_move_sync(auto& accessor, auto const& boxings, std::size_t const i);
    static void deposit(auto& accessor, auto const& boxings, std::size_t const i);

    double dt_ = 0;

private:
    void updateAndDepositDomain_(auto& accessor,
                                 std::unordered_map<std::string, Boxing_t> const& boxings);
    void updateAndDepositAll_(auto& accessor,
                              std::unordered_map<std::string, Boxing_t> const& boxings);

    static auto n_tiles_(auto& view)
    {
        auto& ions = std::get<0>(view.args);
        return ions.size() ? tile_count(ions[0].particleDensity()) : std::size_t{0};
    }

    static void post_move_sync_(auto& view, auto const& boxings);
    static void deposit_tile_(auto& view, std::size_t const t);
};


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::reset_moments(auto& accessor)
{
    if constexpr (use_main_thread)
        for (std::size_t i = 0; i < accessor.size(); ++i)
        {
            auto view      = accessor[i];
            auto [ions, _] = view.args;
            resetMoments(ions);
        }
    else
        ThreadPool::INSTANCE().tile_exec(
            accessor.size(), [&](auto const i) { return accessor[i]; },
            [](auto& view) { return n_tiles_(view); },
            [](auto& view, auto const t) { resetMoments(std::get<0>(view.args), t); });
}


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::post_move_sync_(auto& view,
                                                                      auto const& boxings)
{
    auto& [ions, _]           = view.args;
    auto const& patch_boxings = boxings.at(view.patchID());

    for (auto& pop : ions)
        core::post_move_sync(pop.domainParticles(), pop.patchGhostParticles(),
                             pop.levelGhostParticles(), patch_boxings);
}

template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::post_move_sync(auto& accessor,
                                                                     auto const& boxings,
                                                                     std::size_t const i)
{
    auto view = accessor[i];
    post_move_sync_(view, boxings);
}


// only writes the moments of tile t, domain and patch ghosts of a tile in the same task
template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::deposit_tile_(auto& view,
                                                                    std::size_t const t)
{
    Interpolating_t interp;
    for (auto& pop : std::get<0>(view.args))
    {
        interp.particleToMesh_tile(t, pop.domainParticles(), pop.particleDensity(), pop.flux());
        interp.particleToMesh_tile(t, pop.patchGhostParticles(), pop.particleDensity(),
                                   pop.flux());
    }
}

template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::deposit(auto& accessor, auto const& boxings,
                                                              std::size_t const i)
{
    auto view            = accessor[i];
    auto [ions, _]       = view.args;
    auto const patch_id  = view.patchID();
    auto const& boxing_i = boxings.at(patch_id);
    Interpolating_t interp;
    for (auto& pop : ions)
    {
        interp.particleToMesh(pop.domainParticles(), boxing_i.layout, pop.particleDensity(),
                              pop.flux());
        interp.particleToMesh(pop.patchGhostParticles(), boxing_i.layout, pop.particleDensity(),
                              pop.flux());
    }
}


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updatePopulations(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings, double const& dt,
    UpdaterMode mode)
{
    PHARE_LOG_SCOPE(2, "IonUpdater::updatePopulations");

    reset_moments(accessor);
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
    {
        PHARE_LOG_SCOPE(1, "IonUpdater::updateAndDepositAll_::multi_boris");
        MultiBoris<Accessor_t, Interpolator_t>{dt_, accessor}.move(boxings);
    }

    static constexpr bool per_tile_deposit
        = any_in(ParticleArray_t::layout_mode, LayoutMode::AoSPCTS, LayoutMode::AoSCMTS);

    if constexpr (use_main_thread)
    {
        for (std::size_t i = 0; i < accessor.size(); ++i)
            post_move_sync(accessor, boxings, i);
        for (std::size_t i = 0; i < accessor.size(); ++i)
            deposit(accessor, boxings, i);
    }
    else
    {
        auto& tp = ThreadPool::INSTANCE();
        {
            PHARE_LOG_SCOPE(1, "IonUpdater::updateAndDepositAll_::post_move_sync");
            tp.patch_exec( // moves particles between tiles
                accessor.size(), [&](auto const i) { return accessor[i]; },
                [&](auto& view) { post_move_sync_(view, boxings); });
        }
        {
            PHARE_LOG_SCOPE(1, "IonUpdater::updateAndDepositAll_::deposit");
            if constexpr (per_tile_deposit)
                tp.tile_exec(
                    accessor.size(), [&](auto const i) { return accessor[i]; },
                    [](auto& view) { return n_tiles_(view); },
                    [](auto& view, auto const t) { deposit_tile_(view, t); });
            else
            {
                for (std::size_t i = 0; i < accessor.size(); ++i)
                    tp.async([&, i] { deposit(accessor, boxings, i); });
                tp.sync();
            }
        }
    }
}


} // namespace PHARE::core


#if PHARE_HAVE_MKN_GPU

namespace PHARE::core::mkn_xyz
{

/**
 * @brief mkn_xyz::ParallelIonUpdater: core::ParallelIonUpdater steps dispatched with mkn.gpu
 * ThreadedStreamLauncher (one host thread + stream per patch), and GPU deposit kernels
 */
template<typename ParticleArray_t, typename GridLayout>
class ParallelIonUpdater : public core::ParallelIonUpdater<ParticleArray_t, GridLayout>
{
    using Super = core::ParallelIonUpdater<ParticleArray_t, GridLayout>;

    static constexpr auto dimension     = GridLayout::dimension;
    static constexpr auto interp_order  = GridLayout::options.interp_order;
    static constexpr bool atomic_interp = ParticleArray_t::alloc_mode == AllocatorMode::GPU_UNIFIED;
    using Interpolator_t                = Interpolator<dimension, interp_order, atomic_interp>;
    using Interpolating_t               = Interpolating<dimension, interp_order, atomic_interp>;

public:
    using Boxing_t = Super::Boxing_t;

    auto constexpr static use_main_thread = Super::use_main_thread;

    void updatePopulations(auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings,
                           double const& dt, UpdaterMode mode = UpdaterMode::all);

private:
    void updateAndDepositDomain_(auto& accessor,
                                 std::unordered_map<std::string, Boxing_t> const& boxings);
    void updateAndDepositAll_(auto& accessor,
                              std::unordered_map<std::string, Boxing_t> const& boxings);
};


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updatePopulations(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings, double const& dt,
    UpdaterMode mode)
{
    PHARE_LOG_SCOPE(2, "mkn_xyz::IonUpdater::updatePopulations");

    Super::reset_moments(accessor);
    this->dt_ = dt;
    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(accessor, boxings);
    else
        updateAndDepositAll_(accessor, boxings);
}


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updateAndDepositDomain_(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings)
{
    PHARE_LOG_SCOPE(1, "mkn_xyz::IonUpdater::updateAndDepositDomain_");

    if (accessor.size() == 0)
        return;

    using Accessor_t = std::remove_reference_t<decltype(accessor)>;
    MultiBoris<Accessor_t, Interpolator_t> in{this->dt_, accessor};
    in.template move<MultiBorisMode::COPY>(boxings);

    if constexpr (not use_main_thread)
    {
        in.streamer.join();
        in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositDomain.txt");
    }
}


template<typename ParticleArray_t, typename GridLayout>
void ParallelIonUpdater<ParticleArray_t, GridLayout>::updateAndDepositAll_(
    auto& accessor, std::unordered_map<std::string, Boxing_t> const& boxings)
{
    PHARE_LOG_SCOPE(1, "mkn_xyz::IonUpdater::updateAndDepositAll_");

    if (accessor.size() == 0)
        return;

    using Accessor_t = std::remove_reference_t<decltype(accessor)>;
    using Tile_vt    = Accessor_t::Model_t::electromag_type::vecfield_type::field_type::value_type;

    MultiBoris<Accessor_t, Interpolator_t> in{this->dt_, accessor};
    in.move(boxings);

    auto post_move_sync
        = [&](auto const i) mutable { Super::post_move_sync(accessor, boxings, i); };

    if constexpr (use_main_thread)
        for (std::size_t i = 0; i < accessor.size(); ++i)
            post_move_sync(i);
    else
        in.streamer.host(std::move(post_move_sync));

    if constexpr (ParticleArray_t::alloc_mode == AllocatorMode::GPU_UNIFIED)
    {
        auto deposit = [&](auto const i) mutable {
            auto view      = accessor[i];
            auto [ions, _] = view.args;
            for (std::size_t j = 0; j < ions.size(); ++j)
            {
                auto& pop = ions[j];

                if constexpr (ParticleArray_t::layout_mode == LayoutMode::AoSPCTS)
                {
                    // block = (tile, cell), threads stride the cell (see on_tile_cells)
                    auto const deposit_tiles = [&](auto& particles) {
                        using Launcher       = gpu::TileCellLauncher<false>;
                        auto const parts     = *particles;
                        auto const max_cells = Launcher::max_cells(parts);
                        if (max_cells == 0)
                            return;

                        std::size_t const threads = kernel::warp_size();
                        auto pdensity             = pop.particleDensity();
                        auto flux                 = *pop.flux();
                        Launcher{parts().size(), max_cells, threads}.stream(
                            in.streamer.streams[i], [=] __device__() mutable {
                                Interpolating_t::template on_tile_cells<Tile_vt, Launcher>(
                                    parts, flux, pdensity, threads);
                            });
                    };
                    deposit_tiles(pop.domainParticles());
                    deposit_tiles(pop.patchGhostParticles());
                }
                else
                {
                    std::uint32_t const ds = ions.chargeDensity().max_tile_size();
                    auto const domain      = *pop.domainParticles();
                    using Launcher         = gpu::ChunkLauncher<false>;
                    Launcher launcher{1, 0};
                    launcher.b.x = kernel::warp_size();
                    launcher.g.x = domain().size();
                    launcher.ds  = ds * 4 * 8;
                    assert(launcher.ds < 65000);
                    auto pdensity = pop.particleDensity();
                    auto flux     = *pop.flux();
                    launcher.stream(in.streamer.streams[i], [=] __device__() mutable {
                        Interpolating_t::template on_tiles<Tile_vt>(domain, flux, pdensity);
                    });
                }
            }
        };

        if constexpr (use_main_thread)
            for (std::size_t i = 0; i < accessor.size(); ++i)
                deposit(i);
        else
            in.streamer.host(std::move(deposit));
    }
    else
    {
        auto deposit = [&](auto const i) mutable { Super::deposit(accessor, boxings, i); };

        if constexpr (use_main_thread)
            for (std::size_t i = 0; i < accessor.size(); ++i)
                deposit(i);
        else
            in.streamer.host(std::move(deposit));
    }

    if constexpr (not use_main_thread)
    {
        in.streamer.join();
        in.streamer.dump_times(detail::timings_dir_str + "/updateAndDepositAll.txt");
    }


#if PHARE_DEBUG
    for (std::size_t i = 0; i < accessor.size(); ++i)
    {
        auto view                 = accessor[i];
        auto [ions, _]            = view.args;
        auto const patch_id       = view.patchID();
        auto const& patch_boxings = boxings.at(patch_id);

        for (std::size_t j = 0; j < ions.size(); ++j)
        {
            auto& pop = ions[j];
            particle_array_domain_is_valid(pop.domainParticles(), patch_boxings.domainBox);
            particle_array_ghost_is_valid(pop.patchGhostParticles(), patch_boxings.domainBox,
                                          patch_boxings.ghostBox);
            particle_array_ghost_is_valid(pop.levelGhostParticles(), patch_boxings.domainBox,
                                          patch_boxings.ghostBox);
        }
    }
#endif /* PHARE_DEBUG */
}


} // namespace PHARE::core::mkn_xyz

#endif // PHARE_HAVE_MKN_GPU


#endif // ION_UPDATER_HPP
