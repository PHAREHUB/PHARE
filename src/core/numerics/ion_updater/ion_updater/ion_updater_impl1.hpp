#ifndef PHARE_CORE_NUMERICS_ION_UPDATER1_HPP
#define PHARE_CORE_NUMERICS_ION_UPDATER1_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"


namespace PHARE::core
{

template<typename Ions>
class IonUpdater1
{
    using GridLayout_t = Ions::gridlayout_type;

public:
    static constexpr auto dimension    = GridLayout_t::dimension;
    static constexpr auto interp_order = GridLayout_t::interp_order;

    using Box                 = PHARE::core::Box<int, dimension>;
    using Interpolator_t      = Interpolator<dimension, interp_order>;
    using BoundaryCondition_t = BoundaryCondition<dimension, interp_order>;

public:
    IonUpdater1(auto&&... /*dict ?*/) {}

    void updatePopulations(UpdaterMode const mode, auto& ions, auto&&... args);

    void reset() {}

    auto static make_pusher(GridLayout_t const& layout, auto const dt, auto const mass)
    {
        return Pusher{dt, mass, layout};
    }

    class Pusher;

private:
    void update_copy(auto& ions, auto const& em, auto const& boxing, double const dt);
    void update_ref(auto& ions, auto const& em, auto const& boxing, double const dt);

    Interpolator_t interpolator;
};


template<typename Ions>
class IonUpdater1<Ions>::Pusher
{
public:
    template<bool copy_particle = false>
    void move_interpolate_and_sort(auto& pop, auto&&... args) const
    {
        move_domain<copy_particle>(pop, args...);
        move_level_ghost<copy_particle>(pop, args...); // might modify domain
    }

    double const dt;
    double const mass;
    GridLayout_t const layout;

    double const dto2m = 0.5 * dt / mass;
    std::array<double, dimension> const halfDtOverDl
        = for_N_make_array<dimension>([&](auto i) { return 0.5 * dt / layout.meshSize()[i]; });

private:
    // shared phase 1: first half-advance in-place, update cellmap if cell changed
    void pre_push(auto& array) const;

    // shared phase 2: field interpolation + acceleration + second half-advance, update cellmap
    void accel_advance(auto& array, auto const& em, auto& interpolator) const;

    template<bool copy_particle>
    void move_domain(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
    {
        if constexpr (copy_particle)
            move_domain_copy(pop, em, boxing, interpolator);
        else
            move_domain_ref(pop, em, boxing, interpolator);
    }

    // batched 3-phase push for copy mode (no particle array modification)
    void move_domain_copy(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;

    // 2-phase in-place push, then partition domain / patch_ghost, then deposit
    void move_domain_ref(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;

    template<bool copy_particle>
    void move_level_ghost(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
    {
        if constexpr (copy_particle)
            move_level_ghost_copy(pop, em, boxing, interpolator);
        else
            move_level_ghost_ref(pop, em, boxing, interpolator);
    }

    // copy-mode: push a local copy of each level ghost, deposit those entering domain
    void move_level_ghost_copy(auto& pop, auto const& em, auto const& boxing,
                               auto& interpolator) const;

    // 2-phase in-place push, then partition stays / enters-domain / leaves, then deposit
    void move_level_ghost_ref(auto& pop, auto const& em, auto const& boxing,
                              auto& interpolator) const;

    void advance(auto& particle) const;
    void accel_and_advance(auto& particle, auto const& local_em) const;
};


// ---------------------------------------------------------------------------
// IonUpdater1 top-level
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::updatePopulations(UpdaterMode const mode, auto& ions, auto&&... args)
{
    PHARE_LOG_SCOPE(3, "IonUpdater1::updatePopulations");

    resetMoments(ions);

    if (mode == UpdaterMode::copy)
        update_copy(ions, args...);
    else
        update_ref(ions, args...);
}

template<typename Ions>
void IonUpdater1<Ions>::update_copy(auto& ions, auto const& em, auto const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater1::update_copy");

    for (auto& pop : ions)
        Pusher{dt, pop.mass(), boxing.layout}.template //
            move_interpolate_and_sort<true>(pop, em, boxing, interpolator);
}

template<typename Ions>
void IonUpdater1<Ions>::update_ref(auto& ions, auto const& em, auto const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater1::update_ref");

    for (auto& pop : ions)
        Pusher{dt, pop.mass(), boxing.layout}.template //
            move_interpolate_and_sort<false>(pop, em, boxing, interpolator);
}


// ---------------------------------------------------------------------------
// Pusher advance helpers
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::Pusher::advance(auto& particle) const
{
    try { particle.iCell = boris::advance(particle, halfDtOverDl); }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PrePush Particle moved 2 cells with delta/vel: " << e.delta << "/" << e.vel
           << std::endl;
        throw DictionaryException{ss.str()};
    }
}

template<typename Ions>
void IonUpdater1<Ions>::Pusher::accel_and_advance(auto& particle, auto const& local_em) const
{
    try
    {
        boris::accelerate(particle, local_em, dto2m);
        particle.iCell = boris::advance(particle, halfDtOverDl);
    }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PostPush Particle moved 2 cells with delta/vel: " << e.delta << "/" << e.vel
           << std::endl;
        DictionaryException ex{ss.str()};
        auto const& [E, B] = local_em;
        for (std::uint16_t k = 0; k < 3; ++k)
            ex("E_" + std::to_string(k), std::to_string(E[k]));
        for (std::uint16_t k = 0; k < 3; ++k)
            ex("B_" + std::to_string(k), std::to_string(B[k]));
        ex("level", std::to_string(layout.levelNumber()));
        throw ex;
    }
}


// ---------------------------------------------------------------------------
// Pusher shared phase helpers
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::Pusher::pre_push(auto& array) const
{
    for (std::size_t i = 0; i < array.size(); ++i)
    {
        auto& particle     = array[i];
        auto const oldCell = particle.iCell;
        advance(particle);
        if (particle.iCell != oldCell)
            array.change_icell(particle, oldCell, i);
    }
}

template<typename Ions>
void IonUpdater1<Ions>::Pusher::accel_advance(auto& array, auto const& em, auto& interpolator) const
{
    for (std::size_t i = 0; i < array.size(); ++i)
    {
        auto& particle     = array[i];
        auto const oldCell = particle.iCell;

        auto const& local_em = interpolator(particle, em, layout);
        accel_and_advance(particle, local_em);
        if (particle.iCell != oldCell)
            array.change_icell(particle, oldCell, i);
    }
}


// ---------------------------------------------------------------------------
// Pusher::move_domain
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_domain_copy(auto& pop, auto const& em, auto const& boxing,
                                                 auto& interpolator) const
{
    using Particle_t = std::decay_t<decltype(pop.domainParticles()[0])>;

    auto& domain = pop.domainParticles();

    constexpr std::size_t batch_size = 1024;
    std::array<Particle_t, batch_size> batch;
    std::size_t const total = domain.size();

    for (std::size_t base = 0; base < total; base += batch_size)
    {
        std::size_t const n = std::min(batch_size, total - base);

        // phase 1: copy + half advance
        for (std::size_t j = 0; j < n; ++j)
        {
            batch[j] = domain[base + j];
            advance(batch[j]);
        }

        // phase 2: interpolate + accelerate + half advance, compact valid to front
        std::size_t m = 0;
        for (std::size_t j = 0; j < n; ++j)
        {
            auto const& local_em = interpolator(batch[j], em, layout);
            accel_and_advance(batch[j], local_em);
            if (boxing.isInNonLevelGhostBox(batch[j]))
                batch[m++] = batch[j];
        }

        // phase 3: deposit valid particles
        for (std::size_t j = 0; j < m; ++j)
            interpolator.particleToMesh(batch[j], pop.particleDensity(), pop.chargeDensity(),
                                        pop.flux(), layout);
    }
}


template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_domain_ref(auto& pop, auto const& em, auto const& boxing,
                                                auto& interpolator) const
{
    auto& domain      = pop.domainParticles();
    auto& patch_ghost = pop.patchGhostParticles();

    pre_push(domain);
    accel_advance(domain, em, interpolator);

    {
        auto range = makeIndexRange(domain);
        range      = domain.partition(
            range, [&](auto const& cell) { return core::isIn(cell, boxing.domainBox); });

        auto const not_in_domain  = makeRange(domain, range.iend(), domain.size());
        auto const is_patch_ghost = domain.partition(
            not_in_domain, [&](auto const& cell) { return boxing.isInNonLevelGhostBox(cell); });

        patch_ghost.reserve(patch_ghost.size() + is_patch_ghost.size());
        std::copy(is_patch_ghost.begin(), is_patch_ghost.end(), std::back_inserter(patch_ghost));
        domain.erase(not_in_domain);
    }

    for (auto const& p : domain)
        interpolator.particleToMesh(p, pop.particleDensity(), pop.chargeDensity(), pop.flux(),
                                    layout);
    for (auto const& p : patch_ghost)
        interpolator.particleToMesh(p, pop.particleDensity(), pop.chargeDensity(), pop.flux(),
                                    layout);
}


// ---------------------------------------------------------------------------
// Pusher::move_level_ghost
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_level_ghost_copy(auto& pop, auto const& em,
                                                       auto const& boxing,
                                                       auto& interpolator) const
{
    using Particle_t  = std::decay_t<decltype(pop.levelGhostParticles()[0])>;
    auto& level_ghost = pop.levelGhostParticles();

    constexpr std::size_t batch_size = 1024;
    std::array<Particle_t, batch_size> batch;
    std::size_t const total = level_ghost.size();

    for (std::size_t base = 0; base < total; base += batch_size)
    {
        std::size_t const n = std::min(batch_size, total - base);

        // phase 1: copy + half advance (no change_icell — working on batch copies)
        for (std::size_t j = 0; j < n; ++j)
        {
            batch[j] = level_ghost[base + j];
            advance(batch[j]);
        }

        // phase 2: interpolate + accelerate + half advance, compact those entering domain
        std::size_t m = 0;
        for (std::size_t j = 0; j < n; ++j)
        {
            auto const& local_em = interpolator(batch[j], em, layout);
            accel_and_advance(batch[j], local_em);
            if (core::isIn(batch[j].iCell, boxing.domainBox))
                batch[m++] = batch[j];
        }

        // phase 3: deposit those that entered domain
        for (std::size_t j = 0; j < m; ++j)
            interpolator.particleToMesh(batch[j], pop.particleDensity(), pop.chargeDensity(),
                                        pop.flux(), layout);
    }
}


template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_level_ghost_ref(auto& pop, auto const& em, auto const& boxing,
                                                     auto& interpolator) const
{
    auto& domain      = pop.domainParticles();
    auto& level_ghost = pop.levelGhostParticles();

    pre_push(level_ghost);
    accel_advance(level_ghost, em, interpolator);

    {
        auto range = makeIndexRange(level_ghost);

        // partition: stays in level ghost area vs leaves (entered domain or non-level-ghost area)
        auto const stays = level_ghost.partition(range, [&](auto const& cell) {
            return !core::isIn(cell, boxing.domainBox) && !boxing.isInNonLevelGhostBox(cell);
        });

        auto const leaves = makeRange(level_ghost, stays.iend(), level_ghost.size());

        for (auto const& p : leaves)
            if (core::isIn(p.iCell, boxing.domainBox))
            {
                interpolator.particleToMesh(p, pop.particleDensity(), pop.chargeDensity(),
                                            pop.flux(), layout);
                domain.push_back(p);
            }

        level_ghost.erase(leaves);
    }
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
