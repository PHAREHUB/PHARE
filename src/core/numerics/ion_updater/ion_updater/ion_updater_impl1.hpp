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
    static constexpr auto interp_order = GridLayout_t::options.interp_order;

    using Box                 = PHARE::core::Box<int, dimension>;
    using Interpolator_t      = Interpolator<dimension, interp_order>;
    using BoundaryCondition_t = BoundaryCondition<dimension, interp_order>;

public:
    IonUpdater1(auto&&... /*dict ?*/) {}

    void updatePopulations(UpdaterMode const mode, auto& ions, auto&&... args);

    void reset() {}

    auto static make_pusher(GridLayout_t const& layout, auto const dt, auto const charge,
                            auto const mass)
    {
        return Pusher{dt, charge, mass, layout};
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
    // pushes the population's particles and deposits its moments. With copy_particle the
    // particles are pushed as copies and the arrays are unchanged, else they are pushed in
    // place and moved between the domain, patch ghost and level ghost arrays
    template<bool copy_particle = false>
    void move_and_deposit(auto& pop, auto&&... args) const
    {
        move_domain<copy_particle>(pop, args...);
        move_level_ghost<copy_particle>(pop, args...); // might add to domain
    }

    double const dt;
    double const charge;
    double const mass;
    GridLayout_t const layout;

    double const coef1 = charge * 0.5 * dt / mass;
    std::array<double, dimension> const halfDtOverDl
        = for_N_make_array<dimension>([&](auto i) { return 0.5 * dt / layout.meshSize()[i]; });

private:
    // where a particle may move to: domain particles never leave the ghost box, so whatever
    // throws, every domain particle can be deposited. Level ghosts can go anywhere
    static constexpr auto anywhere = [](auto const& /*cell*/) { return true; };
    static auto inGhostBox(auto const& boxing)
    {
        return [&](auto const& cell) { return isIn(Point{cell}, boxing.ghostBox); };
    }

    // first half push of every particle of the array, in place, the cell map follows
    void pre_push(auto& array, auto const& allowed) const;

    // fields to particles, acceleration and second half push, in place, the cell map follows
    void accel_advance(auto& array, auto const& em, auto& interpolator, auto const& allowed) const;

    template<bool copy_particle>
    void move_domain(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
    {
        if constexpr (copy_particle)
            move_domain_copy(pop, em, boxing, interpolator);
        else
            move_domain_ref(pop, em, boxing, interpolator);
    }

    // pushes copies of the domain particles, by batches, and deposits those still in the
    // domain or patch ghost box, the arrays are unchanged
    void move_domain_copy(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;

    // pushes the domain particles in place, then: owned ones stay, those in a neighbour's ghost
    // box are copied to the patch ghost array, the rest are dropped. Deposits domain and patch
    // ghost particles
    void move_domain_ref(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;

    template<bool copy_particle>
    void move_level_ghost(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
    {
        if constexpr (copy_particle)
            move_level_ghost_copy(pop, em, boxing, interpolator);
        else
            move_level_ghost_ref(pop, em, boxing, interpolator);
    }

    // pushes copies of the level ghosts, by batches, and deposits those entering an owned
    // cell, the arrays are unchanged. Those out of the ghost box after the first half push are
    // skipped
    void move_level_ghost_copy(auto& pop, auto const& em, auto const& boxing,
                               auto& interpolator) const;

    // pushes the level ghosts in place, those out of the ghost box after the first half push are
    // dropped, those entering an owned cell are deposited and moved to the domain array, those
    // leaving the level ghost area are dropped
    void move_level_ghost_ref(auto& pop, auto const& em, auto const& boxing,
                              auto& interpolator) const;

    void advance(auto& particle, auto const& allowed) const;
    void accel_and_advance(auto& particle, auto const& local_em, auto const& allowed) const;

    // a particle is moved only if its new cell is allowed, else this throws first
    static void move_to_(auto& particle, auto const& next, auto const& allowed, char const* step);
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
        Pusher{dt, pop.charge(), pop.mass(), boxing.layout}.template //
            move_and_deposit<true>(pop, em, boxing, interpolator);
}

template<typename Ions>
void IonUpdater1<Ions>::update_ref(auto& ions, auto const& em, auto const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater1::update_ref");

    for (auto& pop : ions)
        Pusher{dt, pop.charge(), pop.mass(), boxing.layout}.template //
            move_and_deposit<false>(pop, em, boxing, interpolator);
}


// ---------------------------------------------------------------------------
// Pusher advance helpers
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_to_(auto& particle, auto const& next, auto const& allowed,
                                         char const* step)
{
    auto const& [delta, cell] = next;
    if (!allowed(cell))
    {
        std::stringstream ss;
        ss << step << " domain particle would leave the ghost box, from cell "
           << Point{particle.iCell} << " to " << Point{cell} << std::endl;
        throw DictionaryException{ss.str()};
    }
    particle.delta = delta;
    particle.iCell = cell;
}

template<typename Ions>
void IonUpdater1<Ions>::Pusher::advance(auto& particle, auto const& allowed) const
{
    try
    {
        move_to_(particle, boris::next_position(particle, halfDtOverDl), allowed, "PrePush");
    }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PrePush Particle moved 2 cells with delta/vel: " << e.delta << "/" << e.vel
           << std::endl;
        throw DictionaryException{ss.str()};
    }
}

template<typename Ions>
void IonUpdater1<Ions>::Pusher::accel_and_advance(auto& particle, auto const& local_em,
                                                  auto const& allowed) const
{
    try
    {
        boris::accelerate(particle, local_em, coef1);
        move_to_(particle, boris::next_position(particle, halfDtOverDl), allowed, "PostPush");
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
void IonUpdater1<Ions>::Pusher::pre_push(auto& array, auto const& allowed) const
{
    for (std::size_t i = 0; i < array.size(); ++i)
    {
        auto& particle     = array[i];
        auto const oldCell = particle.iCell;
        advance(particle, allowed);
        if (particle.iCell != oldCell)
            array.change_icell(particle, oldCell, i);
    }
}

template<typename Ions>
void IonUpdater1<Ions>::Pusher::accel_advance(auto& array, auto const& em, auto& interpolator,
                                              auto const& allowed) const
{
    for (std::size_t i = 0; i < array.size(); ++i)
    {
        auto& particle     = array[i];
        auto const oldCell = particle.iCell;

        auto const& local_em = interpolator(particle, em, layout);
        accel_and_advance(particle, local_em, allowed);
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
            advance(batch[j], inGhostBox(boxing));
        }

        // phase 2: interpolate + accelerate + half advance, compact valid to front
        std::size_t m = 0;
        for (std::size_t j = 0; j < n; ++j)
        {
            auto const& local_em = interpolator(batch[j], em, layout);
            accel_and_advance(batch[j], local_em, inGhostBox(boxing));
            if (boxing.isInNonLevelGhostBox(batch[j]))
                batch[m++] = batch[j];
        }

        // phase 3: deposit valid particles
        for (std::size_t j = 0; j < m; ++j)
            interpolator.particleToMesh(batch[j], pop.particleDensity(), pop.flux(), layout);
    }
}


template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_domain_ref(auto& pop, auto const& em, auto const& boxing,
                                                auto& interpolator) const
{
    auto& domain      = pop.domainParticles();
    auto& patch_ghost = pop.patchGhostParticles();

    pre_push(domain, inGhostBox(boxing));
    accel_advance(domain, em, interpolator, inGhostBox(boxing));

    {
        auto range = makeIndexRange(domain);
        range = domain.partition(range, [&](auto const& cell) { return boxing.isOwned(cell); });

        auto const not_in_domain  = makeRange(domain, range.iend(), domain.size());
        auto const is_patch_ghost = domain.partition(
            not_in_domain, [&](auto const& cell) { return boxing.isInNonLevelGhostBox(cell); });

        patch_ghost.reserve(patch_ghost.size() + is_patch_ghost.size());
        std::copy(is_patch_ghost.begin(), is_patch_ghost.end(), std::back_inserter(patch_ghost));
        domain.erase(not_in_domain);
    }

    for (auto const& p : domain)
        interpolator.particleToMesh(p, pop.particleDensity(), pop.flux(), layout);
    for (auto const& p : patch_ghost)
        interpolator.particleToMesh(p, pop.particleDensity(), pop.flux(), layout);
}


// ---------------------------------------------------------------------------
// Pusher::move_level_ghost
// ---------------------------------------------------------------------------

template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_level_ghost_copy(auto& pop, auto const& em, auto const& boxing,
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

        // phase 1: copy + half advance (no change_icell — working on batch copies), keep those
        // still in the ghost box: fields are only interpolated there
        std::size_t in = 0;
        for (std::size_t j = 0; j < n; ++j)
        {
            batch[in] = level_ghost[base + j];
            advance(batch[in], anywhere);
            if (isIn(Point{batch[in].iCell}, boxing.ghostBox))
                ++in;
        }

        // phase 2: interpolate + accelerate + half advance, compact those entering owned domain
        std::size_t m = 0;
        for (std::size_t j = 0; j < in; ++j)
        {
            auto const& local_em = interpolator(batch[j], em, layout);
            accel_and_advance(batch[j], local_em, anywhere);
            if (boxing.isOwned(batch[j].iCell))
                batch[m++] = batch[j];
        }

        // phase 3: deposit those that entered owned domain
        for (std::size_t j = 0; j < m; ++j)
            interpolator.particleToMesh(batch[j], pop.particleDensity(), pop.flux(), layout);
    }
}


template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_level_ghost_ref(auto& pop, auto const& em, auto const& boxing,
                                                     auto& interpolator) const
{
    auto& domain      = pop.domainParticles();
    auto& level_ghost = pop.levelGhostParticles();

    pre_push(level_ghost, anywhere);

    { // those out of the ghost box are dropped before fields are interpolated on them
        auto range        = makeIndexRange(level_ghost);
        auto const inside = level_ghost.partition(
            range, [&](auto const& cell) { return isIn(Point{cell}, boxing.ghostBox); });
        level_ghost.erase(makeRange(level_ghost, inside.iend(), level_ghost.size()));
    }

    accel_advance(level_ghost, em, interpolator, anywhere);

    {
        auto range = makeIndexRange(level_ghost);

        // partition: stays in level ghost area vs leaves (entered domain or non-level-ghost area,
        // or out of the ghost box)
        auto const stays = level_ghost.partition(range, [&](auto const& cell) {
            return isIn(Point{cell}, boxing.ghostBox) && !core::isIn(cell, boxing.domainBox)
                   && !boxing.isInNonLevelGhostBox(cell);
        });

        auto const leaves = makeRange(level_ghost, stays.iend(), level_ghost.size());

        // only the owner of a cell takes a level ghost entering it
        for (auto const& p : leaves)
            if (boxing.isOwned(p.iCell))
            {
                interpolator.particleToMesh(p, pop.particleDensity(), pop.flux(), layout);
                domain.push_back(p);
            }

        level_ghost.erase(leaves);
    }
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
