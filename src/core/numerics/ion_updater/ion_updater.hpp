#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"

#include "initializer/data_provider.hpp"


#include <array>
#include <string>
#include <sstream>
#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include <type_traits>


namespace PHARE::core
{

// pushes the particle arrays in phases (first half push, then fields to particles, acceleration
// and second half push), then moves the particles between the domain, patch ghost and level ghost
// arrays and deposits them. reference::IonUpdater is the updater it replaces
template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;

    using Box          = PHARE::core::Box<int, dimension>;
    using Interpolator = PHARE::core::Interpolator<dimension, interp_order>;

    class Pusher;

    IonUpdater(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict["pusher"]["name"].template to<std::string>() != "modified_boris")
            throw std::runtime_error("Error : Invalid Pusher name");
    }

    template<typename Boxing_t>
    void updatePopulations(UpdaterMode const mode, Ions& ions, Electromag const& em,
                           Boxing_t const& boxing, double const dt);

    void reset() {}

private:
    template<typename Boxing_t>
    void update_copy_(Ions& ions, Electromag const& em, Boxing_t const& boxing, double const dt);

    template<typename Boxing_t>
    void update_in_place_(Ions& ions, Electromag const& em, Boxing_t const& boxing,
                          double const dt);

    Interpolator interpolator_;
};


template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater<Ions, Electromag, GridLayout>::Pusher
{
public:
    Pusher(double const dt, double const charge, double const mass, GridLayout const& layout)
        : layout_{layout}
        , coef1_{charge * 0.5 * dt / mass}
        , halfDtOverDl_{
              for_N_make_array<dimension>([&](auto i) { return 0.5 * dt / layout.meshSize()[i]; })}
    {
    }

    // pushes the population's particles and deposits its moments. With copy_particle the
    // particles are pushed as copies and the arrays are unchanged, else they are pushed in
    // place and moved between the domain, patch ghost and level ghost arrays
    template<bool copy_particle = false>
    void move_and_deposit(auto& pop, auto&&... args) const
    {
        move_domain<copy_particle>(pop, args...);
        move_level_ghost<copy_particle>(pop, args...); // might add to domain
    }

private:
    // first half push of every particle of the array, in place, the cell map follows
    void pre_push(auto& array) const;

    // fields to particles, acceleration and second half push, in place, the cell map follows
    void accel_advance(auto& array, auto const& em, auto& interpolator) const;

    template<bool copy_particle>
    void move_domain(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
    {
        if constexpr (copy_particle)
            move_domain_copy(pop, em, boxing, interpolator);
        else
            move_domain_in_place(pop, em, boxing, interpolator);
    }

    // pushes copies of the domain particles, by batches, and deposits those still in the
    // domain or patch ghost box, the arrays are unchanged
    void move_domain_copy(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;

    // pushes the domain particles in place, then: owned ones stay, those in a neighbour's ghost
    // box are copied to the patch ghost array, the rest are dropped. Deposits domain and patch
    // ghost particles
    void move_domain_in_place(auto& pop, auto const& em, auto const& boxing,
                              auto& interpolator) const;

    template<bool copy_particle>
    void move_level_ghost(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
    {
        if constexpr (copy_particle)
            move_level_ghost_copy(pop, em, boxing, interpolator);
        else
            move_level_ghost_in_place(pop, em, boxing, interpolator);
    }

    // pushes copies of the level ghosts, by batches, and deposits those entering an owned
    // cell, the arrays are unchanged. Those out of the ghost box after the first half push are
    // skipped
    void move_level_ghost_copy(auto& pop, auto const& em, auto const& boxing,
                               auto& interpolator) const;

    // pushes the level ghosts in place, those out of the ghost box after the first half push are
    // dropped, those entering an owned cell are deposited and moved to the domain array, those
    // leaving the level ghost area are dropped
    void move_level_ghost_in_place(auto& pop, auto const& em, auto const& boxing,
                                   auto& interpolator) const;

    // half push of the particle's position, returns its new cell: iCell is left to the caller
    auto advance(auto& particle) const;
    // acceleration then half push, as advance
    auto accel_and_advance(auto& particle, auto const& local_em) const;

    GridLayout const layout_;
    double const coef1_;
    std::array<double, dimension> const halfDtOverDl_;
};


// ---------------------------------------------------------------------------
// IonUpdater top-level
// ---------------------------------------------------------------------------

template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updatePopulations(UpdaterMode const mode, Ions& ions,
                                                                 Electromag const& em,
                                                                 Boxing_t const& boxing,
                                                                 double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updatePopulations");

    resetMoments(ions);

    if (mode == UpdaterMode::copy)
        update_copy_(ions, em, boxing, dt);
    else
        update_in_place_(ions, em, boxing, dt);
}

template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::update_copy_(Ions& ions, Electromag const& em,
                                                            Boxing_t const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::update_copy_");

    for (auto& pop : ions)
        Pusher{dt, pop.charge(), pop.mass(), boxing.layout}.template //
            move_and_deposit</*copy_particle=*/true>(pop, em, boxing, interpolator_);
}

template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::update_in_place_(Ions& ions, Electromag const& em,
                                                                Boxing_t const& boxing,
                                                                double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::update_in_place_");

    for (auto& pop : ions)
        Pusher{dt, pop.charge(), pop.mass(), boxing.layout}.template //
            move_and_deposit</*copy_particle=*/false>(pop, em, boxing, interpolator_);
}


// ---------------------------------------------------------------------------
// Pusher advance helpers
// ---------------------------------------------------------------------------

template<typename Ions, typename Electromag, typename GridLayout>
auto IonUpdater<Ions, Electromag, GridLayout>::Pusher::advance(auto& particle) const
{
    try
    {
        return boris::advance(particle, halfDtOverDl_);
    }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PrePush Particle moved 2 cells with delta/vel: " << e.delta << "/" << e.vel
           << std::endl;
        throw DictionaryException{ss.str()};
    }
}

template<typename Ions, typename Electromag, typename GridLayout>
auto IonUpdater<Ions, Electromag, GridLayout>::Pusher::accel_and_advance(auto& particle,
                                                                         auto const& local_em) const
{
    try
    {
        boris::accelerate(particle, local_em, coef1_);
        return boris::advance(particle, halfDtOverDl_);
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
        ex("level", std::to_string(layout_.levelNumber()));
        throw ex;
    }
}


// ---------------------------------------------------------------------------
// Pusher shared phase helpers
// ---------------------------------------------------------------------------

template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::pre_push(auto& array) const
{
    for (std::size_t i = 0; i < array.size(); ++i)
    {
        auto& particle = array[i];
        if (auto const newCell = advance(particle); newCell != particle.iCell)
            array.change_icell(newCell, i);
    }
}

template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::accel_advance(auto& array, auto const& em,
                                                                     auto& interpolator) const
{
    for (std::size_t i = 0; i < array.size(); ++i)
    {
        auto& particle       = array[i];
        auto const& local_em = interpolator(particle, em, layout_);
        if (auto const newCell = accel_and_advance(particle, local_em); newCell != particle.iCell)
            array.change_icell(newCell, i);
    }
}


// ---------------------------------------------------------------------------
// Pusher::move_domain
// ---------------------------------------------------------------------------

template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::move_domain_copy(auto& pop, auto const& em,
                                                                        auto const& boxing,
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
            batch[j]       = domain[base + j];
            batch[j].iCell = advance(batch[j]);
        }

        // phase 2: interpolate + accelerate + half advance, compact valid to front
        std::size_t m = 0;
        for (std::size_t j = 0; j < n; ++j)
        {
            auto const& local_em = interpolator(batch[j], em, layout_);
            batch[j].iCell       = accel_and_advance(batch[j], local_em);
            if (boxing.isInNonLevelGhostBox(batch[j].iCell))
                batch[m++] = batch[j];
        }

        // phase 3: deposit valid particles
        for (std::size_t j = 0; j < m; ++j)
            interpolator.particleToMesh(batch[j], pop.particleDensity(), pop.flux(), layout_);
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::move_domain_in_place(
    auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
{
    auto& domain      = pop.domainParticles();
    auto& patch_ghost = pop.patchGhostParticles();

    pre_push(domain);
    accel_advance(domain, em, interpolator);

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
        interpolator.particleToMesh(p, pop.particleDensity(), pop.flux(), layout_);
    for (auto const& p : patch_ghost)
        interpolator.particleToMesh(p, pop.particleDensity(), pop.flux(), layout_);
}


// ---------------------------------------------------------------------------
// Pusher::move_level_ghost
// ---------------------------------------------------------------------------

template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::move_level_ghost_copy(
    auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
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
            batch[in]       = level_ghost[base + j];
            batch[in].iCell = advance(batch[in]);
            if (isIn(Point{batch[in].iCell}, boxing.ghostBox))
                ++in;
        }

        // phase 2: interpolate + accelerate + half advance, compact those entering owned domain
        std::size_t m = 0;
        for (std::size_t j = 0; j < in; ++j)
        {
            auto const& local_em = interpolator(batch[j], em, layout_);
            batch[j].iCell       = accel_and_advance(batch[j], local_em);
            if (boxing.isOwned(batch[j].iCell))
                batch[m++] = batch[j];
        }

        // phase 3: deposit those that entered owned domain
        for (std::size_t j = 0; j < m; ++j)
            interpolator.particleToMesh(batch[j], pop.particleDensity(), pop.flux(), layout_);
    }
}


template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::move_level_ghost_in_place(
    auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
{
    auto& domain      = pop.domainParticles();
    auto& level_ghost = pop.levelGhostParticles();

    // not through the cell map: a half push moves a level ghost up to 2 cells, past the map's box
    // (ghost box + 1), where the map neither sees nor can erase it. Remapped once at the end.
    level_ghost.empty_map();
    auto& particles = level_ghost.vector();

    for (auto& particle : particles)
        particle.iCell = advance(particle);

    // those out of the ghost box are dropped before fields are interpolated on them
    std::erase_if(particles, [&](auto const& p) { return !isIn(Point{p.iCell}, boxing.ghostBox); });

    for (auto& particle : particles)
        particle.iCell = accel_and_advance(particle, interpolator(particle, em, layout_));

    // stays in the level ghost area, or leaves it: entered the domain or the non-level-ghost area,
    // or out of the ghost box
    auto const leaves = std::partition(particles.begin(), particles.end(), [&](auto const& p) {
        return isIn(Point{p.iCell}, boxing.ghostBox) && !core::isIn(p.iCell, boxing.domainBox)
               && !boxing.isInNonLevelGhostBox(p.iCell);
    });

    // only the owner of a cell takes a level ghost entering it
    for (auto it = leaves; it != particles.end(); ++it)
        if (boxing.isOwned(it->iCell))
        {
            interpolator.particleToMesh(*it, pop.particleDensity(), pop.flux(), layout_);
            domain.push_back(*it);
        }

    particles.erase(leaves, particles.end());
    level_ghost.map_particles();
}


} // namespace PHARE::core


#endif // PHARE_ION_UPDATER_HPP
