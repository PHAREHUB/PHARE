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

    auto static make_pusher(GridLayout_t const& layout, auto const dt, auto const mass)
    {
        return Pusher{dt, mass, layout};
    }

    class Pusher;

private:
    void updateAndDepositDomain_(auto& ions, auto const& em, auto const& boxing, double const dt);

    void updateAndDepositAll_(auto& ions, auto const& em, auto const& boxing, double const dt);

    Interpolator_t interpolator; // default
};


template<typename Ions>
class IonUpdater1<Ions>::Pusher
{
public:
    void move_particle(auto& particle, auto const& em, auto const& layout,
                       auto& interpolator) const;

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
    template<bool copy_particle>
    void move_domain(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;

    template<bool copy_particle>
    void move_level_ghost(auto& pop, auto const& em, auto const& boxing, auto& interpolator) const;
};


template<typename Ions>
void IonUpdater1<Ions>::updatePopulations(UpdaterMode const mode, auto& ions, auto&&... args)
{
    PHARE_LOG_SCOPE(3, "IonUpdater1::updatePopulations");

    resetMoments(ions);

    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(ions, args...);
    else
        updateAndDepositAll_(ions, args...);
}



/**
 * @brief IonUpdater1<Ions>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
template<typename Ions>
void IonUpdater1<Ions>::updateAndDepositDomain_(auto& ions, auto const& em, auto const& boxing,
                                                double const dt)
{
    bool constexpr copy_particle = true;

    PHARE_LOG_SCOPE(3, "IonUpdater1::updateAndDepositDomain_");

    for (auto& pop : ions)
        Pusher{dt, pop.mass(), boxing.layout}.template //
            move_interpolate_and_sort<copy_particle>(pop, em, boxing, interpolator);
}



template<typename Ions>
void IonUpdater1<Ions>::updateAndDepositAll_( //
    auto& ions, auto const& em, auto const& boxing, double const dt)
{
    bool constexpr copy_particle = false;

    PHARE_LOG_SCOPE(3, "IonUpdater1::updateAndDepositAll_");

    for (auto& pop : ions)
        Pusher{dt, pop.mass(), boxing.layout}.template //
            move_interpolate_and_sort<copy_particle>(pop, em, boxing, interpolator);
}


template<typename Ions>
void IonUpdater1<Ions>::Pusher::move_particle( //
    auto& particle, auto const& em, auto const& layout, auto& interpolator) const
{
    try
    {
        particle.iCell = boris::advance(particle, halfDtOverDl);
    }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PrePush Particle moved 2 cells with delta/vel: ";
        ss << e.delta << "/" << e.vel << std::endl;
        throw DictionaryException{ss.str()};
    }

    auto const& local_em = interpolator(particle, em, layout);
    try
    {
        boris::accelerate(particle, local_em, dto2m);
        particle.iCell = boris::advance(particle, halfDtOverDl);
    }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PostPush Particle moved 2 cells with delta/vel: ";
        ss << e.delta << "/" << e.vel << std::endl;
        DictionaryException ex{ss.str()};

        auto const& [E, B] = local_em;
        for (std::uint16_t i = 0; i < 3; ++i)
            ex("E_" + std::to_string(i), std::to_string(E[i]));
        for (std::uint16_t i = 0; i < 3; ++i)
            ex("B_" + std::to_string(i), std::to_string(B[i]));
        ex("level", std::to_string(layout.levelNumber()));
        throw ex;
    }
};




template<typename Ions>
template<bool copy_particle>
void IonUpdater1<Ions>::Pusher::move_domain( //
    auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
{
    using Particle_t     = std::decay_t<decltype(pop.domainParticles()[0])>;
    using ParticleLoop_t = std::conditional_t<copy_particle, Particle_t, Particle_t&>;


    auto const& layout = boxing.layout;

    auto& domain      = pop.domainParticles();
    auto& patch_ghost = pop.patchGhostParticles();

    for (std::size_t i = 0; i < domain.size(); ++i) // size might change on iteration!
    {
        ParticleLoop_t particle = domain[i];
        auto const oldCell      = particle.iCell;

        move_particle(particle, em, layout, interpolator);

        if (oldCell == particle.iCell)
        {
            interpolator.particleToMesh( //
                particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

            continue;
        }

        bool const isInDomainBox = boxing.isInDomainBox(particle);
        bool const isInNonLevelGhostBox
            = isInDomainBox || boxing.isInNonLevelGhostBox(particle.iCell);
        bool const should_interpolate = isInNonLevelGhostBox;

        if (!should_interpolate)
        {
            if constexpr (!copy_particle)
            {
                domain.swap_last_reduce_by_one(oldCell, i);
                --i; // redo current index as last is now i
            }
            continue;
        }

        interpolator.particleToMesh( //
            particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

        if constexpr (!copy_particle)
            domain.change_icell(particle, oldCell, i);
    }

    // move to patch_ghost
    if constexpr (!copy_particle)
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

    // PHARE_DEBUG_DO({
    //     for (auto const& p : domain)
    //         if (!isIn(p, boxing.domainBox))
    //             throw std::runtime_error("invalid domain");
    //     for (auto const& p : patch_ghost)
    //         if (isIn(p, boxing.domainBox))
    //             throw std::runtime_error("invalid patch ghost");
    // })
}



template<typename Ions>
template<bool copy_particle>
void IonUpdater1<Ions>::Pusher::move_level_ghost( //
    auto& pop, auto const& em, auto const& boxing, auto& interpolator) const
{
    using Particle_t     = std::decay_t<decltype(pop.domainParticles()[0])>;
    using ParticleLoop_t = std::conditional_t<copy_particle, Particle_t, Particle_t&>;

    auto const& layout = boxing.layout;

    auto& domain      = pop.domainParticles();
    auto& level_ghost = pop.levelGhostParticles();

    for (std::size_t i = 0; i < level_ghost.size(); ++i) // size might change on iteration!
    {
        ParticleLoop_t particle = level_ghost[i];
        auto const oldCell      = particle.iCell;

        move_particle(particle, em, layout, interpolator);

        if (oldCell == particle.iCell)
            continue;

        bool const isInDomainBox      = boxing.isInDomainBox(particle);
        bool const should_interpolate = isInDomainBox;

        if (should_interpolate)
            interpolator.particleToMesh( //
                particle, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);

        if constexpr (!copy_particle)
        {
            if (isInDomainBox)
            {
                domain.push_back(particle);
                level_ghost.swap_last_reduce_by_one(oldCell, i);
                --i; // redo current index as last is now i
                continue;
            }

            bool const isInNonLevelGhostBox
                = isInDomainBox || boxing.isInNonLevelGhostBox(particle.iCell);
            bool const isInLevelGhostBox = !isInNonLevelGhostBox;

            if (isInLevelGhostBox)
                level_ghost.change_icell(particle, oldCell, i);
            else
            {
                level_ghost.swap_last_reduce_by_one(oldCell, i);
                --i; // redo current index as last is now i
            }
        }
    }
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
