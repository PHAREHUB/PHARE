#ifndef PHARE_CORE_NUMERICS_ION_UPDATER0_HPP
#define PHARE_CORE_NUMERICS_ION_UPDATER0_HPP

#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"


#include <array>
#include <cstddef>
#include <iterator>
#include <algorithm>


namespace PHARE::core
{

template<typename Ions>
class IonUpdater0
{
    using GridLayout_t = Ions::gridlayout_type;

public:
    static constexpr auto dimension    = GridLayout_t::dimension;
    static constexpr auto interp_order = GridLayout_t::interp_order;
    using Interpolator_t               = Interpolator<dimension, interp_order>;
    using ParticleArray                = Ions::particle_array_type;

public:
    IonUpdater0(auto&&...) {}

    void updatePopulations( //
        UpdaterMode const, Ions& ions, auto const& em, auto const& boxing, double const dt);

    auto static make_pusher(GridLayout_t const& layout, auto const dt, auto const mass)
    {
        return Pusher{dt, mass, layout};
    }

    void reset()
    {
        tmp_particles_ = std::move(ParticleArray{}); // clear memory
    }

    class Pusher;

private:
    void updateAndDepositDomain_(auto& ions, auto const& em, auto const& boxing, double const dt);


    void updateAndDepositAll_(auto& ions, auto const& em, auto const& boxing, double const dt);

    // dealloced on regridding/load balancing coarsest
    Interpolator_t interpolator;
    ParticleArray tmp_particles_{};
};

template<typename Ions>
class IonUpdater0<Ions>::Pusher
{
public:
    auto move(auto const& rangeIn, auto& rangeOut, auto const& emFields, auto& interpolator,
              auto const& firstSelector, auto const& secondSelector) const;


    double const dt;
    double const mass;
    GridLayout_t const layout;

    double const dto2m = 0.5 * dt / mass;
    std::array<double, dimension> const halfDtOverDl
        = for_N_make_array<dimension>([&](auto i) { return 0.5 * dt / layout.meshSize()[i]; });

private:
    void prePushStep_(auto const& rangeIn, auto& rangeOut) const;
    void postPushStep_(auto& range, std::size_t const idx) const;
};


template<typename Ions>
void IonUpdater0<Ions>::updatePopulations(UpdaterMode mode, Ions& ions, auto const& em,
                                          auto const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater0::updatePopulations");

    resetMoments(ions);

    if (mode == UpdaterMode::domain_only)
        updateAndDepositDomain_(ions, em, boxing, dt);

    else
        updateAndDepositAll_(ions, em, boxing, dt);
}



/**
 * @brief IonUpdater0<Ions>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
template<typename Ions>
void IonUpdater0<Ions>::updateAndDepositDomain_( //
    auto& ions, auto const& em, auto const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater0::updateAndDepositDomain_");

    auto const& layout = boxing.layout;

    for (auto& pop : ions)
    {
        auto const pusher = make_pusher(layout, dt, pop.mass());

        auto& domain = (tmp_particles_ = pop.domainParticles()); // make local copy

        // first push all domain particles twice
        // accumulate those inNonLevelGhostBox
        auto outRange = makeIndexRange(domain);
        auto allowed = outRange = pusher.move( //
            outRange, outRange, em, interpolator, boxing.noop, boxing.inNonLevelGhostBox);

        interpolator(allowed, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);


        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inDomainBox

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy

            auto outRange = makeIndexRange(tmp_particles_);

            auto enteredInDomain = pusher.move( //
                outRange, outRange, em, interpolator, boxing.inGhostBox, boxing.inDomainBox);

            interpolator(enteredInDomain, pop.particleDensity(), pop.chargeDensity(), pop.flux(),
                         layout);
        };

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
 * @brief IonUpdater0<Ions>::updateAndDepositDomain_
   evolves moments and particles from time n to n+1
 */
template<typename Ions>
void IonUpdater0<Ions>::updateAndDepositAll_( //
    auto& ions, auto const& em, auto const& boxing, double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater0::updateAndDepositAll_");

    auto const& layout = boxing.layout;

    // push domain particles, erase from array those leaving domain
    // push level ghost particles that are in ghost area (==ghost box without domain)
    // copy ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in non level ghost box are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto const pusher     = make_pusher(layout, dt, pop.mass());
        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeIndexRange(domainParticles);

        auto inDomain = pusher.move( //
            domainPartRange, domainPartRange, em, interpolator, boxing.noop, boxing.inDomainBox);

        auto now_ghosts = makeRange(domainParticles, inDomain.iend(), domainParticles.size());
        auto const not_level_ghosts = boxing.inNonLevelGhostBox(now_ghosts);

        // copy out new patch ghosts
        auto& patchGhost = pop.patchGhostParticles();
        patchGhost.reserve(patchGhost.size() + not_level_ghosts.size());
        std::copy(not_level_ghosts.begin(), not_level_ghosts.end(), std::back_inserter(patchGhost));

        domainParticles.erase(now_ghosts); // drop all ghosts

        if (pop.levelGhostParticles().size())
        {
            auto particleRange     = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange = pusher.move(particleRange, particleRange, em, interpolator,
                                                 boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(
                domainParticles, [&](auto const& cell) { return isIn(cell, boxing.domainBox); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        interpolator( //
            domainParticles, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);
        interpolator( //
            patchGhost, pop.particleDensity(), pop.chargeDensity(), pop.flux(), layout);
    }
}



template<typename Ions>
auto IonUpdater0<Ions>::Pusher::move(auto const& rangeIn, auto& rangeOut, auto const& emFields,
                                     auto& interpolator, auto const& firstSelector,
                                     auto const& secondSelector) const
{
    PHARE_LOG_SCOPE(3, "Boris::move_no_bc");


    // push the particles of half a step
    // rangeIn : t=n, rangeOut : t=n+1/2
    // Do not partition on this step - this is to keep all domain and ghost
    //   particles consistent. see: https://github.com/PHAREHUB/PHARE/issues/571
    prePushStep_(rangeIn, rangeOut);
    rangeOut = firstSelector(rangeOut);

    for (auto idx = rangeOut.ibegin(); idx < rangeOut.iend(); ++idx)
    {
        auto& currPart = rangeOut.array()[idx];

        //  get electromagnetic fields interpolated on the particles of rangeOut stop at newEnd.
        //  get the particle velocity from t=n to t=n+1
        auto const& local_em = interpolator(currPart, emFields, layout);
        boris::accelerate(currPart, local_em, dto2m);

        // now advance the particles from t=n+1/2 to t=n+1 using v_{n+1} just calculated
        // and get a pointer to the first leaving particle
        try
        {
            postPushStep_(rangeOut, idx);
        }
        catch (DictionaryException const& bex)
        {
            auto ex            = bex;
            auto const& [e, b] = local_em;
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("E_" + std::to_string(i), std::to_string(e[i]));
            for (std::uint16_t i = 0; i < 3; ++i)
                ex("B_" + std::to_string(i), std::to_string(b[i]));
            ex("level", std::to_string(layout.levelNumber()));
            throw ex;
        }
    }

    return secondSelector(rangeOut);
}



template<typename Ions>
void IonUpdater0<Ions>::Pusher::prePushStep_(auto const& rangeIn, auto& rangeOut) const
{
    auto& inParticles  = rangeIn.array();
    auto& outParticles = rangeOut.array();
    for (auto inIdx = rangeIn.ibegin(), outIdx = rangeOut.ibegin(); inIdx < rangeIn.iend();
         ++inIdx, ++outIdx)
    {
        // in the first push, this is the first time
        // we push to rangeOut, which contains crap
        // the push will only touch the particle position
        // but the next step being the acceleration of
        // rangeOut, we need to copy rangeIn weights, charge
        // and velocity. This is done here although
        // not strictly speaking this function's business
        // to take advantage that we're already looping
        // over rangeIn particles.

        outParticles[outIdx].charge = inParticles[inIdx].charge;
        outParticles[outIdx].weight = inParticles[inIdx].weight;
        outParticles[outIdx].v      = inParticles[inIdx].v;

        try
        {
            auto const newCell = boris::advance(outParticles[outIdx], halfDtOverDl);
            if (newCell != inParticles[inIdx].iCell)
                outParticles.change_icell(newCell, outIdx);
        }
        catch (boris::MoveTwoCellException const& e)
        {
            std::stringstream ss;
            ss << "PrePush Particle moved 2 cells with delta/vel: ";
            ss << e.delta << "/" << e.vel << std::endl;
            DictionaryException ex{"cause", ss.str()};
            throw ex;
        }
    }
}

template<typename Ions>
void IonUpdater0<Ions>::Pusher::postPushStep_(auto& range, std::size_t const idx) const
{
    try
    {
        auto& particles    = range.array();
        auto const newCell = boris::advance(particles[idx], halfDtOverDl);
        if (newCell != particles[idx].iCell)
            particles.change_icell(newCell, idx);
    }
    catch (boris::MoveTwoCellException const& e)
    {
        std::stringstream ss;
        ss << "PostPush Particle moved 2 cells with delta/vel: ";
        ss << e.delta << "/" << e.vel << std::endl;
        throw DictionaryException{}("cause", ss.str());
    }
}


} // namespace PHARE::core


#endif // ION_UPDATER_HPP
