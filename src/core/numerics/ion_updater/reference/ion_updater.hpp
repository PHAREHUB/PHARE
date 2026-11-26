#ifndef PHARE_CORE_NUMERICS_ION_UPDATER_REFERENCE_HPP
#define PHARE_CORE_NUMERICS_ION_UPDATER_REFERENCE_HPP


#include "core/errors.hpp"
#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/ion_updater/ion_updater_def.hpp"

#include "initializer/data_provider.hpp"


#include <array>
#include <string>
#include <sstream>
#include <cstddef>
#include <stdexcept>


// the updater as it was before core::IonUpdater pushed the particle arrays in phases, kept to
// compare against
namespace PHARE::core::reference
{

template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
    using This = IonUpdater<Ions, Electromag, GridLayout>;

public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;

    using Box           = PHARE::core::Box<int, dimension>;
    using Interpolator  = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField      = Ions::vecfield_type;
    using ParticleArray = Ions::particle_array_type;
    using Particle_t    = ParticleArray::Particle_t;
    using PartIterator  = ParticleArray::iterator;
    using ParticleRange = IndexRange<ParticleArray>;

    class Pusher;

private:
    Interpolator interpolator_;

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
    {
        if (dict["pusher"]["name"].template to<std::string>() != "modified_boris")
            throw std::runtime_error("Error : Invalid Pusher name");
    }

    template<typename Boxing_t>
    void updatePopulations(UpdaterMode const mode, Ions& ions, Electromag const& em,
                           Boxing_t const& boxing, double const dt);


    void reset()
    {
        // clear memory
        tmp_particles_ = std::move(ParticleArray{Box{}});
    }


private:
    template<typename Boxing_t>
    void updateAndDepositDomain_(Ions& ions, Electromag const& em, Boxing_t const& boxing,
                                 double const dt);

    template<typename Boxing_t>
    void updateAndDepositAll_(Ions& ions, Electromag const& em, Boxing_t const& boxing,
                              double const dt);


    // dealloced on regridding/load balancing coarsest
    ParticleArray tmp_particles_{Box{}}; //{std::make_unique<ParticleArray>(Box{})};
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

    // This move function should be considered when being used so that all particles are pushed
    // twice - see: https://github.com/PHAREHUB/PHARE/issues/571
    auto move(auto const& rangeIn, auto& rangeOut, auto const& emFields, auto& interpolator,
              auto const& firstSelector, auto const& secondSelector) const;

private:
    void prePushStep_(auto const& rangeIn, auto& rangeOut) const;
    void postPushStep_(auto& range, std::size_t const idx) const;

    GridLayout const layout_;
    double const coef1_;
    std::array<double, dimension> const halfDtOverDl_;
};




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
    {
        updateAndDepositDomain_(ions, em, boxing, dt);
    }
    else
    {
        updateAndDepositAll_(ions, em, boxing, dt);
    }
}

};

/**
 * @brief IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_
   evolves moments from time n to n+1 without updating particles, which stay at time n
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::updateAndDepositDomain_(Ions& ions,
                                                                       Electromag const& em,
                                                                       Boxing_t const& boxing,
                                                                       double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositDomain_");

    auto const& layout = boxing.layout;

    for (auto& pop : ions)
    {
        Pusher const pusher{dt, pop.charge(), pop.mass(), layout};

        auto& domain = (tmp_particles_ = pop.domainParticles()); // make local copy

        // first push all domain particles twice
        // accumulate those inNonLevelGhostBox
        auto outRange = makeIndexRange(domain);
        auto allowed = outRange = pusher.move(outRange, outRange, em, interpolator_, boxing.noop,
                                              boxing.inNonLevelGhostBox);

        interpolator_(allowed, pop.particleDensity(), pop.flux(), layout);


        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inOwnedDomain

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy

            auto outRange = makeIndexRange(tmp_particles_);

            auto enteredInDomain = pusher.move(outRange, outRange, em, interpolator_,
                                               boxing.inGhostBox, boxing.inOwnedDomain);

            interpolator_(enteredInDomain, pop.particleDensity(), pop.flux(), layout);
        };

        // After this function is done domain particles overlaping ghost layers of neighbor patches
        // are sent to these neighbor's patchghost particle array.
        // After being pushed, some patch ghost particles may enter the domain. These need to be
        // copied into the domain array so they are transfered to the neighbor patch
        // ghost array and contribute to moments there too.
        // On the contrary level ghost particles entering the domain here do not need to be copied:
        // this pass only computes moments, they are deposited and the particle arrays are left
        // as they are. Level border nodes also receive the contributions of levelghost old and
        // new particles.

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
                                                                    Boxing_t const& boxing,
                                                                    double const dt)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::updateAndDepositAll_");

    auto const& layout = boxing.layout;

    // push domain particles, erase from array those leaving domain
    // push level ghost particles that are in ghost area (==ghost box without domain)
    // copy ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in non level ghost box are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        Pusher const pusher{dt, pop.charge(), pop.mass(), layout};

        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeIndexRange(domainParticles);

        auto inDomain = pusher.move(domainPartRange, domainPartRange, em, interpolator_,
                                    boxing.noop, boxing.inOwnedDomain);

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
                auto const nearbyBox = grow(Box(particle.iCell, particle.iCell), 3);
                for (auto const& xyz : em.E)
                    if (auto const overlap = nearbyBox * layout.AMRGhostBoxFor(xyz))
                        layout.amr_lcl_idx(*overlap, [&](auto const& bix, auto const& lix) {
                            PHARE_LOG_LINE_SS(xyz.name() << " at:" << bix << ":" << xyz(lix));
                        });
            }
            if (outsideGhostBox.size())
                throw core::DictionaryException{}("ID", "Updater::outsideGhostBox");
        })

        domainParticles.erase(now_ghosts); // drop all ghosts

        if (pop.levelGhostParticles().size())
        {
            auto particleRange     = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange = pusher.move(particleRange, particleRange, em, interpolator_,
                                                 boxing.inGhostBox, boxing.inLevelGhostLayer);

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



template<typename Ions, typename Electromag, typename GridLayout>
auto IonUpdater<Ions, Electromag, GridLayout>::Pusher::move(auto const& rangeIn, auto& rangeOut,
                                                            auto const& emFields,
                                                            auto& interpolator,
                                                            auto const& firstSelector,
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
        auto const& local_em = interpolator(currPart, emFields, layout_);
        boris::accelerate(currPart, local_em, coef1_);

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
            ex("level", std::to_string(layout_.levelNumber()));
            throw ex;
        }
    }

    return secondSelector(rangeOut);
}



template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::prePushStep_(auto const& rangeIn,
                                                                    auto& rangeOut) const
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
        // rangeOut, we need to copy rangeIn weights
        // and velocity. This is done here although
        // not strictly speaking this function's business
        // to take advantage that we're already looping
        // over rangeIn particles.

        outParticles[outIdx].weight = inParticles[inIdx].weight;
        outParticles[outIdx].v      = inParticles[inIdx].v;

        try
        {
            // advance from the input position, the output may hold another one
            auto const newCell
                = boris::advance(inParticles[inIdx], outParticles[outIdx], halfDtOverDl_);
            if (newCell != outParticles[outIdx].iCell)
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

template<typename Ions, typename Electromag, typename GridLayout>
void IonUpdater<Ions, Electromag, GridLayout>::Pusher::postPushStep_(auto& range,
                                                                     std::size_t const idx) const
{
    try
    {
        auto& particles    = range.array();
        auto const newCell = boris::advance(particles[idx], halfDtOverDl_);
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


} // namespace PHARE::core::reference


#endif // PHARE_CORE_NUMERICS_ION_UPDATER_REFERENCE_HPP
