#ifndef PHARE_ION_UPDATER_HPP
#define PHARE_ION_UPDATER_HPP


#include "core/logger.hpp"
#include "core/utilities/box/box.hpp"
#include "core/utilities/range/range.hpp"
#include "core/numerics/pusher/boris.hpp"
#include "core/numerics/pusher/pusher.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/pusher/pusher_factory.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/interpolator/cell_deposit.hpp"
#include "core/numerics/boundary_condition/boundary_condition.hpp"

#include "initializer/data_provider.hpp"


#include <memory>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>


namespace PHARE::core
{
enum class UpdaterMode { domain_only = 1, all = 2 };

template<typename Ions, typename Electromag, typename GridLayout>
class IonUpdater
{
    using This = IonUpdater<Ions, Electromag, GridLayout>;

public:
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;

    using Box               = PHARE::core::Box<int, dimension>;
    using Interpolator      = PHARE::core::Interpolator<dimension, interp_order>;
    using VecField          = Ions::vecfield_type;
    using ParticleArray     = Ions::particle_array_type;
    using Particle_t        = ParticleArray::Particle_t;
    using PartIterator      = ParticleArray::iterator;
    using ParticleRange     = IndexRange<ParticleArray>;
    using BoundaryCondition = PHARE::core::BoundaryCondition<dimension, interp_order>;
    using Pusher = PHARE::core::Pusher<dimension, ParticleRange, Electromag, Interpolator,
                                       BoundaryCondition, GridLayout>;
    using Boris  = PHARE::core::BorisPusher<dimension, ParticleRange, Electromag, Interpolator,
                                            BoundaryCondition, GridLayout>;

private:
    constexpr static auto makePusher
        = PHARE::core::PusherFactory::makePusher<dimension, ParticleRange, Electromag, Interpolator,
                                                 BoundaryCondition, GridLayout>;

    std::unique_ptr<Pusher> pusher_;
    // the domain_only update pushes particle by particle, which needs the concrete pusher
    Boris* boris_ = &boris(*pusher_);
    Interpolator interpolator_;

    static Boris& boris(Pusher& pusher)
    {
        if (auto* boris = dynamic_cast<Boris*>(&pusher))
            return *boris;
        throw std::runtime_error("IonUpdater: the pusher must be a BorisPusher");
    }

public:
    IonUpdater(PHARE::initializer::PHAREDict const& dict)
        : pusher_{makePusher(dict["pusher"]["name"].template to<std::string>())}
        , sortParticles_{cppdict::get_value(dict, "particle_sort", sortByDefault)}
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


    // corrector deposit of the domain particles: cell deposit with the cell kernels on, else by
    // blocks or per particle (hyb::blockedDeposit). A bag layout has holes: it is deposited by
    // its dense runs, which keeps the bags for the next sort (no compaction)
    template<typename Field>
    void depositDomain_(ParticleArray const& particles, Field& density, VecField& flux,
                        GridLayout const& layout)
    {
#if PHARE_HYB_XSIMD
        if (cellKernelsOn_)
        {
            cellDep_.depositArray(particles, interpolator_.depositTargets(density, flux),
                                  Interpolator::amrToLocal(layout));
            return;
        }
#endif
        auto deposit = [&](auto const& range) {
            if constexpr (hyb::blockedDeposit<interp_order>)
            {
                hyb::AoSLanes<dimension, std::decay_t<decltype(range)>> lanes{range, range.size()};
                Interpolator::template depositRangeBlocked<hyb::Backend>(
                    lanes, interpolator_.depositTargets(density, flux), layout);
            }
            else
                interpolator_(range, density, flux, layout);
        };
        if (!particles.bagged())
            return deposit(particles);
        particles.forEachLiveRun([&](std::size_t const first, std::size_t const last) {
            deposit(std::span<Particle_t const>{&particles[first], last - first});
        });
    }

    // domain arrays kept sorted by cell (ParticleArray bag layout) by the corrector (dict
    // "particle_sort"). Sorted arrays have a fixed, cache-friendly particle order whatever their
    // age, and their rows are what the cell kernels need; the sort replaces the cell map
    // bookkeeping of the push. Default: on with the cell kernels (xsimd backend), without them
    // the sort has no kernel to pay for it
    static constexpr bool sortByDefault = PHARE_HYB_XSIMD != 0;
    bool sortParticles_                 = sortByDefault;

    static double domainCells_(GridLayout const& layout)
    {
        auto const box = layout.AMRBox();
        double cells   = 1;
        for (std::size_t d = 0; d < dimension; ++d)
            cells *= static_cast<double>(box.upper[d] - box.lower[d] + 1);
        return cells;
    }

    /** cell kernels for this sorted domain array? Below hyb::cellMinPpc particles per cell on
     * average the blocks of a cell are mostly partial and the kernels walking the rows are
     * faster (-21 % at 8 particles per cell). The choice is kept in the array, made again every
     * hyb::cellChoiceEvery correctors and whenever the array is not sorted (new, restarted or
     * regridded array).
     */
    bool cellKernels_(ParticleArray& domain, GridLayout const& layout, bool const corrector)
    {
        if constexpr (!PHARE_HYB_XSIMD)
            return false;
        auto& choice = domain.cellKernels;
        bool choose  = choice.on < 0 or !domain.csrValid();
        if (corrector)
            choose = (choice.correctors++ % hyb::cellChoiceEvery == 0) or choose;
        if (choose)
        {
            double const ppc = static_cast<double>(domain.liveSize()) / domainCells_(layout);
            choice.on        = ppc >= hyb::cellMinPpc ? 1 : 0;
        }
        return choice.on != 0;
    }

    // kernels of one sorted domain array: the cell kernels on or off for the push (Boris) and
    // the deposits (cellKernelsOn_), back to off on scope exit
    struct KernelScope
    {
        KernelScope(IonUpdater& updater, bool const on)
            : updater_{updater}
        {
            set_(on);
        }
        ~KernelScope() { set_(false); }

    private:
        void set_(bool const on)
        {
            updater_.cellKernelsOn_ = on;
            updater_.boris_->cellKernels(on);
        }
        IonUpdater& updater_;
    };
    bool cellKernelsOn_ = false;

#if PHARE_HYB_XSIMD
    // window scratch of the cell deposit, predictor and corrector
    CellDeposit<dimension, interp_order> cellDep_;
#endif

    // sort scratch, reused across populations and patches
    typename ParticleArray::SortScratch sortScratch_;
    std::vector<std::pair<std::size_t, int>> crossers_;

    template<typename Boxing_t>
    void correctorSorted_(Ions& ions, Electromag const& em, Boxing_t const& boxing);

    template<typename Boxing_t>
    void sortDomain_(ParticleArray& domain, Boxing_t const& boxing, std::size_t oldPrefix);

    // sort key of a domain particle: its owned cell, nbrCells() if not owned (leaver)
    template<typename Boxing_t>
    static auto domainKey_(ParticleArray const& domain, Boxing_t const& boxing)
    {
        return [&domain, &boxing, K = domain.nbrCells()](Particle_t const& p) {
            return boxing.isOwned(p.iCell) ? domain.key(p.iCell) : K;
        };
    }

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
    ions.computeChargeDensity();
    ions.computeBulkVelocity();
    ions.computeChargeFlux();
}

// this is to detach how we partition particles from the updater directly
template<typename IonUpdater_t, typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = IonUpdater_t::Box;
    using Selector_t                     = IonUpdater_t::Pusher::ParticleSelector;

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    // domain cells owned by an overlapping same-level patch, usually empty
    std::vector<Box_t> const nonOwnedBoxes{};
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);

    static bool isOwned(auto const& cell, Box_t const& domainBox,
                        std::vector<Box_t> const& nonOwnedBoxes)
    {
        return isIn(cell, domainBox)
               and (nonOwnedBoxes.empty() or !isIn(Point{cell}, nonOwnedBoxes));
    }

    // a particle belongs to this patch if and only if its cell is in the domain and owned
    bool isOwned(auto const& cell) const { return isOwned(cell, domainBox, nonOwnedBoxes); }

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inOwnedDomain
        = [domainBox = domainBox, nonOwnedBoxes = nonOwnedBoxes](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isOwned(cell, domainBox, nonOwnedBoxes);
              });
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
        // push each domain particle on a local copy, leaving the domain array untouched,
        // and deposit it if it ends in the non level ghost box
        auto const targets    = interpolator_.depositTargets(pop.particleDensity(), pop.flux());
        auto const& keepBoxes = boxing.nonLevelGhostBox;
        double const coef1    = boris_->accelerationCoef(pop.charge(), pop.mass());
        // sorted mode: an array not laid out by a corrector (initial, restarted, regridded,
        // load balanced) gets its rows here, so that the predictor does not depend on where the
        // array comes from (a restarted run is bitwise with the uninterrupted one)
        if (sortParticles_ and !pop.domainParticles().csrValid())
            pop.domainParticles().layoutRows(domainKey_(pop.domainParticles(), boxing),
                                             sortScratch_);
        KernelScope kernels{*this, sortParticles_
                                       and cellKernels_(pop.domainParticles(), layout, false)};
        auto const keep = [&](auto const& cell) { return isIn(Point{cell}, keepBoxes); };
#if PHARE_HYB_XSIMD
        if (cellKernelsOn_)
            boris_->pushEachDepositCell(pop.domainParticles(), em, interpolator_, layout, coef1,
                                        targets, cellDep_, keep);
        else
#endif
            if constexpr (hyb::blockedDeposit<interp_order>)
            boris_->pushEachDeposit(pop.domainParticles(), em, interpolator_, layout, coef1,
                                    targets, keep);
        else
            boris_->pushEach(pop.domainParticles(), em, interpolator_, layout, coef1,
                             [&](Particle_t const& pushed) {
                                 if (isIn(Point{pushed.iCell}, keepBoxes))
                                     interpolator_.depositOn(pushed, targets, layout);
                             });


        // push those in the ghostArea (i.e. stop pushing if they're not out of it)
        // deposit moments on those which leave to go inOwnedDomain

        auto pushAndAccumulateGhosts = [&](auto const& inputArray) {
            tmp_particles_ = inputArray; // work on local copy

            auto outRange = makeIndexRange(tmp_particles_);

            auto enteredInDomain
                = pusher_->move(outRange, outRange, em, pop.charge(), pop.mass(), interpolator_,
                                layout, boxing.inGhostBox, boxing.inOwnedDomain);

            interpolator_(enteredInDomain, pop.particleDensity(), pop.flux(), layout);
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

    if (sortParticles_)
        return correctorSorted_(ions, em, boxing);

    auto const& layout = boxing.layout;

    // push domain particles, erase from array those leaving domain
    // push level ghost particles that are in ghost area (==ghost box without domain)
    // copy ghost particles out of ghost area that are in domain, in particle array
    // finally all particles in non level ghost box are to be interpolated on mesh.
    for (auto& pop : ions)
    {
        auto& domainParticles = pop.domainParticles();
        auto domainPartRange  = makeIndexRange(domainParticles);

        auto inDomain = boris_->moveInPlace(domainPartRange, em, pop.charge(), pop.mass(),
                                            interpolator_, layout, boxing.inOwnedDomain);

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
            auto particleRange = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange
                = pusher_->move(particleRange, particleRange, em, pop.charge(), pop.mass(),
                                interpolator_, layout, boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(domainParticles,
                                           [&](auto const& cell) { return boxing.isOwned(cell); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        depositDomain_(domainParticles, pop.particleDensity(), pop.flux(), layout);
        interpolator_( //
            patchGhost, pop.particleDensity(), pop.flux(), layout);
    }
}


/** bag sort of the domain array by owned cell (key), the leavers to sortScratch_.leavers.
 * Starts from the bag layout valid for [0, oldPrefix) before the push and the crossings the
 * push recorded (no key scan of the sorted part); a full layout without one.
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::sortDomain_(ParticleArray& domain,
                                                           Boxing_t const& boxing,
                                                           std::size_t const oldPrefix)
{
    auto const keyOf = domainKey_(domain, boxing);
    crossers_.clear();
    if (oldPrefix > 0 and domain.bagged())
        boris_->forEachCellChange([&](std::size_t idx, auto const& oldCell, auto const&) {
            if (idx < oldPrefix)
                crossers_.emplace_back(idx, domain.key(oldCell));
        });
    domain.sortBags(keyOf, oldPrefix, crossers_, sortScratch_);
    PHARE_DEBUG_DO({
        if (!domain.checkBags(keyOf))
            throw std::runtime_error("IonUpdater: domain particles not sorted after sortBags");
    })
}


/** corrector with the domain arrays kept sorted:
 * 1. push the domain particles in place (no partition, no cell map: the sort makes the rows)
 * 2. push the level ghosts, append those entering the domain (unsorted appends)
 * 3. bag sort the domain: rows of owned cells, leavers taken out of the array
 * 4. copy the leavers in the non level ghost box to the patch ghosts
 * 5. deposit domain and patch ghosts
 */
template<typename Ions, typename Electromag, typename GridLayout>
template<typename Boxing_t>
void IonUpdater<Ions, Electromag, GridLayout>::correctorSorted_(Ions& ions, Electromag const& em,
                                                                Boxing_t const& boxing)
{
    PHARE_LOG_SCOPE(3, "IonUpdater::correctorSorted_");

    auto const& layout = boxing.layout;

    for (auto& pop : ions)
    {
        auto& domain                = pop.domainParticles();
        KernelScope kernels{*this, cellKernels_(domain, layout, true)};
        std::size_t const oldPrefix = domain.csrPrefix();
        domain.dropMap();

        // the push walks the rows of the layout (Boris pushBlocks_), so the CSR is dropped only
        // after it: the in-place push changes iCell
        auto range = makeIndexRange(domain);
        boris_->moveInPlace(range, em, pop.charge(), pop.mass(), interpolator_, layout,
                            boxing.noop);
        domain.invalidateCSR();

        if (pop.levelGhostParticles().size())
        {
            auto particleRange = makeIndexRange(pop.levelGhostParticles());
            auto inGhostLayerRange
                = pusher_->move(particleRange, particleRange, em, pop.charge(), pop.mass(),
                                interpolator_, layout, boxing.inGhostBox, boxing.inGhostLayer);

            auto& particleArray = particleRange.array();
            particleArray.export_particles(domain,
                                           [&](auto const& cell) { return boxing.isOwned(cell); });

            particleArray.erase(
                makeRange(particleArray, inGhostLayerRange.iend(), particleArray.size()));
        }

        sortDomain_(domain, boxing, oldPrefix);

        // the sort took the leavers out of the array
        auto& patchGhost = pop.patchGhostParticles();
        for (auto const& p : sortScratch_.leavers)
            if (isIn(Point{p.iCell}, boxing.nonLevelGhostBox))
                patchGhost.push_back(p);

        PHARE_DEBUG_DO({
            for (auto const& p : sortScratch_.leavers)
                if (!isIn(p.iCell, boxing.ghostBox))
                    throw core::DictionaryException{}("ID", "Updater::outsideGhostBox");
        })

        depositDomain_(domain, pop.particleDensity(), pop.flux(), layout);
        interpolator_( //
            patchGhost, pop.particleDensity(), pop.flux(), layout);
    }
}



} // namespace PHARE::core


#endif // ION_UPDATER_HPP
