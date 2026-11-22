#ifndef PHARE_AMR_PHYSICAL_MODELS_HYBRID_LAZY_MOMENTS_HPP
#define PHARE_AMR_PHYSICAL_MODELS_HYBRID_LAZY_MOMENTS_HPP

#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/numerics/interpolator/interpolator.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/moments/level_ghost_ownership.hpp"

#include "amr/resources_manager/amr_utils.hpp"
#include "mpi/mpi_utils.hpp"
#include "amr/messengers/refiner.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include <SAMRAI/hier/GlobalId.h>
#include <SAMRAI/hier/PatchHierarchy.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <map>
#include <limits>
#include <memory>
#include <ranges>
#include <string>
#include <vector>
#include <stdexcept>

namespace PHARE::amr
{
/**
 * @brief Ion moments the time loop does not need, computed when an output asks for them.
 *
 * Per population: the flux on the moment nodes F_s and the momentum tensor M_s; then the totals:
 * the charge flux on the moment nodes, the bulk velocity and the momentum tensor. F_s and M_s go
 * through the steps the time loop uses for the density and the flux: each patch deposits its
 * domain particles and the level ghost particles it owns (core::LevelGhostOwnership), then the
 * border sums add up the deposits of the patches sharing a node, in place and in the same order on
 * all of them (field_border_sum.hpp). This is valid only when every level is at the same time: the
 * level ghost particles are then levelGhostParticlesOld, which the last fine substep set to the
 * particles at that time. The bulk velocity divides by the mass density the time loop deposited at
 * that time, from the same particles.
 *
 * Results are kept per level and time: a second request at the same time on the same level object
 * does nothing. Filling a level is collective over the MPI ranks (the sum over patch overlaps).
 */
template<typename ResourcesManager, typename Ions, typename GridLayout>
class HybridLazyMoments
{
    static constexpr auto dimension    = GridLayout::dimension;
    static constexpr auto interp_order = GridLayout::options.interp_order;

    using Field = Ions::field_type;

    using FluxSum   = BorderSumRefiner<ResourcesManager, RefinerType::PatchVecFieldBorderSum>;
    using TensorSum = BorderSumRefiner<ResourcesManager, RefinerType::PatchTensorFieldBorderSum>;
    using Ownership = core::LevelGhostOwnership<dimension>;

public:
    struct What
    {
        bool flux           = false;
        bool momentumTensor = false;
    };


    // keyed on the ions' own time stamp: callers may not know the state time (the diagnostics
    // compute before the writer records the dump time)
    void fill(auto& hierarchy, std::shared_ptr<ResourcesManager> const& rmPtr, Ions& ions,
              int const levelNumber, What const what)
    {
        auto& rm         = *rmPtr;
        auto const time  = synchronisedTime_(hierarchy, rm, ions);
        auto const level = hierarchy.getPatchLevel(levelNumber);
        auto& done       = done_[levelNumber];
        // the time is a cache key: same value means same instant, so compare the bits rather than
        // the doubles (CodeQL cpp/equality-on-floats)
        bool const sameLevelAndTime
            = done.level.lock() == level and done.level.lock() != nullptr
              and std::bit_cast<std::uint64_t>(done.time) == std::bit_cast<std::uint64_t>(time);
        if (!sameLevelAndTime)
            done = Done{level, time, What{}};

        What const todo{what.flux and !done.what.flux,
                        what.momentumTensor and !done.what.momentumTensor};
        if (!todo.flux and !todo.momentumTensor)
            return;

        declareSums_(rmPtr, ions);
        auto const& ownership = ownershipOf_(hierarchy, levelNumber);

        core::Interpolator<dimension, interp_order> interpolate;
        core::MomentumTensorInterpolator<dimension, interp_order> interpolateMomentumTensor;

        for (auto& patch : rm.enumerate(*level, ions, density_))
        {
            auto const layout = layoutFromPatch<GridLayout>(*patch);
            auto const* ghostOwnership
                = levelNumber > 0 ? &ownership.at(patch->getGlobalId()) : nullptr;

            for (auto& pop : ions)
            {
                if (todo.flux)
                {
                    // the density goes to a scratch field: the time loop owns the populations'
                    // densities
                    pop.flux().zero();
                    density_.zero();
                    interpolate(pop.domainParticles(), density_, pop.flux(), layout);
                    if (ghostOwnership)
                        interpolate(core::owned(pop.levelGhostParticlesOld(), *ghostOwnership),
                                    density_, pop.flux(), layout);
                }
                if (todo.momentumTensor)
                {
                    pop.momentumTensor().zero();
                    interpolateMomentumTensor(pop.domainParticles(), pop.momentumTensor(), layout,
                                              pop.mass());
                    if (ghostOwnership)
                        interpolateMomentumTensor(
                            core::owned(pop.levelGhostParticlesOld(), *ghostOwnership),
                            pop.momentumTensor(), layout, pop.mass());
                }
            }
        }

        for (std::size_t i = 0; i < ions.size(); ++i)
        {
            if (todo.flux)
                fluxSums_[i].fillTransient(nullptr, level, time);
            if (todo.momentumTensor)
                tensorSums_[i].fillTransient(nullptr, level, time);
        }

        for ([[maybe_unused]] auto& patch : rm.enumerate(*level, ions))
        {
            if (todo.flux)
            {
                ions.computeChargeFluxOnMomentNodes();
                ions.computeBulkVelocity();
            }
            if (todo.momentumTensor)
                ions.computeFullMomentumTensor();
        }

        done.what
            = What{done.what.flux or todo.flux, done.what.momentumTensor or todo.momentumTensor};
    }


private:
    // the time stamp of the ions, which must be the same on every level and every rank.
    // The earliest and latest times are reduced over all ranks before any decision, so every rank
    // returns the same time, or throws: a rank throwing alone would leave the others waiting in
    // the next collective call
    double synchronisedTime_(auto& hierarchy, ResourcesManager& rm, Ions& ions) const
    {
        auto constexpr none = std::numeric_limits<double>::max();
        double earliest = none, latest = -none;
        for (int ilvl = 0; ilvl < hierarchy.getNumberOfLevels(); ++ilvl)
            for (auto& patch : *hierarchy.getPatchLevel(ilvl))
                for (auto const t : rm.getTimes(ions.chargeDensity(), *patch))
                {
                    earliest = std::min(earliest, t);
                    latest   = std::max(latest, t);
                }
        earliest = mpi::min(earliest);
        latest   = -mpi::min(-latest);

        if (earliest != none and latest - earliest > 1e-10 * std::max(1., std::abs(latest)))
            throw std::runtime_error(
                "ion fluxes and momentum tensor are computed only when all levels are at the "
                "same time: the ions are at times from "
                + std::to_string(earliest) + " to " + std::to_string(latest));
        return latest;
    }


    // the border sums of the population flux and momentum tensor. Their schedules are built for one
    // fill and discarded: a kept schedule holds its level, which would keep a level replaced at a
    // regrid alive, with its particles, until the next output on that level number
    void declareSums_(std::shared_ptr<ResourcesManager> const& rm, Ions& ions)
    {
        if (fluxSums_.size() == ions.size())
            return;

        for (auto const& pop : ions)
        {
            auto const& flux = pop.flux().name();
            fluxSums_.emplace_back(
                flux, flux, rm, nullptr,
                std::make_shared<TensorFieldGhostInterpOverlapFillPattern<GridLayout, 1>>());

            auto const& tensor = pop.momentumTensor().name();
            tensorSums_.emplace_back(
                tensor, tensor, rm, nullptr,
                std::make_shared<TensorFieldGhostInterpOverlapFillPattern<GridLayout, 2>>());
        }
    }


    // which level ghost particles each patch of a refined level deposits
    // (core::LevelGhostOwnership), rebuilt when SAMRAI replaces the level; made of boxes only, it
    // does not hold the level
    auto const& ownershipOf_(auto& hierarchy, int const levelNumber)
    {
        auto const level = hierarchy.getPatchLevel(levelNumber);
        auto& entry      = ownership_[levelNumber];
        if (entry.level.lock() == level and level != nullptr)
            return entry.ownership;

        entry = LevelOwnership{level, {}};
        if (levelNumber > 0)
        {
            auto const neighbors = makeLevelGhostNeighbors<GridLayout>(hierarchy, levelNumber);
            for (auto const& patch : *level)
                entry.ownership.emplace(patch->getGlobalId(),
                                        makeLevelGhostOwnershipFor<GridLayout>(*patch, neighbors));
        }
        return entry.ownership;
    }


    struct LevelOwnership
    {
        std::weak_ptr<SAMRAI::hier::PatchLevel> level;
        std::map<SAMRAI::hier::GlobalId, Ownership> ownership;
    };

    struct Done
    {
        std::weak_ptr<SAMRAI::hier::PatchLevel> level;
        double time = 0;
        What what;
    };

    std::vector<FluxSum> fluxSums_;
    std::vector<TensorSum> tensorSums_;
    std::map<int, LevelOwnership> ownership_;
    std::map<int, Done> done_;

    // scratch field registered and allocated by the hybrid messenger, shared with the VTK
    // writer: valid only inside fill. This relies on DiagnosticsManager::dump running every
    // compute before any write.
    Field density_{"PHARE_sumField", core::HybridQuantity::Scalar::rho};
};

} // namespace PHARE::amr

#endif
