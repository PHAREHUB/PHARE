#ifndef PHARE_FK_LEVEL_INITIALIZER_HPP
#define PHARE_FK_LEVEL_INITIALIZER_HPP

#include "core/errors.hpp"
#include "mpi/mpi_utils.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/ampere/ampere_maxwell.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/solvers/fk_electric.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/level_initializer/level_initializer.hpp"

#include "initializer/data_provider.hpp"

#include <cstddef>

namespace PHARE::solver
{
/**
 * @brief FKLevelInitializer initializes a level of a fully kinetic model.
 *
 * At t=0: particles and B (from the user functions on the root, refined from the coarser level on
 * a fine level: split particles, div-free B), E^0 = 0, then E moved half a level step back,
 *     E^{-dt_l/2} = E^0 - (dt_l/2) c^2 (curl B^0 - J_back),
 * with J_back the Esirkepov current of the backward half segment x^0 - (dt_l/2)v^0 -> x^0, so that
 * the leapfrog of SolverFK starts with the Gauss residual of t=0 of the level's own particles.
 * dt_l = dt / ratio^l, the fine time step of SolverFK::maxFinerLevelDt. The level ghost particles
 * deposit their J_back too (P4 plan D5), within the particle ghost width.
 *
 * A level (re)created at t_n > 0 (SAMRAI regrids it after the coarser level's step and
 * synchronization, P5a, P4 plan D6) is built from the old level where it overlaps it and from the
 * coarser level elsewhere: B by the div-free regrid, particles copied or split and the level ghost
 * particles split again (the hybrid regrid), E at its own half step t_n - dt_l/2, copied from the
 * old level, elsewhere and on the level ghosts interpolated in time at 1/4 between the coarser
 * E^{n-1/2} and the coarser E^{n+1/2} the coarser solver replayed just before
 * (ISolver::prepareFinerLevel), then relabelled (the allocation stamped t_n). Its Gauss gap
 * (splitting, interpolation, the old and new boundaries) is reduced by regridCleaningPasses Marder
 * passes, s = 1/2, against rho deposited at x - (dt_l/2) v by the new particles; not with
 * divergence_cleaning = 0 (explicit, measurement runs). A rebalanced root (no coarser level) is
 * copied and relabelled only.
 */
template<typename FKModel>
class FKLevelInitializer : public LevelInitializer<typename FKModel::amr_types>
{
    using amr_types                    = FKModel::amr_types;
    using level_ptr                    = amr_types::level_ptr;
    using hierarchy_ptr                = amr_types::hierarchy_ptr;
    using level_t                      = amr_types::level_t;
    using IPhysicalModelT              = IPhysicalModel<amr_types>;
    using Messenger                    = amr::HybridMessenger<FKModel>;
    using GridLayoutT                  = FKModel::gridlayout_type;
    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::options.interp_order;

    double dt_;

public:
    using model_type = FKModel;

    // Marder passes (s = 1/2) on the E of a level (re)created during the run (P5a); a test sets
    // it to 0 to measure the passes one by one (T14)
    static inline std::size_t regridCleaningPasses = 10;

    explicit FKLevelInitializer(PHARE::initializer::PHAREDict const& dict)
        : dt_{dict["time_step"]["value"].template to<double>()}
    {
    }

    void initialize(hierarchy_ptr const& hierarchy, int levelNumber, level_ptr const& oldLevel,
                    IPhysicalModelT& model, amr::IMessenger<IPhysicalModelT>& messenger,
                    double initDataTime, bool isRegridding) override
    {
        PHARE_LOG_SCOPE(1, "FKLevelInitializer::initialize");

        auto& fkModel     = static_cast<FKModel&>(model);
        auto& level       = amr_types::getLevel(*hierarchy, levelNumber);
        auto& fkMessenger = dynamic_cast<Messenger&>(messenger);
        auto& rm          = *fkModel.resourcesManager;
        auto& ions        = fkModel.state.ions;
        auto& E           = fkModel.state.electromag.E;
        auto& B           = fkModel.state.electromag.B;
        bool const isRoot = levelNumber == 0;
        auto const dt     = dt_ / level.getRatioToLevelZero().max();

        // SAMRAI also creates a level without an old one during the run (a level that comes
        // back): its E needs the coarser one too. A restart does not come here: SAMRAI reads the
        // levels and their data from the file (TimeRefinementIntegrator::initializeHierarchy);
        // with refinement boxes it may regrid a level later where the first run did not (SAMRAI
        // issue 199), a regrid like the others
        bool const duringRun = isRegridding or initDataTime > 0.;

        if (duringRun)
        {
            messenger.regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
            fkMessenger.regridElectric(E, hierarchy, levelNumber, oldLevel,
                                       initDataTime - 0.5 * dt);
            if (isRoot)
                fkMessenger.fillMagneticGhosts(B, level, initDataTime);
        }
        else if (isRoot)
        {
            model.initialize(level);
            messenger.fillRootGhosts(model, level, initDataTime);
            fkMessenger.fillMagneticGhosts(B, level, initDataTime);
        }
        else
            messenger.initLevel(hierarchy, model, level, initDataTime);

        // overlapping patches were each filled over their whole domain,
        // keep particles only on the patch owning their cell
        auto const neighbors = amr::makeSameLevelNeighbors(*hierarchy, levelNumber);
        for (auto& patch : rm.enumerate(level, ions))
            ions.eraseDomainParticlesIn(
                amr::makeNonOwnedBoxesFor<GridLayoutT::dimension>(*patch, neighbors));

        // moments at initDataTime, for diagnostics
        core::Interpolator<dimension, interp_order> interpolate;
        for (auto& patch : rm.enumerate(level, ions))
        {
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            core::resetMoments(ions);
            core::depositParticles(ions, layout, interpolate, core::DomainDeposit{});
        }
        fkMessenger.fillFluxBorders(ions, level, initDataTime);
        fkMessenger.fillDensityBorders(ions, level, initDataTime);
        for (auto& patch : rm.enumerate(level, ions))
        {
            // nodes next to the level boundary, as the hybrid level initializer does
            if (!isRoot)
            {
                auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
                core::depositParticles(ions, layout, interpolate, core::LevelGhostDeposit{});
            }
            ions.computeChargeDensity();
            ions.computeBulkVelocity();
            ions.computeChargeFlux();
        }
        fkMessenger.fillIonBorders(ions, level, initDataTime);

        if (duringRun)
            cleanRegriddedElectric_(*hierarchy, level, fkModel, fkMessenger, dt, initDataTime,
                                    isRoot);
        else
            startElectric_(level, fkModel, fkMessenger, dt, initDataTime, isRoot);

        fkMessenger.prepareStep(fkModel, level, initDataTime);
    }

private:
    // E^{-dt_l/2} from the current of the backward half segment; on a fine level E^0 = 0
    // replaces the coarse E refined by initLevel
    void startElectric_(level_t& level, FKModel& fkModel, Messenger& fkMessenger, double const dt,
                        double const initDataTime, bool const isRoot) const
    {
        auto& rm = *fkModel.resourcesManager;
        auto& E  = fkModel.state.electromag.E;
        auto& B  = fkModel.state.electromag.B;
        auto& J  = fkModel.state.J;

        if (!isRoot)
            for (auto& patch : rm.enumerate(level, E))
                E.zero();

        depositFKCurrent(level, fkModel, fkMessenger, dt, /*back=*/1., /*front=*/0., initDataTime);

        for (auto& patch : rm.enumerate(level, B, J, E))
        {
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            core::AmpereMaxwell{layout, fkModel.c2()}(B, J, E, -0.5 * dt);
        }

        // E is labelled with its own time, half a step before the particles. On a fine level
        // its level ghosts would be interpolated in time between coarse E values that do not
        // bracket t - dt_l/2 (the coarse E^{+dt_c/2} does not exist yet, P4 plan D6), so the
        // ghosts are left at zero, never read: SolverFK fills them at t - dt_l/2 before its
        // cleaning and Ampere update, once the coarse E^{+dt_c/2} exists (P4 plan D1).
        TimeSetter{level, fkModel, initDataTime - 0.5 * dt}(E);
        if (isRoot)
            fkMessenger.fillElectricGhosts(E, level, initDataTime - 0.5 * dt);
    }

    // the regridded E, at its own half step: relabelled, then on a fine level cleaned against
    // the charge of the new particles. Its level ghosts are those of the regrid (the coarser E
    // interpolated at that time), its same-level ghosts are filled before every pass. With
    // boundary_cleaning_mask the passes skip the level's edge nodes, which read those ghosts; the
    // edge of the region its finer level synchronized into the old level is not skipped
    void cleanRegriddedElectric_(SAMRAI::hier::PatchHierarchy const& hierarchy, level_t& level,
                                 FKModel& fkModel, Messenger& fkMessenger, double const dt,
                                 double const initDataTime, bool const isRoot) const
    {
        auto& E             = fkModel.state.electromag.E;
        auto const halfTime = initDataTime - 0.5 * dt;

        TimeSetter{level, fkModel, halfTime}(E);
        if (isRoot)
        {
            fkMessenger.fillElectricGhosts(E, level, halfTime);
            return;
        }
        // rho even without a pass: T14 measures the Gauss residual before the first one
        if (fkModel.divergence_cleaning > 0.)
        {
            depositFKCurrent(level, fkModel, fkMessenger, dt, /*back=*/1., /*front=*/0.,
                             initDataTime);
            auto const skipped = fkModel.boundary_cleaning_mask
                                     ? cleaningSkippedNodes<GridLayoutT>(
                                           hierarchy, level.getLevelNumber(), {})
                                     : PatchBoxes<dimension>{};
            cleanElectric(level, fkModel, fkMessenger, E, dt, 0.5, regridCleaningPasses, halfTime,
                          skipped);
        }
    }
};

} // namespace PHARE::solver

#endif
