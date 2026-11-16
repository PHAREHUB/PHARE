#ifndef PHARE_FK_LEVEL_INITIALIZER_HPP
#define PHARE_FK_LEVEL_INITIALIZER_HPP

#include "core/errors.hpp"
#include "mpi/mpi_utils.hpp"
#include "core/numerics/moments/moments.hpp"
#include "core/numerics/esirkepov/esirkepov.hpp"
#include "core/numerics/ampere/ampere_maxwell.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/level_initializer/level_initializer.hpp"

#include "initializer/data_provider.hpp"

#include <stdexcept>

namespace PHARE::solver
{
/**
 * @brief FKLevelInitializer initializes the root level of a fully kinetic model: particles, B,
 * E^0 = 0, then moves E half a step back,
 *     E^{-1/2} = E^0 - (dt/2) c^2 (curl B^0 - J_back),
 * with J_back the Esirkepov current of the backward half segment x^0 - (dt/2)v^0 -> x^0, so that
 * the leapfrog of SolverFK starts with the Gauss residual of t=0.
 */
template<typename FKModel>
class FKLevelInitializer : public LevelInitializer<typename FKModel::amr_types>
{
    using amr_types                    = FKModel::amr_types;
    using level_ptr                    = amr_types::level_ptr;
    using hierarchy_ptr                = amr_types::hierarchy_ptr;
    using IPhysicalModelT              = IPhysicalModel<amr_types>;
    using Messenger                    = amr::HybridMessenger<FKModel>;
    using GridLayoutT                  = FKModel::gridlayout_type;
    static constexpr auto dimension    = GridLayoutT::dimension;
    static constexpr auto interp_order = GridLayoutT::options.interp_order;

    double c2_;
    double dt_;

public:
    using model_type = FKModel;

    explicit FKLevelInitializer(PHARE::initializer::PHAREDict const& dict)
        : c2_{dict["c_over_va"].template to<double>() * dict["c_over_va"].template to<double>()}
        , dt_{dict["time_step"]["value"].template to<double>()}
    {
    }

    void initialize(hierarchy_ptr const& hierarchy, int levelNumber, level_ptr const& /*oldLevel*/,
                    IPhysicalModelT& model, amr::IMessenger<IPhysicalModelT>& messenger,
                    double initDataTime, bool isRegridding) override
    {
        PHARE_LOG_SCOPE(1, "FKLevelInitializer::initialize");

        if (isRegridding or levelNumber != 0)
            throw std::runtime_error("FKLevelInitializer: only the root level is supported yet");

        auto& fkModel     = static_cast<FKModel&>(model);
        auto& level       = amr_types::getLevel(*hierarchy, levelNumber);
        auto& fkMessenger = dynamic_cast<Messenger&>(messenger);
        auto& rm          = *fkModel.resourcesManager;
        auto& ions        = fkModel.state.ions;
        auto& E           = fkModel.state.electromag.E;
        auto& B           = fkModel.state.electromag.B;
        auto& J           = fkModel.state.J;

        model.initialize(level);
        messenger.fillRootGhosts(model, level, initDataTime);
        fkMessenger.fillMagneticGhosts(B, level, initDataTime);

        // overlapping patches were each filled over their whole domain,
        // keep particles only on the patch owning their cell
        auto const neighbors = amr::makeSameLevelNeighbors(*hierarchy, levelNumber);
        for (auto& patch : rm.enumerate(level, ions))
            ions.eraseDomainParticlesIn(
                amr::makeNonOwnedBoxesFor<GridLayoutT::dimension>(*patch, neighbors));

        // moments at t=0, for diagnostics
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
            ions.computeChargeDensity();
            ions.computeBulkVelocity();
            ions.computeChargeFlux();
        }
        fkMessenger.fillIonBorders(ions, level, initDataTime);

        // E^{-1/2} from the current of the backward half segment
        core::Esirkepov<dimension, interp_order> esirkepov;
        try
        {
            for (auto& patch : rm.enumerate(level, ions, J))
            {
                auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
                J.zero();
                for (auto const& pop : ions)
                {
                    try
                    {
                        esirkepov(pop.domainParticles(), J, layout, pop.charge(), dt_,
                                  /*back=*/1., /*front=*/0.);
                    }
                    catch (core::DictionaryException& ex)
                    {
                        throw ex("population", pop.name());
                    }
                }
            }
        }
        catch (core::DictionaryException const& ex)
        {
            PHARE_LOG_ERROR(ex());
        }
        if (mpi::any_errors())
            throw core::DictionaryException{}("ID", "FKLevelInitializer::initialize");

        fkMessenger.fillCurrentBorders(J, level, initDataTime);

        for (auto& patch : rm.enumerate(level, B, J, E))
        {
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            core::AmpereMaxwell{layout, c2_}(B, J, E, -0.5 * dt_);
        }

        TimeSetter setTime{level, fkModel, initDataTime};
        setTime(E, J);
        fkMessenger.fillElectricGhosts(E, level, initDataTime);

        fkMessenger.prepareStep(fkModel, level, initDataTime);
    }
};

} // namespace PHARE::solver

#endif
