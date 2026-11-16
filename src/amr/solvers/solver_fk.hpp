#ifndef PHARE_SOLVER_FK_HPP
#define PHARE_SOLVER_FK_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/errors.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/numerics/ampere/ampere_maxwell.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/utilities/algorithm.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/solvers/solver.hpp"
#include "amr/solvers/particle_level_boxings.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/solvers/fk_current_deposit.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"

#include <SAMRAI/hier/Patch.h>
#include "SAMRAI/hier/PatchLevel.h"

#include <memory>
#include <string>
#include <stdexcept>
#include <unordered_map>


namespace PHARE::solver
{
/**
 * @brief SolverFK advances a fully kinetic model with an explicit leapfrog: particles (x, v) and
 * B are known at integer times, E at half integer times. One step from n to n+1 is
 *
 *  1. J^n by Esirkepov over x^n - (dt/2)v^n -> x^n + (dt/2)v^n, deposited by owned domain
 *     particles and summed on patch borders
 *  2. E^{n+1/2} = E^{n-1/2} + dt c^2 (curl B^n - J^n)
 *  3. B^{n+1/2} = B^n - (dt/2) curl E^{n+1/2}
 *  4. Boris push with E^{n+1/2} and B^{n+1/2}
 *  5. B^{n+1} = B^n - dt curl E^{n+1/2}
 *
 * The model E is labelled with its own time t_{n+1/2}, J with the time of the particles it
 * comes from, t_n. The deposit only depends on the
 * particles at time n, so anything that changes particles between two steps is accounted for.
 * Single level only for now.
 */
template<typename FKModel, typename AMR_Types>
class SolverFK : public ISolver<AMR_Types>
{
private:
    static constexpr auto dimension    = FKModel::dimension;
    static constexpr auto interp_order = FKModel::interp_order;

    using Electromag       = FKModel::electromag_type;
    using Ions             = FKModel::ions_type;
    using VecFieldT        = FKModel::vecfield_type;
    using GridLayout       = FKModel::gridlayout_type;
    using IPhysicalModel_t = IPhysicalModel<AMR_Types>;
    using IMessenger       = amr::IMessenger<IPhysicalModel_t>;
    using Messenger        = amr::HybridMessenger<FKModel>;

    using Faraday_t    = FaradayLevelTransformer<FKModel>;
    using IonUpdater_t = core::IonUpdater<Ions, Electromag, GridLayout>;

    // E^{n+1/2} and B^{n+1/2}, the fields the particles are pushed with
    Electromag electromagHalf_{"EMHalf"};

    VecFieldT fluxSumE_{this->name() + "_fluxSumE", core::HybridQuantity::Vector::E};

    IonUpdater_t ionUpdater_;


public:
    using patch_t     = AMR_Types::patch_t;
    using level_t     = AMR_Types::level_t;
    using hierarchy_t = AMR_Types::hierarchy_t;
    using level_ptr   = AMR_Types::level_ptr;


    explicit SolverFK(PHARE::initializer::PHAREDict const& dict)
        : ISolver<AMR_Types>{"FK"}
        , ionUpdater_{dict["algo"]["ion_updater"]}
    {
    }

    ~SolverFK() override = default;


    std::string modelName() const override { return FKModel::model_name; }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override
    {
        auto& fkInfo = dynamic_cast<amr::HybridMessengerInfo&>(*info);

        fkInfo.ghostElectric.emplace_back(electromagHalf_.E.name());
        fkInfo.ghostMagnetic.emplace_back(electromagHalf_.B.name());
        fkInfo.initMagnetic.emplace_back(electromagHalf_.B.name());

        // only used with refinement: reflux writes a temporary, never the model E
        fkInfo.refluxElectric  = electromagHalf_.E.name();
        fkInfo.fluxSumElectric = fluxSumE_.name();
    }

    void registerResources(IPhysicalModel_t& model) override
    {
        auto& rm = *dynamic_cast<FKModel&>(model).resourcesManager;
        rm.registerResources(electromagHalf_);
        rm.registerResources(fluxSumE_);
    }

    void allocate(IPhysicalModel_t& model, patch_t& patch, double const allocateTime) const override
    {
        auto& rm = *dynamic_cast<FKModel&>(model).resourcesManager;
        rm.allocate(electromagHalf_, patch, allocateTime);
        rm.allocate(fluxSumE_, patch, allocateTime);
    }

    // B^n stays in the model B until the last Faraday of the step, no copy needed
    void prepareStep(IPhysicalModel_t& /*model*/, level_t& /*level*/,
                     double const /*currentTime*/) override
    {
    }

    void accumulateFluxSum(IPhysicalModel_t& model, level_t& level, double const coef) override
    {
        auto& fkModel = dynamic_cast<FKModel&>(model);
        auto& E       = fkModel.state.electromag.E;

        for (auto& patch : fkModel.resourcesManager->enumerate(level, fluxSumE_, E))
            core::operate<core::PlusEqualsProduct>(fluxSumE_, E, coef);
    }

    void resetFluxSum(IPhysicalModel_t& model, level_t& level) override
    {
        auto& fkModel = dynamic_cast<FKModel&>(model);
        for (auto& patch : fkModel.resourcesManager->enumerate(level, fluxSumE_))
            fluxSumE_.zero();
    }

    void reflux(IPhysicalModel_t& /*model*/, level_t& /*level*/, IMessenger& /*messenger*/,
                double const /*time*/) override
    {
        throw std::runtime_error("SolverFK: refinement is not supported yet");
    }

    double computeStableDt(IPhysicalModel_t& /*model*/, level_t& /*level*/,
                           CFLNumbers const& /*cflNumbers*/) override
    {
        throw std::runtime_error("SolverFK: adaptive time stepping is not supported");
    }

    void advanceLevel(hierarchy_t const& hierarchy, int const levelNumber, IPhysicalModel_t& model,
                      IMessenger& fromCoarserMessenger, double const currentTime,
                      double const newTime) override;

    void onRegrid() override
    {
        boxing.clear();
        ionUpdater_.reset();
    }


    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(electromagHalf_, fluxSumE_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(electromagHalf_, fluxSumE_);
    }


private:
    void moveParticles_(level_t& level, FKModel& model, Messenger& messenger,
                        double const currentTime, double const newTime);


    ParticleLevelBoxings<IonUpdater_t, GridLayout, AMR_Types> boxing;
};



template<typename FKModel, typename AMR_Types>
void SolverFK<FKModel, AMR_Types>::advanceLevel(hierarchy_t const& hierarchy, int const levelNumber,
                                                IPhysicalModel_t& level_model,
                                                IMessenger& fromCoarserMessenger,
                                                double const currentTime, double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverFK::advanceLevel");

    if (levelNumber != 0)
        throw std::runtime_error("SolverFK: only the root level is supported yet");

    auto& model     = dynamic_cast<FKModel&>(level_model);
    auto& messenger = dynamic_cast<Messenger&>(fromCoarserMessenger);
    auto& level     = boxing.setup_level(hierarchy, levelNumber);
    auto& rm        = *model.resourcesManager;
    auto& E         = model.state.electromag.E;
    auto& B         = model.state.electromag.B;
    auto const dt   = newTime - currentTime;

    auto const halfTime = 0.5 * (currentTime + newTime);

    // 1. J^n from the particles at time n
    depositFKCurrent(level, model, messenger, dt, /*back=*/1., /*front=*/1., currentTime);

    // 2. E^{n-1/2} -> E^{n+1/2}
    for (auto& patch : rm.enumerate(level, B, model.state.J, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        core::AmpereMaxwell{layout, model.c2()}(B, model.state.J, E, dt);
    }
    TimeSetter{level, model, halfTime}(E);
    messenger.fillElectricGhosts(E, level, halfTime);

    // 3. B^{n+1/2}, and the E^{n+1/2} copy the particles are pushed with
    auto const halfDt = 0.5 * dt;
    Faraday_t{level, model}(B, E, electromagHalf_.B, halfDt);
    for (auto& patch : rm.enumerate(level, E, electromagHalf_))
        electromagHalf_.E.copyData(E);
    TimeSetter{level, model, halfTime}(electromagHalf_.B, electromagHalf_.E);
    messenger.fillMagneticGhosts(electromagHalf_.B, level, halfTime);

    // 4. particles from n to n+1
    moveParticles_(level, model, messenger, currentTime, newTime);

    // 5. B^n -> B^{n+1}
    Faraday_t{level, model}(B, E, B, dt);
    TimeSetter{level, model, newTime}(B);
    messenger.fillMagneticGhosts(B, level, newTime);
}



template<typename FKModel, typename AMR_Types>
void SolverFK<FKModel, AMR_Types>::moveParticles_(level_t& level, FKModel& model,
                                                  Messenger& messenger, double const currentTime,
                                                  double const newTime)
{
    PHARE_LOG_SCOPE(1, "SolverFK::moveParticles_");

    TimeSetter setTime{level, model, newTime};
    auto& rm                = *model.resourcesManager;
    auto const& levelBoxing = boxing[level.getLevelNumber()].patches;
    auto& ions              = model.state.ions;

    try
    {
        auto const dt = newTime - currentTime;
        for (auto& patch : rm.enumerate(level, ions, electromagHalf_))
            ionUpdater_.updatePopulations(ions, electromagHalf_,
                                          levelBoxing.at(amr::to_string(patch->getGlobalId())), dt,
                                          core::UpdaterMode::all);
    }
    catch (core::DictionaryException const& ex)
    {
        PHARE_LOG_ERROR(ex());
    }
    if (mpi::any_errors())
        throw core::DictionaryException{}("ID", "SolverFK::moveParticles_");

    // this needs to be done before calling the messenger
    setTime(ions);

    // the moments at n+1 are diagnostics only, the fields only use the deposit of step 1.
    // No level ghost moments (fillIonPopMomentGhosts): single level for now
    messenger.fillFluxBorders(ions, level, newTime);
    messenger.fillDensityBorders(ions, level, newTime);

    // the next deposit relies on domain particles being owned ones
    PHARE_DEBUG_DO({
        for (auto& patch : rm.enumerate(level, ions))
        {
            auto const& boxes = levelBoxing.at(amr::to_string(patch->getGlobalId()));
            for (auto const& pop : ions)
                for (auto const& particle : pop.domainParticles())
                    if (!boxes.isOwned(particle.iCell))
                        throw std::runtime_error("domain particle outside owned cells");
        }
    })

    messenger.fillIonGhostParticles(ions, level, newTime);

    // the exchange hands a leaving particle to every patch whose domain contains its
    // cell; on overlapping patches only the owner keeps it
    for (auto& patch : rm.enumerate(level, ions))
        ions.eraseDomainParticlesIn(
            levelBoxing.at(amr::to_string(patch->getGlobalId())).nonOwnedBoxes);

    for (auto& patch : rm.enumerate(level, ions))
        ionUpdater_.updateIons(ions);

    messenger.fillIonBorders(ions, level, newTime);
}



} // namespace PHARE::solver


#endif
