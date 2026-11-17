#ifndef PHARE_SOLVER_FK_HPP
#define PHARE_SOLVER_FK_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/errors.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/numerics/ampere/ampere_maxwell.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/numerics/marder/marder.hpp"
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

#include <cmath>
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
 *     particles and summed on patch borders; with divergence cleaning, also rho^{n-1/2} at
 *     x^n - (dt/2)v^n
 *  2. with divergence cleaning, E^{n-1/2} += d dt grad(div E^{n-1/2} - c^2 rho^{n-1/2}) (Marder),
 *     then E^{n+1/2} = E^{n-1/2} + dt c^2 (curl B^n - J^n)
 *  3. B^{n+1/2} = B^n - (dt/2) curl E^{n+1/2}
 *  4. Boris push with E^{n+1/2} and B^{n+1/2}
 *  5. B^{n+1} = B^n - dt curl E^{n+1/2}
 *  6. on a level with a finer one, E one step ahead: steps 1 and 2 replayed from the particles
 *     and B at n+1 on a copy of E^{n+1/2}, E^{n+3/2}, the newer end of the time interpolation of
 *     the finer level E ghosts past E^{n+1/2} (P4 plan D1)
 *
 * The model E is labelled with its own time t_{n+1/2}, J with the time of the particles it
 * comes from, t_n, rho with its own, t_n - dt/2. The deposit only depends on the
 * particles at time n, so anything that changes particles between two steps is accounted for.
 * The cleaning reads one ghost layer of E^{n-1/2}, filled at the end of the previous step's
 * Ampere update (or by the level initializer on the root), and writes its physical nodes; Ampere
 * only reads E on the nodes it writes, and its ghost fill covers both. On a refined level the E
 * ghosts are filled again at t_n - dt/2 before step 2: the level initializer cannot fill them
 * (no coarse E bracket at t=0), and at the first substep of a coarse step the coarse E has been
 * synchronized since. Charge conservation keeps the Gauss residual between two cleanings: each
 * cleaning multiplies its mode k by 1 - d dt sum_i 4/dx_i^2 sin^2(k_i dx_i / 2). The level
 * initializer does not clean: the residual of E^{-1/2} (E starts at zero) is what the first step
 * cleans.
 *
 * Synchronization with a finer level (P4 plan D4): the messenger coarsens the finer level's time
 * average of E (the flux sum, at t_{n+1/2}) onto the model E over the covered region, then
 * reflux recomputes B^{n+1} = B^n - dt curl E^{n+1/2} on the whole level from a copy of B^n:
 * bitwise step 5 where E did not change. The finer model E, at another time, is not coarsened.
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

    // E^{n+3/2} on a level with a finer one (step 6)
    VecFieldT Eahead_{this->name() + "_Eahead", core::HybridQuantity::Vector::E};

    // B^n and its time, per level: the start of the reflux Faraday
    VecFieldT Bold_{this->name() + "_Bold", core::HybridQuantity::Vector::B};
    std::unordered_map<int, double> oldTime_;

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

        // with refinement: the finer level's flux sum is coarsened onto the model E, at the same
        // time t_{n+1/2} (the model fills its info first, messenger_registration.hpp)
        fkInfo.refluxElectric  = fkInfo.modelElectric;
        fkInfo.fluxSumElectric = fluxSumE_.name();
        fkInfo.aheadElectric   = Eahead_.name();
    }

    void registerResources(IPhysicalModel_t& model) override
    {
        auto& rm = *dynamic_cast<FKModel&>(model).resourcesManager;
        rm.registerResources(electromagHalf_);
        rm.registerResources(fluxSumE_);
        rm.registerResources(Eahead_);
        rm.registerResources(Bold_);
    }

    void allocate(IPhysicalModel_t& model, patch_t& patch, double const allocateTime) const override
    {
        auto& rm = *dynamic_cast<FKModel&>(model).resourcesManager;
        rm.allocate(electromagHalf_, patch, allocateTime);
        rm.allocate(fluxSumE_, patch, allocateTime);
        rm.allocate(Eahead_, patch, allocateTime);
        rm.allocate(Bold_, patch, allocateTime);
    }

    // B^n for the reflux of a level with a finer one (the step itself keeps B^n in the model B
    // until its last Faraday)
    void prepareStep(IPhysicalModel_t& model, level_t& level, double const currentTime) override
    {
        auto& fkModel = dynamic_cast<FKModel&>(model);
        auto& B       = fkModel.state.electromag.B;
        for (auto& patch : fkModel.resourcesManager->enumerate(level, B, Bold_))
            Bold_.copyData(B);
        oldTime_[level.getLevelNumber()] = currentTime;
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

    // the messenger has coarsened the finer level's flux sum onto the model E^{n+1/2}
    void reflux(IPhysicalModel_t& model, level_t& level, IMessenger& messenger,
                double const time) override
    {
        auto& fkModel     = dynamic_cast<FKModel&>(model);
        auto& fkMessenger = dynamic_cast<Messenger&>(messenger);
        auto& B           = fkModel.state.electromag.B;
        auto const dt     = time - oldTime_.at(level.getLevelNumber());

        Faraday_t{level, fkModel}(Bold_, fkModel.state.electromag.E, B, dt);
        TimeSetter{level, fkModel, time}(B);
        fkMessenger.fillMagneticGhosts(B, level, time);
    }

    double computeStableDt(IPhysicalModel_t& /*model*/, level_t& /*level*/,
                           CFLNumbers const& /*cflNumbers*/) override
    {
        throw std::runtime_error("SolverFK: adaptive time stepping is not supported");
    }

    // the light CFL c dt / dx is the same on every level with dt / ratio
    NO_DISCARD double maxFinerLevelDt(double const coarseDt, int const ratio) const override
    {
        return coarseDt / ratio;
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
        return std::forward_as_tuple(electromagHalf_, fluxSumE_, Eahead_, Bold_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(electromagHalf_, fluxSumE_, Eahead_, Bold_);
    }


private:
    void advanceElectric_(level_t& level, FKModel& model, Messenger& messenger, VecFieldT& E,
                          double const dt, double const time);

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

    auto& model     = dynamic_cast<FKModel&>(level_model);
    auto& messenger = dynamic_cast<Messenger&>(fromCoarserMessenger);
    auto& level     = boxing.setup_level(hierarchy, levelNumber);
    auto& rm        = *model.resourcesManager;
    auto& E         = model.state.electromag.E;
    auto& B         = model.state.electromag.B;
    auto const dt   = newTime - currentTime;

    // E^{n-1/2} was built with this dt: by the level initializer from the input time step, or by
    // the previous step. A different dt would silently break the half step staggering.
    for (auto& patch : level)
        for (auto const time : rm.getTimes(E, *patch))
            if (std::abs(time - (currentTime - 0.5 * dt)) > 1e-6 * dt)
                throw std::runtime_error("SolverFK: E is not at t - dt/2, the time step changed");

    auto const halfTime = 0.5 * (currentTime + newTime);

    // E^{n-1/2} level ghosts, from the coarser level's E bracketing t_n - dt/2
    if (levelNumber > 0)
        messenger.fillElectricGhosts(E, level, currentTime - 0.5 * dt);

    // 1. and 2.
    advanceElectric_(level, model, messenger, E, dt, currentTime);
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

    // 6. E^{n+3/2}, before the finer level's substeps (P4 plan D1). On the root it is the next
    // step's E^{n+3/2} where the synchronization does not change E^{n+1/2} and B^{n+1}
    if (hierarchy.finerLevelExists(levelNumber))
    {
        for (auto& patch : rm.enumerate(level, E, Eahead_))
            Eahead_.copyData(E); // with its ghosts, read by the cleaning
        advanceElectric_(level, model, messenger, Eahead_, dt, newTime);
        TimeSetter{level, model, newTime + 0.5 * dt}(Eahead_);
    }
}



// 1. J and rho from the particles at time, then 2. E (at time - dt/2) cleaned, d from
// s = d dt sum_i 2/dx_i^2 on this level, and advanced to time + dt/2
template<typename FKModel, typename AMR_Types>
void SolverFK<FKModel, AMR_Types>::advanceElectric_(level_t& level, FKModel& model,
                                                    Messenger& messenger, VecFieldT& E,
                                                    double const dt, double const time)
{
    auto& rm  = *model.resourcesManager;
    auto& B   = model.state.electromag.B;
    auto& J   = model.state.J;
    auto& rho = model.state.rho;

    depositFKCurrent(level, model, messenger, dt, /*back=*/1., /*front=*/1., time);

    if (model.divergence_cleaning > 0.)
        for (auto& patch : rm.enumerate(level, rho, E))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            double sum        = 0.;
            for (auto const dl : layout.meshSize())
                sum += 2. / (dl * dl);
            auto const d = model.divergence_cleaning / (dt * sum);
            core::Marder{layout, model.c2(), d}(E, rho, dt);
        }
    for (auto& patch : rm.enumerate(level, B, J, E))
    {
        auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
        core::AmpereMaxwell{layout, model.c2()}(B, J, E, dt);
    }
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

    // the moments at n+1 are diagnostics only (and coarsened at synchronization), the fields
    // only use the deposit of step 1. On a refined level the level ghost particles complete them
    messenger.fillFluxBorders(ions, level, newTime);
    messenger.fillDensityBorders(ions, level, newTime);
    if (level.getLevelNumber() > 0)
        messenger.fillIonPopMomentGhosts(ions, level, newTime);

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
