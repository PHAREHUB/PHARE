#ifndef PHARE_SOLVER_FK_HPP
#define PHARE_SOLVER_FK_HPP

#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/errors.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/numerics/ion_updater/ion_updater.hpp"
#include "core/utilities/algorithm.hpp"
#include "mpi/mpi_utils.hpp"

#include "amr/samrai.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/solvers/particle_level_boxings.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/solvers/solver_field_evolvers.hpp"
#include "amr/solvers/fk_electric.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"

#include <SAMRAI/hier/Patch.h>
#include "SAMRAI/hier/PatchLevel.h"
#include <SAMRAI/tbox/Database.h>
#include <SAMRAI/tbox/Serializable.h>

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
 *
 * Partial reflux after each finer substep but the last (P4 T9 residual, option B-x,
 * docs/superpowers/specs/p4/t9_residual.md): the next substep's B ghosts are filled from the
 * coarse B at its start, B^n - (k+1) dt_f curl E^, where E^ is the coarse E^{n+1/2} outside the
 * finer level and, on the covered edges, the mean of the finer E of the substeps done,
 * extrapolated to the coarse E time with the last two (t_n + dt/2: 1.5 E_f^(0) - 0.5 E_f^(-1) at
 * ratio 2). Interpolated in time towards the predicted coarse B^{n+1} instead, the ghost read a
 * value the reflux then replaced, and the interface grew. Outside the finer level the result is
 * the former interpolation; the half step ghosts of the push keep the interpolation.
 *
 * With boundary_cleaning_mask, the cleanings of steps 2 and 6 skip (F = 0) the nodes on the edge of
 * a refined level and on the edge of the cells the finer level covered at the last
 * synchronization of E^{n-1/2} (fk_electric.hpp, regionEdgeNodes). That region is the finer
 * level's at the end of the previous step, kept across a regrid, which comes after the
 * synchronization: none on a level whose finer level is new, the old one for one more step on a
 * level whose finer level was removed. The skipped nodes are built when the level, or the kept
 * region, changes. Without the option nothing is built and nothing skipped.
 *
 * On a refined level, after each fill of E^{n-1/2} (step 1) and E^{n+1/2} (after step 2), the
 * first faces of E outside the level take Gauss's law on its edge nodes, with the charge at E's
 * time (fk_electric.hpp, gaussOnLevelBorder): the fills give them the coarser E, while the charge
 * of those nodes follows the fine particles (BORDER-GAUSS, docs/superpowers/specs/p4/
 * border_gauss.md). No state: the next fill replaces them, the charge sets them again. A band
 * refresh (lastStep) moves the charge of the edge nodes until the next step's fill.
 */
template<typename FKModel, typename AMR_Types>
class SolverFK : public ISolver<AMR_Types>, public SAMRAI::tbox::Serializable
{
private:
    static constexpr auto dimension    = FKModel::dimension;
    static constexpr auto interp_order = FKModel::interp_order;

    using Electromag       = FKModel::electromag_type;
    using Ions             = FKModel::ions_type;
    using VecFieldT        = FKModel::vecfield_type;
    using FieldT           = FKModel::field_type;
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

    // partial reflux (B-x). Finer level: its E at the start of the substep, then the covered E
    // coarsened. Coarser level: its E with the covered edges from the finer level, and the B its
    // Faraday gives at the end of the finer substep, read by the finer level's next B ghosts
    VecFieldT Epartial_{this->name() + "_Epartial", core::HybridQuantity::Vector::E};
    VecFieldT Ehat_{this->name() + "_Ehat", core::HybridQuantity::Vector::E};
    VecFieldT Bmid_{this->name() + "_Bmid", core::HybridQuantity::Vector::B};

    // the model rho of step 1 (the start of the segments, which the diagnostics dump), kept while
    // the charge at their end sets the faces of E^{n+1/2} outside a refined level
    FieldT rhoStart_{this->name() + "_rhoStart", core::HybridQuantity::Scalar::rho};
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
        auto& restarts = *SamraiLifeCycle::getRestartManager();
        if (restarts.isFromRestart())
            if (auto const root = restarts.getRootDatabase(); root->keyExists(restartName_))
                coveredFromRestart_(*root->getDatabase(restartName_));
        restarts.registerRestartItem(restartName_, this);
    }

    ~SolverFK() override
    {
        SamraiLifeCycle::getRestartManager()->unregisterRestartItem(restartName_);
    }

    // in the restart file, the cells each level's finer level synchronized at the end of the
    // last step (boundary_cleaning_mask), which a regrid at the end of that step has replaced
    // by the time the file is written. Every other record of the solver is rewritten by the
    // first step before it is read
    void putToRestart(std::shared_ptr<SAMRAI::tbox::Database> const& db) const override
    {
        for (auto const& [levelNumber, mask] : cleaningMasks_)
        {
            std::vector<int> bounds;
            for (auto const& box : mask.covered)
                for (auto const& point : {box.lower, box.upper})
                    for (std::size_t d = 0; d < dimension; ++d)
                        bounds.push_back(point[d]);
            if (!bounds.empty()) // SAMRAI errors on writing 0 size arrays
                db->putIntegerVector(coveredKey_ + std::to_string(levelNumber), bounds);
        }
    }


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

        fkInfo.partialFluxElectric = Epartial_.name();
        fkInfo.partialElectric     = Ehat_.name();
        fkInfo.midMagnetic         = Bmid_.name();
    }

    void registerResources(IPhysicalModel_t& model) override
    {
        auto& rm = *dynamic_cast<FKModel&>(model).resourcesManager;
        rm.registerResources(electromagHalf_);
        rm.registerResources(fluxSumE_);
        rm.registerResources(Eahead_);
        rm.registerResources(Bold_);
        rm.registerResources(Epartial_);
        rm.registerResources(Ehat_);
        rm.registerResources(Bmid_);
        rm.registerResources(rhoStart_);
    }

    void allocate(IPhysicalModel_t& model, patch_t& patch, double const allocateTime) const override
    {
        auto& rm = *dynamic_cast<FKModel&>(model).resourcesManager;
        rm.allocate(electromagHalf_, patch, allocateTime);
        rm.allocate(fluxSumE_, patch, allocateTime);
        rm.allocate(Eahead_, patch, allocateTime);
        rm.allocate(Bold_, patch, allocateTime);
        rm.allocate(Epartial_, patch, allocateTime);
        rm.allocate(Ehat_, patch, allocateTime);
        rm.allocate(Bmid_, patch, allocateTime);
        rm.allocate(rhoStart_, patch, allocateTime);
    }

    // B^n for the reflux of a level with a finer one (the step itself keeps B^n in the model B
    // until its last Faraday); on a refined level, E^{n-1/2} for the trend of the partial reflux
    void prepareStep(IPhysicalModel_t& model, level_t& level, double const currentTime) override
    {
        auto& fkModel = dynamic_cast<FKModel&>(model);
        auto& B       = fkModel.state.electromag.B;
        auto& E       = fkModel.state.electromag.E;
        for (auto& patch : fkModel.resourcesManager->enumerate(level, B, Bold_))
            Bold_.copyData(B);
        if (level.getLevelNumber() > 0)
            for (auto& patch : fkModel.resourcesManager->enumerate(level, E, Epartial_))
                Epartial_.copyData(E);
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

    void partialReflux(IPhysicalModel_t& model, level_t& coarseLevel, level_t& fineLevel,
                       IMessenger& messenger, std::size_t const done, std::size_t const ratio,
                       double const midTime) override;

    // a finer level (re)created at time, after this level's step and synchronization (P5a): step
    // 6 replayed from the synchronized B^n and the particles at n, E^{n+1/2} one step ahead of
    // the model E^{n-1/2}. The Eahead_ of the step, if any, predates the synchronization, and a
    // level that had no finer one has none. With boundary_cleaning_mask the cleaning skips the
    // level's nodes for the region the finer level synchronized in this step (keepCoveredRegion_)
    void prepareFinerLevel(hierarchy_t const& hierarchy, IPhysicalModel_t& model, level_t& level,
                           IMessenger& messenger, double const time) override
    {
        auto& fkModel     = dynamic_cast<FKModel&>(model);
        auto& fkMessenger = dynamic_cast<Messenger&>(messenger);
        auto& E           = fkModel.state.electromag.E;
        auto const dt     = time - oldTime_.at(level.getLevelNumber());

        for (auto& patch : fkModel.resourcesManager->enumerate(level, E, Eahead_))
            Eahead_.copyData(E); // with its ghosts, read by the cleaning
        advanceFKElectric(level, fkModel, fkMessenger, Eahead_, dt, time,
                          skippedNodes_(hierarchy, level.getLevelNumber(), fkModel),
                          borderNodes_(hierarchy, level.getLevelNumber()));
        TimeSetter{level, fkModel, time + 0.5 * dt}(Eahead_);
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
        return std::forward_as_tuple(electromagHalf_, fluxSumE_, Eahead_, Bold_, Epartial_, Ehat_,
                                     Bmid_, rhoStart_);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(electromagHalf_, fluxSumE_, Eahead_, Bold_, Epartial_, Ehat_,
                                     Bmid_, rhoStart_);
    }


private:
    void moveParticles_(level_t& level, FKModel& model, Messenger& messenger,
                        double const currentTime, double const newTime);

    PatchBoxes<dimension> const& skippedNodes_(hierarchy_t const& hierarchy, int const levelNumber,
                                               FKModel const& model);

    // the level border nodes of a level (levelBorderNodesOf), built once per level object
    PatchBorderNodes<dimension> const& borderNodes_(hierarchy_t const& hierarchy,
                                                    int const levelNumber)
    {
        auto const level      = hierarchy.getPatchLevel(levelNumber);
        auto& [object, nodes] = borderNodeCache_[levelNumber];
        if (object.lock() != level)
        {
            object = level;
            nodes  = levelBorderNodesOf<GridLayout>(hierarchy, levelNumber);
        }
        return nodes;
    }
    std::unordered_map<int, std::pair<std::weak_ptr<level_t>, PatchBorderNodes<dimension>>>
        borderNodeCache_;
    void keepCoveredRegion_(hierarchy_t const& hierarchy, int const levelNumber,
                            FKModel const& model);


    ParticleLevelBoxings<IonUpdater_t, GridLayout, AMR_Types> boxing;

    // the nodes the cleaning skips on a level (boundary_cleaning_mask): those of the level object
    // they were built for, from the cells the finer level covered (`covered`, in level cells, not
    // per patch) at the end of the last step, taken from the level `finer`
    struct CleaningMask
    {
        std::weak_ptr<level_t> level;
        std::weak_ptr<level_t> finer;
        std::vector<core::Box<int, dimension>> covered;
        PatchBoxes<dimension> skipped;
        bool stale = true;
    };
    std::unordered_map<int, CleaningMask> cleaningMasks_;

    // the masks of the restart file: their covered cells, the level objects are new (stale)
    void coveredFromRestart_(SAMRAI::tbox::Database& db)
    {
        for (auto const& key : db.getAllKeys())
            if (key.starts_with(coveredKey_))
            {
                auto const bounds = db.getIntegerVector(key);
                auto& covered = cleaningMasks_[std::stoi(key.substr(coveredKey_.size()))].covered;
                for (std::size_t i = 0; i + 2 * dimension <= bounds.size(); i += 2 * dimension)
                {
                    core::Box<int, dimension> box;
                    for (std::size_t d = 0; d < dimension; ++d)
                    {
                        box.lower[d] = bounds[i + d];
                        box.upper[d] = bounds[i + dimension + d];
                    }
                    covered.push_back(box);
                }
            }
    }

    static inline std::string const restartName_ = "SolverFK";
    static inline std::string const coveredKey_  = "covered_level_";
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

    // E^{n-1/2} and B^n level ghosts, from the coarser level bracketing t_n - dt/2 and t_n. The B
    // ghosts were filled at the end of the previous substep; at the first substep of a coarse step
    // that fill predates the coarser level's reflux, which changed its B^n: refilled here, the
    // Ampere update reads the refluxed B^n (Bold_ of the messenger, at t_n: no interpolation), and
    // later substeps the same values as their previous fill
    if (levelNumber > 0)
    {
        messenger.fillElectricGhosts(E, level, currentTime - 0.5 * dt);
        messenger.fillMagneticGhosts(B, level, currentTime);
    }

    // the nodes the cleanings of this step skip (none without boundary_cleaning_mask)
    auto const& skipped = skippedNodes_(hierarchy, levelNumber, model);

    // 1. and 2.
    auto const& border = borderNodes_(hierarchy, levelNumber);
    advanceFKElectric(level, model, messenger, E, dt, currentTime, skipped, border);
    TimeSetter{level, model, halfTime}(E);
    messenger.fillElectricGhosts(E, level, halfTime);
    if (!border.empty()) // with the charge at the end of the current's segments, E's time
    {
        auto& rho = model.state.rho;
        for (auto& patch : rm.enumerate(level, rho, rhoStart_))
            rhoStart_.copyData(rho);
        depositFKCharge(level, model, messenger, dt, /*shift=*/1., halfTime);
        gaussOnLevelBorderOf(level, model, E, border);
        for (auto& patch : rm.enumerate(level, rho, rhoStart_))
            rho.copyData(rhoStart_);
        TimeSetter{level, model, currentTime - 0.5 * dt}(rho);
    }

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
        advanceFKElectric(level, model, messenger, Eahead_, dt, newTime, skipped, border);
        TimeSetter{level, model, newTime + 0.5 * dt}(Eahead_);
    }

    // the cells the finer level synchronizes at the end of this step, for the next step's mask
    keepCoveredRegion_(hierarchy, levelNumber, model);
}



// built once per level object and kept region; the level identity as in ParticleLevelBoxings
template<typename FKModel, typename AMR_Types>
auto SolverFK<FKModel, AMR_Types>::skippedNodes_(hierarchy_t const& hierarchy,
                                                 int const levelNumber, FKModel const& model)
    -> PatchBoxes<dimension> const&
{
    static PatchBoxes<dimension> const none;
    if (!model.boundary_cleaning_mask)
        return none;

    auto const level = hierarchy.getPatchLevel(levelNumber);
    auto& mask       = cleaningMasks_[levelNumber];
    // a new level object (a rebalanced root, a regridded level): the cells the finer level
    // synchronized at the end of the last step are the same, the patches are new
    if (mask.level.lock() != level)
    {
        mask.level = level;
        mask.stale = true;
    }
    if (mask.stale)
    {
        mask.skipped = cleaningSkippedNodes<GridLayout>(hierarchy, levelNumber, mask.covered);
        mask.stale   = false;
    }
    return mask.skipped;
}



// at the end of advanceLevel: the finer level advanced next is the one synchronized
template<typename FKModel, typename AMR_Types>
void SolverFK<FKModel, AMR_Types>::keepCoveredRegion_(hierarchy_t const& hierarchy,
                                                      int const levelNumber, FKModel const& model)
{
    if (!model.boundary_cleaning_mask)
        return;

    auto& mask       = cleaningMasks_[levelNumber];
    auto const finer = hierarchy.finerLevelExists(levelNumber)
                           ? hierarchy.getPatchLevel(levelNumber + 1)
                           : std::shared_ptr<level_t>{};
    if (mask.finer.lock() == finer and (finer or mask.covered.empty()))
        return;
    mask.finer   = finer;
    mask.covered = finerLevelCoverage<dimension>(hierarchy, levelNumber);
    mask.stale   = true;
}



template<typename FKModel, typename AMR_Types>
void SolverFK<FKModel, AMR_Types>::partialReflux(IPhysicalModel_t& level_model,
                                                 level_t& coarseLevel, level_t& fineLevel,
                                                 IMessenger& fromCoarserMessenger,
                                                 std::size_t const done, std::size_t const ratio,
                                                 double const midTime)
{
    PHARE_LOG_SCOPE(1, "SolverFK::partialReflux");

    auto& model     = dynamic_cast<FKModel&>(level_model);
    auto& messenger = dynamic_cast<Messenger&>(fromCoarserMessenger);
    auto& rm        = *model.resourcesManager;
    auto& E         = model.state.electromag.E;

    if (done == 0 or done >= ratio)
        throw std::runtime_error("SolverFK::partialReflux: not after a finer substep but the last");

    // finer level: the mean of its E over the substeps done (the flux sum weighs each 1 / ratio),
    // at t_n + done dt_f / 2, extrapolated to the coarse E time t_n + ratio dt_f / 2 with the
    // last two (Epartial_ holds the E of the substep start). Ghosts included, as the flux sum
    auto const mean  = static_cast<double>(ratio) / static_cast<double>(done);
    auto const trend = 0.5 * static_cast<double>(ratio - done);
    for (auto& patch : rm.enumerate(fineLevel, E, fluxSumE_, Epartial_))
        for (std::size_t c = 0; c < Epartial_.size(); ++c)
        {
            auto* partial    = Epartial_[c].data();
            auto const* sum  = fluxSumE_[c].data();
            auto const* last = E[c].data();
            for (std::size_t i = 0; i < Epartial_[c].size(); ++i)
                partial[i] = mean * sum[i] + trend * (last[i] - partial[i]);
        }

    // coarser level: its E^{n+1/2}, the covered edges from the finer level
    for (auto& patch : rm.enumerate(coarseLevel, E, Ehat_))
        Ehat_.copyData(E);
    messenger.partialReflux(coarseLevel.getLevelNumber(), fineLevel.getLevelNumber(), midTime);

    // the coarse B at midTime, from B^n (one coarsening, one Faraday, no ghost fill: the finer
    // level's B ghosts read it through the coarse interpolation of their fill)
    auto const dt = midTime - oldTime_.at(coarseLevel.getLevelNumber());
    Faraday_t{coarseLevel, model}(Bold_, Ehat_, Bmid_, dt);
    TimeSetter{coarseLevel, model, midTime}(Bmid_);
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

    // the particles that left the level go into every same-level band holding their cell, as
    // their own patch kept them (ion_updater.hpp), for the next step's band deposit. Out of the
    // level, they are not in a nonLevelGhostBox: the erase below leaves them
    if (level.getLevelNumber() > 0)
        messenger.fillBandLeavers(ions, level, newTime);

    // a level ghost copy pushed into a same-level neighbor's domain is that neighbor's domain
    // particle now: kept here, it would deposit twice (P4 plan D5: the neighbor's deposit reaches
    // this patch through the border sum) and could later be exported a second time
    for (auto& patch : rm.enumerate(level, ions))
        ions.eraseLevelGhostParticlesIn(
            levelBoxing.at(amr::to_string(patch->getGlobalId())).nonLevelGhostBox);

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
