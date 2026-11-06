#ifndef PHARE_MULTIPHYSICS_INTEGRATOR_HPP
#define PHARE_MULTIPHYSICS_INTEGRATOR_HPP

#include "phare_core.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep

#include "amr/amr_constants.hpp"
#include "amr/level_initializer/level_initializer.hpp"
#include "amr/messenger_registration.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/tagging/tagger.hpp"

#include "core/logger.hpp"
#include "core/utilities/algorithm.hpp"
#include "core/utilities/timestamps.hpp"

#include "mpi/mpi_utils.hpp"

#include "load_balancing/load_balancer_manager.hpp"

#include "SAMRAI/algs/TimeRefinementLevelStrategy.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/tbox/RestartManager.h"

#include <algorithm>
#include <map>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>


namespace PHARE
{
namespace solver
{
    struct LevelDescriptor
    {
        static int constexpr NOT_SET = -1;
        int modelIndex               = NOT_SET;
        int solverIndex              = NOT_SET;
        int resourcesManagerIndex    = NOT_SET;
        int taggerIndex              = NOT_SET;
        std::string messengerName;
    };


    inline bool isRootLevel(int levelNumber)
    {
        return levelNumber == 0;
    }


    /**
     * @brief The MultiPhysicsIntegrator is given to the SAMRAI GriddingAlgorithm and is used by
     * SAMRAI to manipulate data of patch levels for application-dependent tasks such as
     * initializing a level or advancing a level. It inherits from
     * SAMRAI::algs::TimeRefinementLevelStrategy and SAMRAI::mesh::StandardTagAndInitStrategy and
     * implement the needed application-dependent operations.
     *
     * The MultiPhysicsIntegrator uses a pool of IMessenger, ISolver and IPhysicalModel to
     * initialize, advance and synchronize a given level. It manipulate these abstractions and
     * therefore does not know any details of the specific models, solvers and messengers used to
     * manipulate a level. Each time SAMRAI calls the MultiPhysicsIntegrator to perform an action at
     * a specific level, the MultiPhysicsIntegrator usually first has to figure out which of the
     * IMessenger, IModel and ISolver to use.
     *
     * The interface of the MultiPhysicsIntegrator is composed of:
     *
     * - the implementations of the StandardTagAndInitStrategy and TimeRefinementLevelStrategy
     * - methods to register the models, messengers and solvers for a given PatchHierarchy.
     *
     * these last methods are to be called in the following order:
     *
     * - registerModel() : used to register a model for a range of levels
     * - registerAndInitSolver(): used to register and initialize a solver for a range of levels
     * where the compatible model has been registered
     * - registerAndSetupMessengers() : used to register the messengers associated with the
     * registered IPhysicalModel and ISolver objects
     *
     */
    template<typename MessengerFactory, typename LevelInitializerFactory, typename AMR_Types>
    class MultiPhysicsIntegrator : public SAMRAI::mesh::StandardTagAndInitStrategy,
                                   public SAMRAI::algs::TimeRefinementLevelStrategy
    {
        using SimFunctorParams = typename core::PHARE_Sim_Types::SimFunctorParams;
        using SimFunctors      = typename core::PHARE_Sim_Types::SimulationFunctors;

    public:
        static constexpr auto dimension = MessengerFactory::dimension;

        // model comes with its variables already registered to the manager system
        MultiPhysicsIntegrator(PHARE::initializer::PHAREDict const& dict,
                               SimFunctors const& simFuncs)
            : nbrOfLevels_{dict["AMR"]["max_nbr_levels"].template to<int>()}
            , levelDescriptors_(dict["AMR"]["max_nbr_levels"].template to<int>())
            , simFuncs_{simFuncs}
            , dict_{dict}
            , load_balancer_manager_{nullptr}
            , isAdaptive_{cppdict::get_value(dict, "time_step/mode", core::TimeStepType::constant)
                          == core::TimeStepType::adaptive}
            , constantDt_{cppdict::get_value(dict, "time_step/value", 0.)}
            , levelDts_(nbrOfLevels_, std::numeric_limits<double>::quiet_NaN())
            , cfl_{cppdict::get_value(dict, "time_step/cfl_wave", 0.),
                   cppdict::get_value(dict, "time_step/cfl_diffusive", 0.)}
            , fineLevelSubstepsMax_{
                  cppdict::get_value(dict, "time_step/fine_level_substeps_max", 1000)}

        {
            // auto mhdSolver = std::make_unique<SolverMHD<ResourcesManager>>(resourcesManager_);
            // solvers.push_back(std::move(mhdSolver));

            // auto hybridSolver = std::make_unique<
            //    SolverPPC<ResourcesManager, decltype(hybridModel.electromag), HybridState>>(
            //    hybridModel, resourcesManager_);

            // solvers.push_back(std::move(hybridSolver));



            //@TODO - chaque modele utilisé doit register ses variables aupres du ResourcesManager
            //@TODO - chaque solveur utilisé doit register ses variables aupres du ResourcesManager
        }



        /* -------------------------------------------------------------------------
                           MultiPhysicsIntegrator proper interface
           ------------------------------------------------------------------------- */


        auto nbrOfLevels() const { return nbrOfLevels_; }


        void registerTagger(int coarsestLevel, int finestLevel,
                            std::unique_ptr<PHARE::amr::Tagger<AMR_Types>> tagger)
        {
            if (!validLevelRange_(coarsestLevel, finestLevel))
            {
                throw std::runtime_error("invalid range level");
            }
            if (existTaggerOnRange_(coarsestLevel, finestLevel))
            {
                throw std::runtime_error(
                    "error - level range contains levels with a registered tagger");
            }
            addTagger_(std::move(tagger), coarsestLevel, finestLevel);
        }



        /**
         * @brief registerModel registers the model to the multiphysics integrator for a given level
         * range. The level index for the coarsest and finest must be greater or equal to zero, less
         * than the nbrOfLevels(). Once this method is called, the MultiPhysicsIntegrator will use
         * this model for all levels in the given range.
         *
         * @param coarsestLevel is the index of the coarsest level using the model
         * @param finestLevel is the index of the finest level using the model. finestLevel >
         * coarsestModel
         * @param model the model to be registered to the MultiphysicsIntegrator. The model must not
         * have been registered already.
         */
        void registerModel(int coarsestLevel, int finestLevel,
                           std::shared_ptr<IPhysicalModel<AMR_Types>> model)
        {
            if (!validLevelRange_(coarsestLevel, finestLevel))
            {
                throw std::runtime_error("invalid range level");
            }

            if (existModelOnRange_(coarsestLevel, finestLevel))
            {
                throw std::runtime_error(
                    "error - level range contains levels with registered model");
            }


            addModel_(model, coarsestLevel, finestLevel);
        }




        /**
         * @brief registerAndInitSolver registers and initialize a given solver to the
         * MultiphysicsIntegrator
         *
         * Once this method is called :
         *
         * - the MultiPhysicsIntegrator will use this solver to advance the
         * hierarchy for levels between coarsestLevel and finestLevel (included). The level index
         * for the coarsest and finest must be greater or equal to zero, less than the
         * nbrOfLevels().
         *
         * - the solver will have registered its variables to the ResourcesManager of the model in
         * the given range.
         *
         * @param coarsestLevel is the index of the coarsest level using the model
         * @param finestLevel is the index of the finest level using the model. finestLevel >
         * coarsestModel
         * @param solver is the ISolver to register to the MultiPhysicsIntegrator. This solver must
         * be compatible with all models in the given range and not yet registered.
         */
        void registerAndInitSolver(int coarsestLevel, int finestLevel,
                                   std::unique_ptr<ISolver<AMR_Types>> solver)
        {
            if (!validLevelRange_(coarsestLevel, finestLevel))
            {
                throw std::runtime_error("invalid level range");
            }

            if (!canBeRegistered_(coarsestLevel, finestLevel, *solver))
            {
                throw std::runtime_error(
                    solver->name() + " is not compatible with model on specified level range");
            }


            auto& model = getModel_(coarsestLevel);
            solver->registerResources(model);


            addSolver_(std::move(solver), coarsestLevel, finestLevel);
        }




        /**
         * @brief registerAndSetupMessengers registers and setup the messengers for the hierarchy
         *
         * This method create all necessary messengers for levels in the hierarchy to exchange data
         * knowing the sequence of models. This method must thus be called *after* the
         * registerModel().
         *
         *
         * @param messengerFactory is used to create the appropriate messenger
         */
        void registerAndSetupMessengers(MessengerFactory& messengerFactory)
        {
            registerMessengers_(messengerFactory);

            // now setup all messengers we've just created

            registerQuantitiesAllLevels_();
        }




        std::string solverName(int const iLevel) const { return getSolver_(iLevel).name(); }




        std::string modelName(int const iLevel) const { return getModel_(iLevel).name(); }




        std::string messengerName(int iLevel)
        {
            auto& messenger = getMessengerWithCoarser_(iLevel);
            return messenger.name();
        }



        // -----------------------------------------------------------------------------------------------
        //
        //                          SAMRAI StandardTagAndInitStrategy interface
        //
        // -----------------------------------------------------------------------------------------------




        /**
         * @brief see SAMRAI documentation. This function initializes the data on the given level.
         *
         * The method first checks wether allocation of patch data must be performed.
         * If it does, all objects using resources on patches must see their allocate() function
         * called.
         *
         * This is:
         * - the model (by definition the model has data defined on patches)
         * - the solver (Some solvers has internal data that needs to exist on patches)
         * - the messenger (the messenger has data defined on patches for internal reasons)
         *
         *
         * then the level needs to be registered to the messenger.
         *
         * Then data initialization per se begins and one can be on one of the following cases:
         *
         * - regridding
         * - initialization of the root level
         * - initialization of a new level from scratch (not a regridding)
         *
         */
        void initializeLevelData(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                 int const levelNumber, double const initDataTime,
                                 bool const canBeRefined, bool const initialTime,
                                 std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel
                                 = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                                 bool const allocateData = true) override
        {
            auto& model            = getModel_(levelNumber);
            auto& solver           = getSolver_(levelNumber);
            auto& messenger        = getMessengerWithCoarser_(levelNumber);
            auto& levelInitializer = getLevelInitializer(model.name());

            bool const isRegridding = oldLevel != nullptr;
            auto level              = hierarchy->getPatchLevel(levelNumber);


            PHARE_LOG_LINE_SS("init level " << levelNumber << " with regriding = " << isRegridding);

            PHARE_LOG_START(3, "initializeLevelData::allocate block");
            if (allocateData)
            {
                for (auto patch : *level)
                {
                    model.allocate(*patch, initDataTime);
                    solver.allocate(model, *patch, initDataTime);
                    messenger.allocate(*patch, initDataTime);
                    load_balancer_manager_->allocate(*patch, initDataTime);
                }
            }

            PHARE_LOG_STOP(3, "initializeLevelData::allocate block");

            bool const isRegriddingL0 = isRegridding and levelNumber == 0;

            {
                PHARE_LOG_SCOPE(1, "Multiphys::initializeLevelData::registerLevel");
                messenger.registerLevel(hierarchy, levelNumber);
            }

            if (isRegridding and not isRegriddingL0)
                solver.onRegrid();

            levelInitializer.initialize(hierarchy, levelNumber, oldLevel, model, messenger,
                                        initDataTime, isRegridding);

            if (isRegriddingL0)
                solver.onRegrid();
            else
                load_balancer_manager_->estimate(*level, model);
        }



        void
        resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                    int const coarsestLevel, int const finestLevel) override
        {
            // handle samrai restarts / schedule creation
            // allocation of patch datas which may not want to be saved to restart files will
            // likely need to go here somehow https://github.com/PHAREHUB/PHARE/issues/664
            if (!restartInitialized_ and SamraiLifeCycle::getRestartManager()->isFromRestart())
            {
                auto& messenger = getMessengerWithCoarser_(coarsestLevel);
                for (auto ilvl = coarsestLevel; ilvl <= finestLevel; ++ilvl)
                {
                    messenger.registerLevel(hierarchy, ilvl);

                    auto level = hierarchy->getPatchLevel(ilvl);
                    for (auto& patch : *level)
                    {
                        auto time = dict_["restarts"]["restart_time"].template to<double>();
                        load_balancer_manager_->allocate(*patch, time);
                    }
                    // if load balance on restart advance
                    load_balancer_manager_->estimate(*level, getModel_(ilvl));
                }
                restartInitialized_ = true;
            }

            PHARE_LOG_SCOPE(1, "Multiphys::resetHierarchyConfiguration::onHierarchyChanged");
            // coarsestLevel and finestLevel cannot bound the levels to notify: after a level 0
            // rebalance SAMRAI passes (0, 0) while finer levels read level 0 directly or through
            // recursive coarse interpolation. Each messenger works out what is stale.
            for (auto ilvl = 0; ilvl <= hierarchy->getFinestLevelNumber(); ++ilvl)
                getMessengerWithCoarser_(ilvl).onHierarchyChanged(hierarchy, ilvl);
        }



        void applyGradientDetector(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   int const levelNumber, double const /*error_data_time*/,
                                   int const tag_index, bool const /*initialTime*/,
                                   bool const /*usesRichardsonExtrapolationToo*/) override
        {
            PHARE_LOG_LINE_STR("apply gradient detector on level " + std::to_string(levelNumber));

            auto level = hierarchy->getPatchLevel(levelNumber);
            for (auto& patch : *level)
            {
                auto& model  = getModel_(levelNumber);
                auto& tagger = getTagger_(levelNumber);
                tagger.tag(model, *patch, tag_index);
            }
        }


        // -----------------------------------------------------------------------------------------------
        //
        //                          SAMRAI TimeRefinementLevelStrategy interface
        //
        // -----------------------------------------------------------------------------------------------


        /**
         * @brief Called by SAMRAI first in TimeRefinementIntegrator::initializeHierarchy, at
         * startup and on restart.
         *
         * Currently it only precomputes each level's dt in constant time-stepping mode.
         */
        void initializeLevelIntegrator(
            std::shared_ptr<SAMRAI::mesh::GriddingAlgorithmStrategy> const& /*griddingAlg*/)
            override
        {
            if (!isAdaptive_)
                levelDts_ = constantLevelDts_();
        }

        /**
         * @brief Return the dt of a level created while SAMRAI builds the hierarchy at startup.
         *
         * Only called by SAMRAI at startup, in
         * TimeRefinementIntegrator::initializeRefinedTimesteppingLevelData, for the levels that
         * exist at that time (never on restart). For a level > 0, SAMRAI keeps the min between
         * this value and getMaxFinerLevelDt.
         *
         * - adaptive mode: computes the level's stable dt on the initial state with
         *   levelStableDt_ and stores it in levelDts_.
         * - constant mode: reads levelDts_, already filled by initializeLevelIntegrator.
         *
         * @param level the level SAMRAI has just created and initialized
         * @param time the simulation start time
         * @return the dt of the level
         */
        double getLevelDt(std::shared_ptr<SAMRAI::hier::PatchLevel> const& level, double const time,
                          bool const /*initialTime*/) override
        {
            auto const levelNumber = level->getLevelNumber();
            if (isAdaptive_)
                levelDts_.at(levelNumber) = levelStableDt_(levelNumber, *level, cfl_, time);
            return levelDts_.at(levelNumber);
        }

        /**
         * @brief Return the dt bound of a finer level from the dt of its coarser level:
         * coarseDt / r^2.
         *
         * Only place that implements the time ratio logic between levels. It mostly applies to
         * constant time-stepping mode, but also has some effect in adaptive time-stepping (see
         * below), so be careful when modifying it.
         *
         * It is called by:
         * - constantLevelDts_, in initializeLevelIntegrator, to build the dt of each level in
         *   constant time-stepping mode.
         * - SAMRAI's TimeRefinementIntegrator::initializeRefinedTimesteppingLevelData, at
         *   simulation startup, for the levels > 0 that exist at that time: their first dt is
         *   the min between getMaxFinerLevelDt and getLevelDt.
         * - SAMRAI's TimeRefinementIntegrator::advanceRecursivelyForRefinedTimestepping, at the
         *   first advance of a level created by a regrid: getLevelDt is not called, and the
         *   first dt of the new level is bounded by getMaxFinerLevelDt (and by the coarser dt).
         *   So even in adaptive mode, the first substep of a freshly created level is set by
         *   this fixed time ratio rather than by its state; afterwards the dt returned by
         *   advanceLevel takes over. Basing it on the state would require calling levelStableDt_
         *   here, but as only the level number is passed, the MultiPhysicsIntegrator would need
         *   the hierarchy as a member.
         *
         * In the future, the current hardcoded rule dt_finer = dt_coarser / r^2 could be
         * replaced by a call to the finer level's solver, responsible for deciding, from its own
         * model and from the coarser level's model and solver, how the finer level's dt is
         * bounded by the coarser one. Both are needed because the Hall term is a solver option,
         * not a model property: e.g. ideal MHD over ideal MHD could use r, Hall MHD over Hall MHD
         * needs r^2, and a hybrid level should only sit below Hall MHD.
         *
         * @param coarseDt the dt of the coarser level
         * @param ratio the refinement ratio between the two levels
         * @return the maximum dt allowed on the finer level
         */
        double getMaxFinerLevelDt(int const /*finerLevelNumber*/, double const coarseDt,
                                  SAMRAI::hier::IntVector const& ratio) override
        {
            // whistler waves require the dt ~ dx^2
            // so dividing the mesh size by ratio means dt
            // needs to be divided by ratio^2.
            return coarseDt / (ratio.max() * ratio.max());
        }


        double computeStableDt(AMR_Types::hierarchy_t& hierarchy, CFLNumbers const& stability,
                               double const time)
        {
            return levelDts_.at(0)
                   = levelStableDt_(0, *hierarchy.getPatchLevel(0), stability, time);
        }


        void dump_(int iLevel)
        {
            // we skip first/last as that's done via regular diag dump mechanism
            bool fineDumpsActive = simFuncs_.at("pre_advance").count("fine_dump") > 0;
            bool notCoarsestTime = subcycleEndTimes_[iLevel] != subcycleEndTimes_[0];
            bool shouldDump      = fineDumpsActive and notCoarsestTime and iLevel > 0;

            if (shouldDump)
            {
                auto const& dump_functor = simFuncs_.at("pre_advance").at("fine_dump");
                SimFunctorParams fineDumpParams;
                fineDumpParams["level_nbr"] = iLevel;
                fineDumpParams["timestamp"] = subcycleEndTimes_[iLevel];

                dump_functor(fineDumpParams);
            }
        }


        /**
         * @brief advanceLevel is derived from the abstract method of TimeRefinementLevelStrategy
         *
         * In this method, the MultiPhysicsIntegrator needs to get the model, solver and messenger
         * necessary to advance the given level. It then forwards the call to the solver's
         * advanceLevel method, passing it the model and messenger, among others.
         *
         * If it is the first step of the subcycling the Messenger may have something to do before
         * calling the solver's advanceLevel(). Typically it may need to grab coarser level's data
         * at the coarser future time so that all subsequent subcycles may use time interpolation
         * without involving communications. This is done by calling messenger.firstStep()
         *
         *
         * At the last step of the subcycle, the Messenger may also need to perform some actions,
         * like working on its internal data for instance. messenger.lastStep()
         *
         * The returned value is the dt bound SAMRAI uses for the next substeps of the level:
         * - constant mode: levelDts_[L], the constant dt of the level.
         * - adaptive mode, level 0: the step just taken. The next level-0 step is set by the
         *   simulator (computeStableDt), so this value only updates SAMRAI's bookkeeping.
         * - adaptive mode, level > 0: the min between the coarser level's substep and the stable
         *   dt of the level at newTime.
         *
         * For a level > 0, SAMRAI's TimeRefinementIntegrator::findNextDtAndStepsRemaining then
         * chooses the next substep:
         * - bound = min(returned value, previous bound * grow_dt). grow_dt is a SAMRAI param >=1
         *   set in the Integrator constructor.
         * - n = the number of substeps that fit in the time remaining until the end of the
         *   coarser step: int(time remaining / bound), plus one if more than a round-off
         *   remainder is left. For a level that can be refined, n is rounded up to a multiple of
         *   the regrid interval.
         * - next substep = time remaining / n, so the substeps left in a coarser step are equal
         *   and at most the bound.
         *
         *
         * @return the dt bound for the next substeps of the level
         */

        double advanceLevel(std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                            std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            double const currentTime, double const newTime, bool const firstStep,
                            bool const lastStep, bool const regridAdvance = false) override
        {
            PHARE_LOG_SCOPE(1, "Multiphys::advanceLevel");

            if (regridAdvance)
                throw std::runtime_error("Error - regridAdvance must be False and is True");

            auto iLevel = level->getLevelNumber();

            if (iLevel == 0 && !(newTime > currentTime))
            {
                std::ostringstream ss;
                ss << std::setprecision(17) << "zero-length step on level 0: t = " << currentTime
                   << ", newTime = " << newTime;
                throw std::runtime_error(ss.str());
            }

            PHARE_LOG_LINE_SS("advanceLevel " << iLevel << " with dt = " << newTime - currentTime
                                              << " from t = " << currentTime
                                              << " to t = " << newTime);

            auto& solver      = getSolver_(iLevel);
            auto& model       = getModel_(iLevel);
            auto& fromCoarser = getMessengerWithCoarser_(iLevel);


            subcycleStartTimes_[iLevel] = currentTime;
            subcycleEndTimes_[iLevel]   = newTime;

            if (firstStep)
            {
                PHARE_LOG_SCOPE(3, "Multiphys::advanceLevel.firstStep");
                fromCoarser.firstStep(model, *level, currentTime, subcycleStartTimes_[iLevel - 1],
                                      subcycleEndTimes_[iLevel - 1]);

                solver.resetFluxSum(model, *level);
            }

            solver.prepareStep(model, *level, currentTime);
            fromCoarser.prepareStep(model, *level, currentTime);

            solver.advanceLevel(*hierarchy, iLevel, model, fromCoarser, currentTime, newTime);

            if (lastStep)
            {
                PHARE_LOG_SCOPE(3, "Multiphys::advanceLevel.lastStep");
                fromCoarser.lastStep(model, *level);
            }

            if (iLevel == hierarchy->getFinestLevelNumber())
            {
                dump_(iLevel);
            }

            // If we are on the finest level, we cannot do the accumulation in the
            // standardLevelSynchronization, so we need to do it per substep here
            if (iLevel != 0 && !hierarchy->finerLevelExists(iLevel))
            {
                solver.accumulateFluxSum(model, *level, substepWeight_(iLevel));
            }

            load_balancer_manager_->estimate(*level, model);

            if (!isAdaptive_)
                return levelDts_.at(iLevel);

            if (iLevel == 0)
                return newTime - currentTime;

            auto const coarseDt
                = subcycleEndTimes_.at(iLevel - 1) - subcycleStartTimes_.at(iLevel - 1);
            levelDts_.at(iLevel) = levelStableDt_(iLevel, *level, cfl_, newTime);
            return std::min(coarseDt, levelDts_[iLevel]);
        }




        void
        standardLevelSynchronization(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                     int const coarsestLevel, int const finestLevel,
                                     double const syncTime,
                                     std::vector<double> const& /*oldTimes*/) override
        {
            PHARE_LOG_SCOPE(1, "Multiphys::standardLevelSynchronization");

            // TODO use messengers to sync with coarser
            for (auto ilvl = finestLevel; ilvl > coarsestLevel; --ilvl)
            {
                auto& toCoarser = getMessengerWithCoarser_(ilvl);
                auto& fineLevel = *hierarchy->getPatchLevel(ilvl);
                toCoarser.synchronize(fineLevel);

                // refluxing
                auto& fineSolver   = getSolver_(ilvl);
                auto iCoarseLevel  = ilvl - 1;
                auto& coarseLevel  = *hierarchy->getPatchLevel(iCoarseLevel);
                auto& coarseSolver = getSolver_(iCoarseLevel);
                auto& coarseModel  = getModel_(iCoarseLevel);

                toCoarser.reflux(iCoarseLevel, ilvl, syncTime);
                coarseSolver.reflux(coarseModel, coarseLevel, toCoarser, syncTime);

                // Now the fluxSum includes the contributions of the finer levels thanks to
                // toCoarser.reflux(). We can now accumulate the fluxSum that will be used for the
                // next coarser reflux.
                // It is not needed to do this accumulation on the coarsest level.
                if (iCoarseLevel != 0)
                {
                    coarseSolver.accumulateFluxSum(coarseModel, coarseLevel,
                                                   substepWeight_(iCoarseLevel));
                }

                // recopy (patch) ghosts
                toCoarser.postSynchronize(coarseModel, coarseLevel, syncTime);

                // advancing all but the finest includes synchronization of the finer
                dump_(iCoarseLevel);
            }
        }




        void
        synchronizeNewLevels(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                             int const /*coarsestLevel*/, int const /*finestLevel*/,
                             double const /*syncTime*/, bool const /*initialTime*/) override
        {
        }


        void resetTimeDependentData(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/,
                                    double const /*newTime*/, bool const /*canBeRefined*/) override
        {
        }

        void resetDataToPreadvanceState(
            std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/) override
        {
        }

        bool usingRefinedTimestepping() const override { return true; }


        void setLoadBalancerManager(std::unique_ptr<amr::LoadBalancerManager<dimension>> lbm)
        {
            load_balancer_manager_ = std::move(lbm);
        }



    private:
        bool restartInitialized_ = false;
        int nbrOfLevels_;
        std::unordered_map<std::size_t, double> subcycleEndTimes_;
        std::unordered_map<std::size_t, double> subcycleStartTimes_;
        using IMessengerT       = amr::IMessenger<IPhysicalModel<AMR_Types>>;
        using LevelInitializerT = LevelInitializer<AMR_Types>;
        std::vector<LevelDescriptor> levelDescriptors_;
        std::vector<std::unique_ptr<ISolver<AMR_Types>>> solvers_;
        std::vector<std::shared_ptr<IPhysicalModel<AMR_Types>>> models_;

        std::vector<std::shared_ptr<PHARE::amr::Tagger<AMR_Types>>> taggers_;
        std::map<std::string, std::unique_ptr<IMessengerT>> messengers_;
        std::map<std::string, std::unique_ptr<LevelInitializerT>> levelInitializers_;
        SimFunctors const& simFuncs_;
        PHARE::initializer::PHAREDict const& dict_;
        std::unique_ptr<amr::LoadBalancerManager<dimension>> load_balancer_manager_;
        bool const isAdaptive_;
        double const constantDt_;
        /**
         * @brief The dt of each level.
         *
         * - constant mode: each level's dt, built by initializeLevelIntegrator.
         * - adaptive mode: each level's last stable dt (NaN until computed), read by
         *   levelStableDt_ to check that level L's stable dt is not more than
         *   fine_level_substeps_max times smaller than level L-1's.
         */
        std::vector<double> levelDts_;
        CFLNumbers const cfl_;
        int const fineLevelSubstepsMax_; //!< limit on the number of substeps of a fine level


        /**
         * @brief Compute the list of level dts in constant time-stepping mode, based on the
         * inter-level time ratio logic implemented in getMaxFinerLevelDt.
         */
        std::vector<double> constantLevelDts_()
        {
            SAMRAI::hier::IntVector const ratio{SAMRAI::tbox::Dimension(dimension),
                                                static_cast<int>(amr::refinementRatio)};
            std::vector<double> dts{constantDt_};
            for (int iLevel = 1; iLevel < nbrOfLevels_; ++iLevel)
                dts.push_back(
                    MultiPhysicsIntegrator::getMaxFinerLevelDt(iLevel, dts.back(), ratio));
            return dts;
        }

        /**
         * @brief Compute the level's stable dt, reduced across MPI ranks.
         *
         * For a level > 0, throw if it is more than fine_level_substeps_max times smaller than
         * the coarser level's stable dt, or if that one was never computed (NaN). The caller
         * stores the result in levelDts_.
         *
         * @param iLevel the level number
         * @param level the level
         * @param cfl the CFL numbers
         * @param time the time of the level's data, used in the error message
         * @return the stable dt of the level
         */
        double levelStableDt_(int const iLevel, SAMRAI::hier::PatchLevel& level,
                              CFLNumbers const& cfl, double const time)
        {
            auto const dt
                = mpi::min(getSolver_(iLevel).computeStableDt(getModel_(iLevel), level, cfl));
            if (iLevel > 0)
            {
                auto const coarserDt = levelDts_.at(iLevel - 1);
                if (!(coarserDt <= fineLevelSubstepsMax_ * dt))
                {
                    std::ostringstream ss;
                    ss << std::setprecision(17) << "stable dt on level " << iLevel
                       << " at t = " << time << " is " << dt
                       << ", more than fine_level_substeps_max = " << fineLevelSubstepsMax_
                       << " times smaller than on level " << iLevel - 1 << " (" << coarserDt << ")";
                    throw std::runtime_error(ss.str());
                }
            }
            return dt;
        }


        double substepWeight_(int const iLevel) const
        {
            auto const coarseDuration
                = subcycleEndTimes_.at(iLevel - 1) - subcycleStartTimes_.at(iLevel - 1);
            if (!(coarseDuration > 0.))
                return 0.;
            return (subcycleEndTimes_.at(iLevel) - subcycleStartTimes_.at(iLevel)) / coarseDuration;
        }


        bool validLevelRange_(int coarsestLevel, int finestLevel)
        {
            if (coarsestLevel < 0 || finestLevel >= nbrOfLevels_ || finestLevel < coarsestLevel)
            {
                return false;
            }
            return true;
        }



        bool existTaggerOnRange_(int coarsestLevel, int finestLevel)
        {
            for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
            {
                if (levelDescriptors_[iLevel].taggerIndex != LevelDescriptor::NOT_SET)
                {
                    return true;
                }
            }
            return false;
        }



        bool existModelOnRange_(int coarsestLevel, int finestLevel)
        {
            bool hasModel = true;

            for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
            {
                if (levelDescriptors_[iLevel].modelIndex != LevelDescriptor::NOT_SET)
                {
                    return hasModel;
                }
            }
            return !hasModel;
        }


        void addTagger_(std::unique_ptr<PHARE::amr::Tagger<AMR_Types>> tagger, int coarsestLevel,
                        int finestLevel)
        {
            if (core::notIn(tagger, taggers_))
            {
                taggers_.push_back(std::move(tagger));
                int taggerIndex = taggers_.size() - 1;
                for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
                {
                    levelDescriptors_[iLevel].taggerIndex = taggerIndex;
                }
            }
        }



        void addModel_(std::shared_ptr<IPhysicalModel<AMR_Types>> model, int coarsestLevel,
                       int finestLevel)
        {
            if (core::notIn(model, models_))
            {
                levelInitializers_[model->name()]
                    = LevelInitializerFactory::create(model->name(), dict_);
                models_.push_back(std::move(model));
                int modelIndex = models_.size() - 1;

                for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
                {
                    levelDescriptors_[iLevel].modelIndex = modelIndex;
                }
            }
            else
            {
                throw std::runtime_error("model " + model->name() + " already registered");
            }
        }




        void addSolver_(std::unique_ptr<ISolver<AMR_Types>> solver, int coarsestLevel,
                        int finestLevel)
        {
            if (core::notIn(solver, solvers_))
            {
                solvers_.push_back(std::move(solver)); // check that solver exist

                for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
                {
                    levelDescriptors_[iLevel].solverIndex = solvers_.size() - 1;
                }
            }
            else
            {
                throw std::runtime_error("solver " + solver->name() + " already registered");
            }
        }




        bool canBeRegistered_(int coarsestLevel, int finestLevel, ISolver<AMR_Types> const& solver)
        {
            bool itCan = true;

            for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
            {
                if (auto& model = getModel_(iLevel); !areCompatible(model, solver))
                {
                    return false;
                }
            }
            return itCan;
        }




        void registerMessengers_(MessengerFactory& messengerFactory)
        {
            for (auto iLevel = 0; iLevel < nbrOfLevels_; ++iLevel)
            {
                auto coarseLevelNumber = iLevel - 1;
                auto fineLevelNumber   = iLevel;

                if (iLevel == 0)
                {
                    coarseLevelNumber = fineLevelNumber;
                }

                auto& coarseModel = getModel_(coarseLevelNumber);
                auto& fineModel   = getModel_(fineLevelNumber);

                registerMessenger_(messengerFactory, coarseModel, fineModel, iLevel);
            }
        }




        void registerMessenger_(MessengerFactory const& factory,
                                IPhysicalModel<AMR_Types> const& coarseModel,
                                IPhysicalModel<AMR_Types> const& fineModel, int iLevel)
        {
            if (auto messengerName = factory.name(coarseModel, fineModel); messengerName)
            {
                auto foundMessenger = messengers_.find(*messengerName);
                if (foundMessenger == messengers_.end())
                {
                    messengers_[*messengerName]
                        = std::move(factory.create(*messengerName, coarseModel, fineModel, iLevel));
                }

                levelDescriptors_[iLevel].messengerName = *messengerName;
            }
            else
            {
                throw std::runtime_error("No viable messenger found");
            }
        }




        void registerQuantities_(int iLevel, IMessengerT& messenger)
        {
            auto coarseLevelNumber = iLevel - 1;
            auto fineLevelNumber   = iLevel;

            if (iLevel == 0)
            {
                coarseLevelNumber = fineLevelNumber;
            }

            auto& coarseModel = getModel_(coarseLevelNumber);
            auto& fineModel   = getModel_(fineLevelNumber);
            auto& solver      = getSolver_(iLevel);

            MessengerRegistration::registerQuantities(messenger, coarseModel, fineModel, solver);
        }




        void registerQuantitiesAllLevels_()
        {
            auto lastMessengerName = levelDescriptors_[0].messengerName;
            registerQuantities_(0, getMessengerWithCoarser_(0));

            for (auto iLevel = 1; iLevel < nbrOfLevels_; ++iLevel)
            {
                auto currentMessengerName = levelDescriptors_[iLevel].messengerName;
                if (lastMessengerName != currentMessengerName)
                {
                    lastMessengerName = currentMessengerName;
                    registerQuantities_(iLevel, getMessengerWithCoarser_(iLevel));
                }
            }
        }




        ISolver<AMR_Types>& getSolver_(int iLevel)
        {
            return const_cast<ISolver<AMR_Types>&>(
                const_cast<std::remove_pointer_t<decltype(this)> const*>(this)->getSolver_(iLevel));
        }




        ISolver<AMR_Types> const& getSolver_(int iLevel) const
        {
            auto& descriptor = levelDescriptors_[iLevel];
            return *solvers_[descriptor.solverIndex];
        }



        IPhysicalModel<AMR_Types>& getModel_(int iLevel)
        {
            return const_cast<IPhysicalModel<AMR_Types>&>(
                const_cast<std::remove_pointer_t<decltype(this)> const*>(this)->getModel_(iLevel));
        }




        IPhysicalModel<AMR_Types> const& getModel_(int iLevel) const
        {
            auto& descriptor = levelDescriptors_[iLevel];
            if (models_[descriptor.modelIndex] == nullptr)
                throw std::runtime_error("Error - no model assigned to level "
                                         + std::to_string(iLevel));
            return *models_[descriptor.modelIndex];
        }


        amr::Tagger<AMR_Types> const& getTagger(int iLevel) const
        {
            auto& descriptor = levelDescriptors_[iLevel];
            if (taggers_[descriptor.taggerIndex] == nullptr)
            {
                throw std::runtime_error("Error - no tagger assigned to level "
                                         + std::to_string(iLevel));
            }
            return *taggers_[descriptor.taggerIndex];
        }



        amr::Tagger<AMR_Types>& getTagger_(int iLevel)
        {
            auto& descriptor = levelDescriptors_[iLevel];
            if (taggers_[descriptor.taggerIndex] == nullptr)
            {
                throw std::runtime_error("Error - no tagger assigned to level "
                                         + std::to_string(iLevel));
            }
            return *taggers_[descriptor.taggerIndex];
        }




        LevelInitializerT& getLevelInitializer(std::string modelName)
        {
            assert(levelInitializers_.count(modelName));
            return *levelInitializers_[modelName];
        }



        IMessengerT& getMessengerWithCoarser_(int iLevel)
        {
            auto& descriptor = levelDescriptors_[iLevel];
            auto messenger   = messengers_[descriptor.messengerName].get();

            if (messenger)
                return *messenger;
            else
                throw std::runtime_error("no found messenger");
        }
    };
} // namespace solver
} // namespace PHARE
#endif
