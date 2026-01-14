#ifndef PHARE_MULTIPHYSICS_INTEGRATOR_HPP
#define PHARE_MULTIPHYSICS_INTEGRATOR_HPP

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>
#include <cstdint>
#include <unordered_map>


#include "core/def/phare_mpi.hpp" // IWYU pragma: keep

#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/hier/PatchDataRestartManager.h"


#include "amr/messengers/messenger.hpp"
#include "amr/tagging/tagger.hpp"
#include "amr/physical_models/hybrid_model.hpp"
#include "amr/physical_models/mhd_model.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"
#include "amr/messenger_registration.hpp"
#include "amr/level_initializer/level_initializer.hpp"
#include "amr/solvers/solver_mhd.hpp"
#include "amr/solvers/solver_ppc.hpp"

#include "core/logger.hpp"
#include "core/utilities/algorithm.hpp"

#include "load_balancing/load_balancer_manager.hpp"
#include "load_balancing/load_balancer_estimator.hpp"
#include "phare_core.hpp"


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
                            std::unique_ptr<PHARE::amr::Tagger> tagger)
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


            PHARE_LOG_LINE_SS("init level " << levelNumber << " with regriding = " << isRegridding
                                            << " and initial time = " << initialTime);

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

            if (isRegriddingL0)
            {
                messenger.registerLevel(hierarchy, levelNumber);
            }
            else if (isRegridding)
            {
                // regriding the current level has broken schedules for which
                // this level is the source or destination
                // we therefore need to rebuild them
                auto const finestLvlNbr = hierarchy->getFinestLevelNumber();
                auto nextFiner = (levelNumber == finestLvlNbr) ? levelNumber : levelNumber + 1;

                for (auto ilvl = levelNumber; ilvl <= nextFiner; ++ilvl)
                {
                    messenger.registerLevel(hierarchy, ilvl);
                }
                solver.onRegrid();
            }
            else
            {
                // we're not regriding, just making a new level
                messenger.registerLevel(hierarchy, levelNumber);
            }

            levelInitializer.initialize(hierarchy, levelNumber, oldLevel, model, messenger,
                                        initDataTime, isRegridding);

            if (isRegriddingL0)
            {
                for (auto ilvl = 1; ilvl <= hierarchy->getFinestLevelNumber(); ++ilvl)
                    messenger.registerLevel(hierarchy, ilvl);

                solver.onRegrid();
            }
            else
                load_balancer_manager_->estimate(*level, model);

            if (static_cast<std::size_t>(levelNumber) == model_views_.size())
                model_views_.push_back(solver.make_view(*level, model));
            else
                model_views_[levelNumber] = solver.make_view(*level, model);
        }



        void
        resetHierarchyConfiguration(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                    int const coarsestLevel, int const finestLevel) override
        {
            // handle samrai restarts / schedule creation
            // allocation of patch datas which may not want to be saved to restart files will
            // likely need to go here somehow https://github.com/PHAREHUB/PHARE/issues/664
            if (!restartInitialized_
                and SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
            {
                auto& messenger = getMessengerWithCoarser_(coarsestLevel);
                for (auto ilvl = coarsestLevel; ilvl <= finestLevel; ++ilvl)
                {
                    messenger.registerLevel(hierarchy, ilvl);

                    model_views_.push_back(getSolver_(ilvl).make_view(
                        AMR_Types::getLevel(*hierarchy, ilvl), getModel_(ilvl)));

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



        void initializeLevelIntegrator(
            std::shared_ptr<SAMRAI::mesh::GriddingAlgorithmStrategy> const& /*griddingAlg*/)
            override
        {
        }

        double getLevelDt(std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*level*/,
                          double const dtTime, bool const /*initialTime*/) override
        {
            return dtTime;
        }


        double getMaxFinerLevelDt(int const /*finerLevelNumber*/, double const coarseDt,
                                  SAMRAI::hier::IntVector const& ratio) override
        {
            // whistler waves require the dt ~ dx^2
            // so dividing the mesh size by ratio means dt
            // needs to be divided by ratio^2.
            // we multiply that by a constant < 1 for safety.
            return coarseDt / (ratio.max() * ratio.max()) /** 0.4*/;
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
         */

        double advanceLevel(std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                            std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            double const currentTime, double const newTime, bool const firstStep,
                            bool const lastStep, bool const regridAdvance = false) override
        {
            PHARE_LOG_SCOPE(3, "Multiphys::advanceLevel");

            if (regridAdvance)
                throw std::runtime_error("Error - regridAdvance must be False and is True");

            auto iLevel = level->getLevelNumber();

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
                fromCoarser.firstStep(model, *level, hierarchy, currentTime,
                                      subcycleStartTimes_[iLevel - 1],
                                      subcycleEndTimes_[iLevel - 1]);

                solver.resetFluxSum(model, *level);
            }

            solver.prepareStep(model, *level, currentTime);
            fromCoarser.prepareStep(model, *level, currentTime);

            solver.advanceLevel(*hierarchy, iLevel, getModelView_(iLevel), fromCoarser, currentTime,
                                newTime);

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
                auto ratio = (level->getRatioToCoarserLevel()).max();
                auto coef  = 1. / (ratio * ratio);
                solver.accumulateFluxSum(model, *level, coef);
            }

            load_balancer_manager_->estimate(*level, model);

            return newTime;
        }




        void
        standardLevelSynchronization(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                     int const coarsestLevel, int const finestLevel,
                                     double const syncTime,
                                     std::vector<double> const& /*oldTimes*/) override
        {
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
                    auto ratio = (coarseLevel.getRatioToCoarserLevel()).max();
                    auto coef  = 1. / (ratio * ratio);
                    coarseSolver.accumulateFluxSum(coarseModel, coarseLevel, coef);
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

        std::vector<std::shared_ptr<ISolverModelView>> model_views_;

        std::vector<std::shared_ptr<PHARE::amr::Tagger>> taggers_;
        std::map<std::string, std::unique_ptr<IMessengerT>> messengers_;
        std::map<std::string, std::unique_ptr<LevelInitializerT>> levelInitializers_;
        SimFunctors const& simFuncs_;
        PHARE::initializer::PHAREDict const& dict_;
        std::unique_ptr<amr::LoadBalancerManager<dimension>> load_balancer_manager_;


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


        void addTagger_(std::unique_ptr<PHARE::amr::Tagger> tagger, int coarsestLevel,
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


        auto& getModelView_(int iLevel) { return *model_views_[iLevel]; }


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


        amr::Tagger const& getTagger(int iLevel) const
        {
            auto& descriptor = levelDescriptors_[iLevel];
            if (taggers_[descriptor.taggerIndex] == nullptr)
            {
                throw std::runtime_error("Error - no tagger assigned to level "
                                         + std::to_string(iLevel));
            }
            return *taggers_[descriptor.taggerIndex];
        }



        amr::Tagger& getTagger_(int iLevel)
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
