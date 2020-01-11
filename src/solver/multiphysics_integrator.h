
#ifndef PHARE_MULTIPHYSICS_INTEGRATOR_H
#define PHARE_MULTIPHYSICS_INTEGRATOR_H

#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>


#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>


#include "amr/messengers/messenger.h"

#include "solver/physical_models/hybrid_model.h"
#include "solver/physical_models/mhd_model.h"
#include "solver/physical_models/physical_model.h"
#include "solver/solvers/solver.h"
#include "solver/messenger_registration.h"
#include "solver/level_initializer/level_initializer.h"
#include "solvers/solver_mhd.h"
#include "solver/solvers/solver_ppc.h"

#include "core/utilities/algorithm.h"



namespace PHARE
{
namespace solver
{
    struct LevelDescriptor
    {
        static const int NOT_SET  = -1;
        int modelIndex            = NOT_SET;
        int solverIndex           = NOT_SET;
        int resourcesManagerIndex = NOT_SET;
        std::string messengerName;
    };


    inline bool isRootLevel(int levelNumber) { return levelNumber == 0; }


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
    template<typename MessengerFactory, typename LevelnitializerFactory, typename AMR_Types>
    class MultiPhysicsIntegrator : public SAMRAI::mesh::StandardTagAndInitStrategy,
                                   public SAMRAI::algs::TimeRefinementLevelStrategy
    {
    public:
        static constexpr auto dimension = MessengerFactory::dimension;

        // model comes with its variables already registered to the manager system
        MultiPhysicsIntegrator(int nbrOfLevels)
            : nbrOfLevels_{nbrOfLevels}
            , levelDescriptors_(nbrOfLevels)

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




        std::string solverName(int iLevel) const { return getSolver_(iLevel).name(); }




        std::string modelName(int iLevel) const { return getModel_(iLevel).name(); }




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
        virtual void
        initializeLevelData(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                            const int levelNumber, const double initDataTime,
                            const bool canBeRefined, const bool initialTime,
                            const std::shared_ptr<SAMRAI::hier::PatchLevel>& oldLevel
                            = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                            const bool allocateData = true) override
        {
            auto& model            = getModel_(levelNumber);
            auto& solver           = getSolver_(levelNumber);
            auto& messenger        = getMessengerWithCoarser_(levelNumber);
            auto& levelInitializer = getLevelInitializer(model.name());

            const bool isRegridding = oldLevel != nullptr;
            auto level              = hierarchy->getPatchLevel(levelNumber);


            if (allocateData)
            {
                for (auto patch : *level)
                {
                    model.allocate(*patch, initDataTime);
                    solver.allocate(model, *patch, initDataTime);
                    messenger.allocate(*patch, initDataTime);
                }
            }

            messenger.registerLevel(hierarchy, levelNumber);

            levelInitializer.initialize(hierarchy, levelNumber, oldLevel, model, messenger,
                                        initDataTime, isRegridding);
        }



        virtual void
        resetHierarchyConfiguration(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                    const int coarsestLevel, const int finestLevel) override
        {
        }



        void applyGradientDetector(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                   const int levelNumber, const double error_data_time,
                                   const int tag_index, const bool initialTime,
                                   const bool usesRichardsonExtrapolationToo) override
        {
        }


        // -----------------------------------------------------------------------------------------------
        //
        //                          SAMRAI TimeRefinementLevelStrategy interface
        //
        // -----------------------------------------------------------------------------------------------



        virtual void initializeLevelIntegrator(
            const std::shared_ptr<SAMRAI::mesh::GriddingAlgorithmStrategy>& griddingAlg) override
        {
        }

        virtual double getLevelDt(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                  const double dtTime, const bool initialTime) override
        {
            return dtTime;
        }


        virtual double getMaxFinerLevelDt(const int finerLevelNumber, const double coarseDt,
                                          const SAMRAI::hier::IntVector& ratio) override
        {
            return std::pow(coarseDt, ratio.max());
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
        virtual double advanceLevel(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                    const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                    const double currentTime, const double newTime,
                                    const bool firstStep, const bool lastStep,
                                    const bool regridAdvance = false) override
        {
            if (regridAdvance)
                throw std::runtime_error("Error - regridAdvance must be False and is True");


            auto iLevel       = level->getLevelNumber();
            auto& solver      = getSolver_(iLevel);
            auto& model       = getModel_(iLevel);
            auto& fromCoarser = getMessengerWithCoarser_(iLevel);


            if (firstStep)
            {
                fromCoarser.firstStep(model, *level, hierarchy, currentTime);
            }


            fromCoarser.prepareStep(model, *level);


            solver.advanceLevel(hierarchy, iLevel, model, fromCoarser, currentTime, newTime);


            if (lastStep)
            {
                fromCoarser.lastStep(model, *level);
            }


            return newTime;
        }




        virtual void
        standardLevelSynchronization(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                     const int coarsestLevel, const int finestLevel,
                                     const double syncTime,
                                     const std::vector<double>& oldTimes) override
        {
            // TODO use messengers to sync with coarser
            auto& toCoarser = getMessengerWithCoarser_(finestLevel);
            auto level      = hierarchy->getPatchLevel(finestLevel);
            toCoarser.synchronize(*level);
        }




        virtual void
        synchronizeNewLevels(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                             const int coarsestLevel, const int finestLevel, const double syncTime,
                             const bool initialTime) override
        {
        }


        virtual void resetTimeDependentData(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                            const double newTime, const bool canBeRefined) override
        {
        }

        virtual void
        resetDataToPreadvanceState(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level) override
        {
        }

        virtual bool usingRefinedTimestepping() const override { return true; }




    private:
        int nbrOfLevels_;
        using IMessengerT       = amr::IMessenger<IPhysicalModel<AMR_Types>>;
        using LevelInitializerT = LevelInitializer<AMR_Types>;
        std::vector<LevelDescriptor> levelDescriptors_;
        std::vector<std::unique_ptr<ISolver<AMR_Types>>> solvers_;
        std::vector<std::shared_ptr<IPhysicalModel<AMR_Types>>> models_;
        std::map<std::string, std::unique_ptr<IMessengerT>> messengers_;
        std::map<std::string, std::unique_ptr<LevelInitializerT>> levelInitializers_;


        bool validLevelRange_(int coarsestLevel, int finestLevel)
        {
            if (coarsestLevel < 0 || finestLevel >= nbrOfLevels_ || finestLevel < coarsestLevel)
            {
                return false;
            }
            return true;
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




        void addModel_(std::shared_ptr<IPhysicalModel<AMR_Types>> model, int coarsestLevel,
                       int finestLevel)
        {
            if (core::notIn(model, models_))
            {
                levelInitializers_[model->name()] = LevelnitializerFactory::create(model->name());
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



        LevelInitializerT& getLevelInitializer(std::string modelName)
        {
            assert(levelInitializers_.count(modelName));
            return *levelInitializers_[modelName];
        }



        IMessengerT& getMessengerWithCoarser_(int iLevel)
        {
            auto& descriptor = levelDescriptors_[iLevel];
            auto messenger   = messengers_[descriptor.messengerName].get();
            auto s           = messenger->name();

            if (messenger)
                return *messenger;
            else
                throw std::runtime_error("no found messenger");
        }
    };
} // namespace solver
} // namespace PHARE
#endif
