
#ifndef PHARE_MULTIPHYSICS_INTEGRATOR_H
#define PHARE_MULTIPHYSICS_INTEGRATOR_H

#include <map>
#include <memory>
#include <type_traits>
#include <vector>

#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#include <SAMRAI/mesh/StandardTagAndInitStrategy.h>

#include "evolution/solvers/solver.h"
#include "evolution/solvers/solver_mhd.h"
#include "evolution/solvers/solver_ppc.h"

#include "physical_models/hybrid_model.h"
#include "physical_models/mhd_model.h"
#include "physical_models/physical_model.h"

#include "evolution/transactions/hybrid_transaction.h"
#include "evolution/transactions/mhd_transaction.h"
#include "evolution/transactions/transaction.h"
#include "evolution/transactions/transaction_initializer.h"
#include "evolution/transactions/transaction_manager.h"


#include "utilities/algorithm.h"



namespace PHARE
{
struct LevelDescriptor
{
    static const int NOT_SET = -1;
    int modelIndex           = NOT_SET;
    int solverIndex          = NOT_SET;
    // std::unique_ptr<ITransaction> fromCoarser;
    std::string transactionName;
};



template<typename TransactionFactory>
class MultiPhysicsIntegrator : public SAMRAI::mesh::StandardTagAndInitStrategy,
                               public SAMRAI::algs::TimeRefinementLevelStrategy
{
public:
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


    auto nbrOfLevels() const { return nbrOfLevels_; }


    /**
     * @brief registerModel registers the model to the multiphysics integrator for a given level
     * range. The level index for the coarsest and finest must be greater or equal to zero, less
     * than the nbrOfLevels().
     *
     * Once this method is called, the MultiPhysicsIntegrator will use this model for all levels in
     * the given range.
     * @param coarsestLevel is the index of the coarsest level using the model
     * @param finestLevel is the index of the finest level using the model. finestLevel >
     * coarsestModel
     * @param model the model to be registered to the MultiphysicsIntegrator. The model must not
     * have been registered already.
     */
    void registerModel(int coarsestLevel, int finestLevel, std::shared_ptr<PhysicalModel> model)
    {
        if (!validLevelRange_(coarsestLevel, finestLevel))
        {
            throw std::runtime_error("invalid range level");
        }

        for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
        {
            if (levelDescriptors_[iLevel].modelIndex != LevelDescriptor::NOT_SET)
            {
                throw std::runtime_error("error - model already set for level "
                                         + std::to_string(iLevel));
            }
        }

        if (notIn(model, models_))
        {
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




    /**
     * @brief registerAndInitSolver registers and initialize a given solver to the
     * MultiphysicsIntegrator
     *
     * Once this method is called :
     *
     * - the MultiPhysicsIntegrator will use this solver to advance the
     * hierarchy for levels between coarsestLevel and finestLevel (included). The level index for
     * the coarsest and finest must be greater or equal to zero, less than the nbrOfLevels().
     *
     * - the solver will have registered its variables to the ResourcesManager of the model in the
     * given range.
     *
     * @param coarsestLevel is the index of the coarsest level using the model
     * @param finestLevel is the index of the finest level using the model. finestLevel >
     * coarsestModel
     * @param solver is the ISolver to register to the MultiPhysicsIntegrator. This solver must be
     * compatible with all models in the given range and not yet registered.
     */
    void registerAndInitSolver(int coarsestLevel, int finestLevel, std::unique_ptr<ISolver> solver)
    {
        if (!validLevelRange_(coarsestLevel, finestLevel))
        {
            throw std::runtime_error("invalid level range");
        }


        for (auto iLevel = coarsestLevel; iLevel <= finestLevel; ++iLevel)
        {
            if (auto& model = getModel_(iLevel); !areCompatible(model, *solver))
            {
                throw std::runtime_error(model.name() + " is not compatible with " + solver->name()
                                         + ", (expecting " + solver->modelName() + ")");
            }
        }


        auto& model = getModel_(coarsestLevel);
        solver->registerResources(model);


        if (notIn(solver, solvers_))
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



    /**
     * @brief registerAndSetupTransactions registers and setup the transactions for the hierarchy
     *
     * This method create all necessary transactions for levels in the hierarchy to exchange data
     * knowing the sequence of models. This method must thus be called *after* the registerModel().
     *
     *
     * @param transactionFactory is used to create the appropriate transaction
     */
    void registerAndSetupTransactions(TransactionFactory& transactionFactory)
    {
        // for each pair of level register the correct transaction
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

            registerTransaction_(transactionFactory, coarseModel, fineModel, iLevel);
        }

        // now setup all transactions we've just created


        auto doSetup = [this](auto iLevel, auto& transaction) {
            auto coarseLevelNumber = iLevel - 1;
            auto fineLevelNumber   = iLevel;

            if (iLevel == 0)
            {
                coarseLevelNumber = fineLevelNumber;
            }

            auto& coarseModel = getModel_(coarseLevelNumber);
            auto& fineModel   = getModel_(fineLevelNumber);
            auto& solver      = getSolver_(iLevel);

            TransactionInitializer::setup(transaction, coarseModel, fineModel, solver);
        };



        auto transactionName = levelDescriptors_[0].transactionName;
        doSetup(0, getTransactionWithCoarser_(0));

        for (auto iLevel = 1; iLevel < nbrOfLevels_; ++iLevel)
        {
            if (transactionName != levelDescriptors_[iLevel].transactionName)
            {
                transactionName = levelDescriptors_[iLevel].transactionName;
                doSetup(iLevel, getTransactionWithCoarser_(iLevel));
            }
        }
    }




    std::string solverName(int iLevel) const { return getSolver_(iLevel).name(); }


    std::string modelName(int iLevel) const { return getModel_(iLevel).name(); }


    std::string transactionName(int iLevel)
    {
        auto& transaction = getTransactionWithCoarser_(iLevel);
        return transaction.name();
    }



    // -----------------------------------------------------------------------------------------------
    //
    //                          SAMRAI StandardTagAndInitStrategy interface
    //
    // -----------------------------------------------------------------------------------------------


    virtual void initializeLevelData(const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                     const int levelNumber, const double initDataTime,
                                     const bool canBeRefined, const bool initialTime,
                                     const std::shared_ptr<SAMRAI::hier::PatchLevel>& oldLevel
                                     = std::shared_ptr<SAMRAI::hier::PatchLevel>(),
                                     const bool allocateData = true) override
    {
        auto& model       = getModel_(levelNumber);
        auto& solver      = getSolver_(levelNumber);
        auto& transaction = getTransactionWithCoarser_(levelNumber);


        // here we need to allocate PatchDatas for
        // - the model
        // - the solver
        if (allocateData)
        {
            auto level = hierarchy->getPatchLevel(levelNumber);
            for (auto patch : *level)
            {
                model.allocate(*patch, initDataTime);
                solver.allocate(model, *patch, initDataTime);
                transaction.allocate(*patch, initDataTime);
            }
            // TODO: transactions may need to allocate data in initializeLevelData too
        }



        transaction.setLevel(hierarchy, levelNumber);

        // on est en train de changer la hierarchy soit en créant un nouveau niveau (finest)
        // soit en regriddant un niveau.
        // du coup tous les schedules concernant ce niveau sont devenus invalides
        // en gros on doit refaire les memes en passant le pointeur sur le newLevel

        if (oldLevel)
        {
            // in case of a regrid we need to make a bunch of temporary regriding schedules
            // using the init algorithms and actually perform the .fillData() for all of them
            transaction.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
        }


        else // we're creating a brand new finest level in the hierarchy
        {
            if (levelNumber == 0)
            {
                // here we are either starting the simulation and building the root level
                // or building from a restart
                // either way it's not our business here, and we use the initializer
                // we where kindy given

                // initializer.init(model);
            }
            else
            {
                transaction.initLevel(levelNumber, initDataTime);
            }
        }
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
     * In this method, the MultiPhysicsIntegrator needs to get the model, solver and transaction
     * necessary to advance the given level. It then forwards the call to the solver's advanceLevel
     * method, passing it the model and transaction, among others.
     *
     * If it is the first step of the subcycling the Transaction may have something to do before
     * calling the solver's advanceLevel(). Typically it may need to grab coarser level's data at
     * the coarser future time so that all subsequent subcycles may use time interpolation without
     * involving communications. This is done by calling transaction.firstStep()
     *
     *
     * At the last step of the subcycle, the Transaction may also need to perform some actions, like
     * working on its internal data for instance. transaction.lastStep()
     */
    virtual double advanceLevel(const std::shared_ptr<SAMRAI::hier::PatchLevel>& level,
                                const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                const double currentTime, const double newTime,
                                const bool firstStep, const bool lastStep,
                                const bool regridAdvance = false) override
    {
        //
        if (regridAdvance)
            throw std::runtime_error("Error - regridAdvance must be False and is True");


        // TODO transaction needs to copy stuff from the model into internal 'old' variables


        auto iLevel = level->getLevelNumber();

        auto& solver      = getSolver_(iLevel);
        auto& model       = getModel_(iLevel);
        auto& fromCoarser = getTransactionWithCoarser_(iLevel);


        if (firstStep)
        {
            fromCoarser.firstStep(model);
        }

        // TODO give firstStep and lastStep bools to the transaction
        // and levelNumber, so that it knows which schedules to apply later on.
        // fromCoarser.setStepInfo(firstStep, lastStep, levelNumber)

        // solver msut have a view on the model from its init
        solver.advanceLevel(hierarchy, iLevel, model, fromCoarser, currentTime, newTime);


        if (lastStep)
        {
            fromCoarser.lastStep(model);
        }


        return newTime;
    }




    virtual void standardLevelSynchronization(
        const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy, const int coarsestLevel,
        const int finestLevel, const double syncTime, const std::vector<double>& oldTimes) override
    {
        // TODO use transactions to sync with coarser
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
    std::vector<LevelDescriptor> levelDescriptors_;
    std::vector<std::unique_ptr<ISolver>> solvers_;
    std::vector<std::shared_ptr<PhysicalModel>> models_;
    std::map<std::string, std::unique_ptr<ITransaction>> transactions_;


    bool validLevelRange_(int coarsestLevel, int finestLevel)
    {
        if (coarsestLevel < 0 || finestLevel >= nbrOfLevels_ || finestLevel <= coarsestLevel)
        {
            return false;
        }
        return true;
    }


    void registerTransaction_(TransactionFactory const& transactions,
                              PhysicalModel const& coarseModel, PhysicalModel const& fineModel,
                              int iLevel)
    {
        if (auto transactionName = transactions.name(coarseModel, fineModel); transactionName)
        {
            auto foundTransaction = transactions_.find(*transactionName);
            if (foundTransaction == transactions_.end())
            {
                transactions_[*transactionName] = std::move(
                    transactions.create(*transactionName, coarseModel, fineModel, iLevel));
            }

            levelDescriptors_[iLevel].transactionName = *transactionName;
        }
        else
        {
            throw std::runtime_error("No viable transaction found");
        }
    }



    ISolver& getSolver_(int iLevel)
    {
        return const_cast<ISolver&>(
            const_cast<std::remove_pointer_t<decltype(this)> const*>(this)->getSolver_(iLevel));
    }


    ISolver const& getSolver_(int iLevel) const
    {
        auto& descriptor = levelDescriptors_[iLevel];
        return *solvers_[descriptor.solverIndex];
    }


    PhysicalModel& getModel_(int iLevel)
    {
        return const_cast<PhysicalModel&>(
            const_cast<std::remove_pointer_t<decltype(this)> const*>(this)->getModel_(iLevel));
    }

    PhysicalModel const& getModel_(int iLevel) const
    {
        auto& descriptor = levelDescriptors_[iLevel];
        return *models_[descriptor.modelIndex];
    }

    ITransaction& getTransactionWithCoarser_(int iLevel)
    {
        auto& descriptor = levelDescriptors_[iLevel];
        auto transaction = transactions_[descriptor.transactionName].get();
        auto s           = transaction->name();

        if (transaction)
            return *transaction;
        else
            throw std::runtime_error("no found transaction");
    }


}; // namespace PHARE
} // namespace PHARE
#endif
