
#ifndef PHARE_HYBRID_HYBRID_TRANSACTION_STRATEGY_H
#define PHARE_HYBRID_HYBRID_TRANSACTION_STRATEGY_H

#include "data/field/coarsening/field_data_coarsen.h"
#include "data/field/refine/field_refine_operator.h"
#include "data/field/time_interpolate/field_linear_time_interpolate.h"
#include "data/particles/refine/particles_data_split.h"
#include "evolution/transactions/hybrid_transaction_info.h"
#include "evolution/transactions/hybrid_transaction_strategy.h"
#include "physical_models/physical_model.h"
#include "tools/resources_manager_utilities.h"

#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>


#include <utility>


namespace PHARE
{
/** \brief An HybridTransaction purpose is to manage transaction from an to : the electric field,
 * magnetic field, and the ions.
 *
 *
 * TODO: we would also want to add rho soon
 * TODO: we have a lot of get* send* maybe we should use some enum instead,so that we can make
 * combinaison with other parameter easier
 *
 */
template<typename HybridModel>
class HybridHybridTransactionStrategy : public HybridTransactionStrategy<HybridModel>
{
    using IonsT     = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

public:
    static const std::string stratName;


    HybridHybridTransactionStrategy(
        std::shared_ptr<typename HybridModel::resources_manager_type> manager, int const firstLevel)
        : HybridTransactionStrategy<HybridModel>{stratName}
        , resourcesManager_{std::move(manager)}
        , firstLevel_{firstLevel}
    {
        resourcesManager_->registerResources(EM_old_.E);
        resourcesManager_->registerResources(EM_old_.B);
    }

    virtual ~HybridHybridTransactionStrategy() = default;



    /**
     * @brief allocate the transaction strategy internal variables to the model resourceManager
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
    {
        // hybModel.resourcesManager->allocate(EM_old_.E, patch, allocateTime);
        // hybModel.resourcesManager->allocate(EM_old_.B, patch, allocateTime);
        // hybModel.resourcesManager->allocate(EM_old_, patch, allocateTime);
        resourcesManager_->allocate(EM_old_, patch, allocateTime);
    }



    /**
     * @brief setup creates all SAMRAI algorithms to communicate data involved in a transaction
     * between the coarse and fine levels.
     *
     * This method creates the SAMRAI algorithms for communications associated between pairs of
     * variables. The function does not create the SAMRAI schedules since they depend on the levels
     */
    virtual void setup(std::unique_ptr<ITransactionInfo> fromCoarserInfo,
                       [[maybe_unused]] std::unique_ptr<ITransactionInfo> fromFinerInfo) override
    {
        std::unique_ptr<HybridTransactionInfo> hybridInfo{
            dynamic_cast<HybridTransactionInfo*>(fromCoarserInfo.release())};


        // TODO here we need to use both coarse and fine since solvers may differ

        setupEandBGhostTransactions_(hybridInfo);
        setupEandBInitTransactions_(hybridInfo);

        // on a tous les algos pour remplir les ghosts de
        // - model E, B vers model E, B
        // - model E, B vers solver list of E_internal and list of B_internal

        // TODO setup particle transactions
    }



    /**
     * @brief setLevel creates SAMRAI schedules for all variables relevant for hybrid to hybrid
     * communications, whether they are from model, solver or transaction internals
     *
     * The method takes maps associating quantityName to the pair (algo, map of schedules)
     * for each algo, it creates a bunch of associated schedules and store them in the vector
     */
    virtual void setLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                          int const levelNumber) override
    {
        // update the schedules for ghost communications of the electric and magnetic field

        auto level = hierarchy->getPatchLevel(levelNumber);
        createGhostSchedules_(magneticGhostsRefine_, hierarchy, level);
        createGhostSchedules_(electricGhostsRefine_, hierarchy, level);

        // now update electric and magnetic initialization
        // root level is not initialized with a schedule using coarser level data
        // so we don't create these schedules if root level
        if (levelNumber != 0)
        {
            createInitSchedules_(magneticInitRefine_, hierarchy, level);
            createInitSchedules_(electricInitRefine_, hierarchy, level);
        }
    }



    /**
     * @brief regrid performs the regriding communications for Hybrid to Hybrid transactions
     *
     * This routine must create and execute schedules for:
     *
     * - filling the model magnetic field
     */
    virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        const int levelNumber,
                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                        double const initDataTime) override
    {
        // create temporary schedule from algorithms

        auto doRegrid = [&hierarchy, levelNumber, oldLevel, initDataTime](auto& refiners) //
        {
            for (auto& [key, refiner] : refiners)
            {
                auto& algo = refiner.algo;

                // here 'nullptr' is for 'oldlevel' which is always nullptr in this function
                // the regriding schedule for which oldptr is not nullptr is handled in another
                // function
                auto const& level = hierarchy->getPatchLevel(levelNumber);

                auto schedule = algo->createSchedule(
                    level, oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy);

                schedule->fillData(initDataTime);
            }
        };

        // root level is not initialized with a schedule using coarser level data
        // so we don't create these schedules if root level

        doRegrid(magneticInitRefine_);
        doRegrid(electricInitRefine_);
    }




    virtual std::string fineModelName() const override { return HybridModel::model_name; }

    virtual std::string coarseModelName() const override { return HybridModel::model_name; }




    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromCoarser() override
    {
        return std::make_unique<HybridTransactionInfo>();
    }

    virtual std::unique_ptr<ITransactionInfo> emptyInfoFromFiner() override
    {
        return std::make_unique<HybridTransactionInfo>();
    }




    virtual void initLevel(int const levelNumber, double const initDataTime) const override
    {
        applyInitSchedules_(levelNumber, initDataTime, magneticInitRefine_);
        applyInitSchedules_(levelNumber, initDataTime, electricInitRefine_);
    }



#if 0
    HybridHybridTransactionStrategy(
        HybridState const &hybridModel, ResourcesManager &resources,
        std::shared_ptr<SAMRAI::hier::RefineOperator> const &fieldRefineOp,
        std::shared_ptr<SAMRAI::hier::CoarsenOperator> const &fieldCoarsenOp,
        std::shared_ptr<SAMRAI::hier::RefineOperator> const &particlesRefineOp,
        int const startLevel)
        : model_{hybridModel}
        , electricNew_{model_.electromag.E.name() + "_new", HybridQuantity::Vector{}}
        , magneticNew_{model_.electromag.B.name() + "_new", HybridQuantity::Vector{}}
        , resources_{resources}
        , fieldRefineOp_{fieldRefineOp}
        , fieldCoarsenOp_{fieldCoarsenOp}
        , particlesRefineOp_{particlesRefineOp}
        , startLevel_{startLevel}
        , dimension_{SAMRAI::tbox::Dimension{HybridState::dimension}}
    {
        resources_.registerResources(electricNew_);
        resources_.registerResources(magneticNew_);
    }
    // here we have two comportement:
    // either we already have a schedule and in this case we simply use it
    // or we don't have one and we create the schedule
    // note that we will have the current level , and the hierarchy access
    // also we know the level number of the source and the destination
    //
#endif


    virtual void fillMagneticGhosts(VecFieldT& B, int const levelNumber,
                                    double const fillTime) override
    {
        std::cout << "perform the magnetic ghost fill\n";

        auto name = B.name();

        auto algoAndSchedulesItem = magneticGhostsRefine_.find(name);
        if (algoAndSchedulesItem != std::end(magneticGhostsRefine_))
        {
            auto& schedules   = algoAndSchedulesItem->second.schedules;
            auto scheduleItem = schedules.find(levelNumber);
            if (scheduleItem != std::end(schedules))
            {
                auto& schedule = scheduleItem->second;
                schedule->fillData(fillTime);
            }
            else
            {
                throw std::runtime_error(
                    "Error - cannot find schedules for Magnetic Ghosts filling");
            }
        }
        else
        {
            throw std::runtime_error("Error - cannot find algo for Magnetic Ghosts filling");
        }

        // auto name      = B.name();
        // auto algoID    = getAlgo(name);
        // auto schedules = ghostSchedules(algoID, name);

        // getSchedules()
        // for (auto& schedule : schedules_)
        //{
        // schedule.fillData();
        //}

#if 0
        showAction_<withTemporal, Enum, fillType>();

        auto findIt = algorithmsMagneticIn_.find(B.name());
        if (findIt != std::end(algorithmsMagneticIn_))
        {
            auto &algo = findIt->second;
            // here it all depend on fillType
            if (algo->refineSchedule.size() <= static_cast<std::size_t>(relativeLevel_))
            {
                algo->refineSchedule.push_back(algo->refine.createSchedule(
                    *currentLevel_, currentLevelNumber_ - 1, *hierarchy_));
            }

            algo->refineSchedule[relativeLevel_]->fillData(fillTime);
        }
#endif
    }

    virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber,
                                    double const fillTime) override
    {
        std::cout << "perform the electric ghost fill\n";

        auto name = E.name();

        auto algoAndSchedulesItem = electricGhostsRefine_.find(name);
        if (algoAndSchedulesItem != std::end(electricGhostsRefine_))
        {
            auto& schedules   = algoAndSchedulesItem->second.schedules;
            auto scheduleItem = schedules.find(levelNumber);
            if (scheduleItem != std::end(schedules))
            {
                auto& schedule = scheduleItem->second;
                schedule->fillData(fillTime);
            }
            else
            {
                throw std::runtime_error(
                    "Error - cannot find schedules for electric Ghosts filling");
            }
        }
        else
        {
            throw std::runtime_error("Error - cannot find algo for electric Ghosts filling");
        }
    }



    virtual void fillIonGhostParticles(IonsT& ions, int const levelNumber,
                                       double const fillTime) override
    {
        std::cout << "perform the ghost particle fill\n";

        // tout ce qu'on fait ici c'est d'laler récup les ghosts des patchs
        // voisins sur le meme level, pour pouvoir projeter la densite/flux
        // manquants sur les points bleus "dans" le level
        // et puis avoir les ghosts particuels ready pour le next step
        for (auto& pop : ions)
        {
            // ici on appelle le schedule qui va chercher sur le meme niveau
            // le particules des patchs voisins qui sont dans ma ghost zone
            // pour aller dans mes ghosts
            // algo.registerRefine(pop_particleArray_ID // dest
            //  , pop_particleArray_ID // src
            //  , pop_particleArray_ID // scratch
            // nullptr )// pas de refineOp car same level

            // on a trouve l'algo 'algo'
            // le schedule  a deja ete cree (dans initiliazeLevelData)
            // schedule =  algo.schedules[relatveLevel]
            // schedule.fillData(fillTime);
        }
    }


    virtual void fillIonMomentGhosts(IonsT& ions, int const levelNumber,
                                     double const fillTime) override
    {
        std::cout << "perform the moments ghosts fill\n";

        // recuperer coarse2fine2 de la transction depuis le coarser level
        // et calculer alpha*coarse2fine2 + (1-alpha)*coarse2fine1 sur les
        // density et flux
    }




    // synchronization/coarsening methods
    virtual void syncMagnetic(VecFieldT& B) override
    {
        //
        std::cout << "perform coarse magnetic sync to the next coarser level\n";
    }


    virtual void syncElectric(VecFieldT& E) override
    {
        //
        std::cout << "perform coarse electric sync to the next coarser level\n";
    }


    virtual void syncIonMoments(IonsT& ions) override
    {
        //
        std::cout << "perform coarse moments sync to the next coarser level\n";
    }



private:
    /**
     * @brief The Refiner struct encapsulate the algorithm and its associated
     * schedules
     *
     * We have several of those object, one per level, that is used to retrieve which schedule to
     * use for a given transaction communication, and which algorithm to use to re-create schedules
     * when initializing a level
     */
    struct Refiner
    {
        Refiner()
            : algo{std::make_unique<SAMRAI::xfer::RefineAlgorithm>()}
        {
        }

        std::unique_ptr<SAMRAI::xfer::RefineAlgorithm> algo; // this part is set in setup

        // this part is created in initializeLevelData()
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> schedules;
    };

    using Refiners = std::map<std::string, Refiner>;



    /**
     * @brief setupEandBGhostTransactions_ creates the SAMRAI algorithms to transfer
     * electromagnetic fields on ghost nodes of level patches and register all refine operations.
     *
     * For each variable to register, we use the registerRefine overload that allows time
     * interpolation.
     *
     * The operation must use 'src_told' and 'src_tnew' versions of each field and
     * time/space interpolate that onto the ghost nodes of the same variable onto the destination
     * level patches.
     *
     * The schedules originating from these algorithms will be executed in the solver
     * to fill the ghosts on a given level. At this point, the 'src_tnew' IDs are just model
     * electric field IDs since when advancing a level, the next model of the next coarser level is
     * already advanced. The "src_told" IDs however, is found in the Transaction itself as the
     * EM_old_ electromagnetic field This field is a copy of the model electromagnetic field at t=n
     * *before* the solver advanced the model to t=n+1 on the coarser level.
     *
     *
     * Each operation thus takes:
     *  - src       : variable IDs in model (not sure what it is used for here)
     *  - dest      : variable IDs in the TransactionInfo.ghostXXXXnames (could be from model, could
     * be from solver)
     *  - src_told  : variable IDs in Transaction
     *  - src_tnew  : variable IDs in model
     *
     *
     * The SAMRAI algorithms created by this function are stored in RefineAlgosAndSchedule
     * objects, themselves stored in an appropriate map with a key provided by the
     * HybridTransactionInfo
     *
     */
    void setupEandBGhostTransactions_(std::unique_ptr<HybridTransactionInfo> const& info)
    {
        auto const& Eold = EM_old_.E;
        auto const& Bold = EM_old_.B;


        // register a refine operation for all components in src/dest componentNames
        // using the components of the vecField_old as the old source
        auto registerRefine = [this](auto const& srcComponentNames, auto const& destComponentNames,
                                     auto const& names_old) //
        {
            Refiner algsAndSched;
            algsAndSched.algo = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();


            // there is one refine operation to register per component of the VecField
            // since they are each a specific variable
            for (auto componentIndex = 0u; componentIndex < 3; ++componentIndex)
            {
                auto src_id  = resourcesManager_->getID(srcComponentNames[componentIndex]);
                auto dest_id = resourcesManager_->getID(destComponentNames[componentIndex]);
                auto old_id  = resourcesManager_->getID(names_old[componentIndex]);

                if (src_id && dest_id && old_id)
                {
                    // dest, src, old, new, scratch
                    algsAndSched.algo->registerRefine(
                        *dest_id, // dest
                        *src_id,  // source at same time
                        *old_id,  // source at past time (for time interp)
                        *src_id,  // source at future time (for time interp)
                        *dest_id, // scratch
                        fieldRefineOp_, fieldTimeOp_);
                }
            }
            return algsAndSched;
        };



        auto algosToMap = [&registerRefine](VecFieldNames const& modelVec,
                                            std::vector<VecFieldNames> const& ghostVec,
                                            auto const& oldVecNames, auto& map) //
        {
            auto nbrVariables = ghostVec.size();
            for (auto i = 0u; i < nbrVariables; ++i)
            {
                std::array<std::string, 3> destComponentNames{
                    {ghostVec[i].xName, ghostVec[i].yName, ghostVec[i].zName}};

                std::array<std::string, 3> srcComponentNames{
                    {modelVec.xName, modelVec.yName, modelVec.zName}};

                map[ghostVec[i].vecName]
                    = registerRefine(srcComponentNames, destComponentNames, oldVecNames);
            }
        };


        algosToMap(info->modelMagnetic, info->ghostMagnetic, extractNames(Bold),
                   magneticGhostsRefine_);
        algosToMap(info->modelElectric, info->ghostElectric, extractNames(Eold),
                   electricGhostsRefine_);
    }




    void setupEandBInitTransactions_(std::unique_ptr<HybridTransactionInfo> const& info)
    {
        Refiner initB;
        Refiner initE;

        auto registerRefine = [this](auto& srcComponentNames, auto& destComponentNames,
                                     auto& algo) //
        {
            for (auto componentIndex = 0u; componentIndex < 3; ++componentIndex)
            {
                auto src_id  = resourcesManager_->getID(srcComponentNames[componentIndex]);
                auto dest_id = resourcesManager_->getID(destComponentNames[componentIndex]);

                auto id1 = *src_id;
                auto id2 = *dest_id;
                if (src_id && dest_id)
                {
                    // dest, src, old, new, scratch
                    algo->registerRefine(*dest_id, // dest
                                         *src_id,  // source at same time
                                         *dest_id, // scratch
                                         fieldRefineOp_);
                }
                else
                {
                    throw std::runtime_error("invalid IDs");
                }
            }
        };


        auto algosToMap = [&registerRefine](std::vector<VecFieldNames> const& names, auto& map,
                                            auto& algsAndSchedules) //
        {
            auto nbrVariables = names.size();
            for (auto i = 0u; i < nbrVariables; ++i)
            {
                std::array<std::string, 3> componentNames{
                    {names[i].xName, names[i].yName, names[i].zName}};

                registerRefine(componentNames, componentNames, algsAndSchedules.algo);
                map[names[i].vecName] = std::move(algsAndSchedules);
            }
        };


        algosToMap(info->initMagnetic, magneticInitRefine_, initB);
        algosToMap(info->initElectric, electricInitRefine_, initE);
    }




    void createGhostSchedules_(Refiners& refiners,
                               std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
    {
        for (auto& [key, refiner] : refiners)
        {
            auto& algo    = refiner.algo;
            auto schedule = algo->createSchedule(level, level->getNextCoarserHierarchyLevelNumber(),
                                                 hierarchy);

            refiner.schedules[level->getLevelNumber()] = std::move(schedule);
        }
    }




    void createInitSchedules_(Refiners& refiners,
                              std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        for (auto& [key, refiner] : refiners)
        {
            auto& algo       = refiner.algo;
            auto levelNumber = level->getLevelNumber();

            // note that here we must take that createsSchedule() overload and put nullptr as src
            // since we want to take from coarser level everywhere.
            // using the createSchedule overload that takes level, next_coarser_level only
            // would result in interior ghost nodes to be filled with interior of neighbor patches
            // but there is nothing there.
            refiner.schedules[levelNumber]
                = algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy);
        }
    }



    void applyInitSchedules_(int levelNumber, double initDataTime, Refiners const& refiners) const
    {
        for (auto& [key, refiner] : refiners)
        {
            if (refiner.algo == nullptr)
            {
                throw std::runtime_error("Algorithm is nullptr");
            }

            auto schedule = refiner.schedules.find(levelNumber);
            if (schedule != std::end(refiner.schedules))
            {
                schedule->second->fillData(initDataTime);
            }
            else
            {
                throw std::runtime_error("Error - schedule cannot be found for this level");
            }
        }
    }




    using Electromag = decltype(std::declval<HybridModel>().state.electromag);

    using gridlayout_type = typename HybridModel::gridLayout_type;

    using field_type =
        typename decltype(std::declval<HybridModel>().state.electromag)::vecfield_type::field_type;

    // field data refine op
    std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_{
        std::make_shared<FieldRefineOperator<gridlayout_type, field_type>>()};

    // field data time op
    std::shared_ptr<FieldLinearTimeInterpolate<gridlayout_type, field_type>> fieldTimeOp_{
        std::make_shared<FieldLinearTimeInterpolate<gridlayout_type, field_type>>()};



    //! keeps a copy of the model electromagnetic field at t=n
    Electromag EM_old_{"EM_old"}; // TODO needs to be allocated somewhere and
                                  // updated to t=n before advanceLevel()


    std::shared_ptr<typename HybridModel::resources_manager_type> resourcesManager_;


    int const firstLevel_;




    // the following maps store the algorithm and the associated schedules for different quantities
    // and different operations. The key of the map is the name of the quantity, e.g.
    // "HybridModel_EM_B" and the value will be the {algo, schedules[]}
    // there may be several maps for the same quantity if that quantity needs several algorithms
    // this is typically the case for the magnetic field for instance, which needs to spatially
    // interpolated at initialization phase (spatial only registerRefine()) AND to be spatial/time
    // interpolated at advance phase when filling ghosts


    //! stores the algo and its associated schedules for the getting ghost magnetic field
    Refiners magneticGhostsRefine_;

    //! stores the algo and its associated schedules for getting magnetic field at initialization
    Refiners magneticInitRefine_;

    // same as for the magnetic field
    Refiners electricGhostsRefine_;
    Refiners electricInitRefine_;


    // algo and schedule used to initialize domain particles
    // from coarser level using particleRefineOp<domain>
    Refiner particleInitRefine_;

    // keys : model particles (initialization and 2nd push), temporaryParticles (firstPush)
    Refiners particleGhostExchange_;

    // at first step of advance:
    // from temporaryParticle of coarseLevel to model PRA1 ( + PRA1 copy into PRA)
    // and from modelParticles of coarserLevel to model PRA2
    // the copy of PRA1 vector to PRA is done after the schedule in a PatchStrategy post truc
    // these are ran before solver->advanceLevel() in the MultiPhysics::advanceLevel()
    // with a method : transaction.firstStepOperation() or somthg like that...
    // solver->advanceLevel() when starting, has all PRAs set correctly

    // keys: PRA1, PRA2 , chosen by transaction
    Refiners particlePRA_;




#if 0
template<typename Variable, bool withTemporal, typename Enum, Enum fillType>
void getElectric(Variable &E, double fillTime, BooleanSelector<withTemporal>,
                 EnumSelector<Enum, fillType>)
{
    std::cout << "perform the electric fill\n";
    showAction_<withTemporal, Enum, fillType>();
}




template<typename Variable, bool withTemporal, typename Enum, Enum fillType>
void getIons(Variable &ions, BooleanSelector<withTemporal>, EnumSelector<Enum, fillType>)
{
    std::cout << "perform ions fill\n";
    showAction_<withTemporal, Enum, fillType>();
}


template<typename Variable, bool withTemporal, typename Enum, Enum fillType>
void sendMagnetic(Variable &B, BooleanSelector<withTemporal>, EnumSelector<Enum, fillType>)
{
    std::cout << "send magnetic to the next finer level\n";
    showAction_<withTemporal, Enum, fillType>();
}




template<typename Variable, bool withTemporal, typename Enum, Enum fillType>
void sendElectric(Variable &E, BooleanSelector<withTemporal>, EnumSelector<Enum, fillType>)
{
    std::cout << "send electric to the next finer level\n";
    showAction_<withTemporal, Enum, fillType>();
}
template<typename Variable, bool withTemporal, typename Enum, Enum fillType>
void sendIons(Variable &ions, BooleanSelector<withTemporal>, EnumSelector<Enum, fillType>)
{
    std::cout << "send ions to the next finer level\n";
    showAction_<withTemporal, Enum, fillType>();
}


template<typename Variable>
void syncMagnetic(Variable &B)
{
    std::cout << "perform coarse sync to the next coarser level\n";
}
template<typename Variable>
void syncElectric(Variable &E)
{
    std::cout << "perform coarse sync to the next coarser level\n";
}



// TODO : we have to take the temporal argument
// so we may have twice more algorithms maps
template<typename Variable, bool withTemporal>
void registerMagneticIn(Variable &B, BooleanSelector<withTemporal>)
{
    std::cout << "register a transaction from magnetic(coarse level) to variable in the "
                 "current level\n";

    showTemporalStatus_<withTemporal>();

    std::cout << "create algorithm : magnetic to " << B.name() << "\n";
    auto destID = resources_.getIDs(B);
    showVec_(destID);

    if constexpr (withTemporal)
    {
        // if we want a temporal sync, we first transfert the data in magneticNew_;
        // and then we time interpolate from model_.electromag.B to magneticNew_
        // into B, finnaly we sync at the same level

        // TODO this should be unique for electric to electric
        //      and magnetic to magnetic
        if (!alreadyRegisterMagneticNew_)
        {
            auto refineAlgInit = makeRefineAlgForField_(model_.electromag.B, magneticNew_);
            algorithmsMagneticIn_.try_emplace(magneticNew_.name(), refineAlgInit);
            alreadyRegisterMagneticNew_ = true;
        }

        // finnaly sync at the same level the destination
        auto refineAlg = makeRefineAlgForField_(B, B);
        algorithmsMagneticIn_.try_emplace(B.name(), refineAlg);
    }
    else
    {
        auto refineAlg = makeRefineAlgForField_(model_.electromag.B, B);
        algorithmsMagneticIn_.try_emplace(B.name(), refineAlg);
    }
}
template<typename Variable, bool withTemporal>
void registerElectricIn(Variable &E, BooleanSelector<withTemporal>)
{
    std::cout << "register a transaction from electric(coarse level) to variable in the "
                 "current level\n";

    showTemporalStatus_<withTemporal>();

    std::cout << "create algorithm : electric to " << E.name() << "\n";

    auto destID = resources_.getIDs(E);
    showVec_(destID);

    if constexpr (withTemporal)
    {
        if (!alreadyRegisterElectricNew_)
        {
            auto refineAlgInit = makeRefineAlgForField_(model_.electromag.E, electricNew_);

            algorithmsElectricIn_.try_emplace(electricNew_.name(), refineAlgInit);

            alreadyRegisterElectricNew_ = true;
        }


        auto refineAlg = makeRefineAlgForField_(E, E);

        algorithmsElectricIn_.try_emplace(E.name(), refineAlg);
    }
    else
    {
        auto refineAlg = makeRefineAlgForField_(model_.electromag.E, E);

        algorithmsElectricIn_.try_emplace(E.name(), refineAlg);
    }
}
template<typename Variable, bool withTemporal>
void registerIonsIn(Variable &ions, BooleanSelector<withTemporal>)
{
    std::cout << "register a transaction from ions(coarse level) to variable in the "
                 "current level\n";

    showTemporalStatus_<withTemporal>();

    for (auto const &pop : ions)
    {
        std::cout << "create algorithm : ionsPop to " << pop.name() << "\n";

        auto destID = resources_.getIDs(pop);
        showVec_(destID);

        algorithmsIonsIn_.try_emplace(pop.name(), std::make_shared<Algorithm>(dimension_));
    }
}

template<typename Variable, bool withTemporal>
void registerMagneticOut(Variable &B, BooleanSelector<withTemporal>)
{
    std::cout << "register a transaction from variable(current level) to magnetic"
                 "of next finer level\n";

    showTemporalStatus_<withTemporal>();
    std::cout << "create algorithm : " << B.name() << " to magnetic\n";

    auto destID = resources_.getIDs(B);
    showVec_(destID);

    auto refineAlg = makeRefineAlgForField_(B, model_.electromag.B);

    algorithmsMagneticOut_.try_emplace(B.name(), refineAlg);
}
template<typename Variable, bool withTemporal>
void registerElectricOut(Variable &E, BooleanSelector<withTemporal>)
{
    std::cout << "register a transaction from variable(current level) to electric"
                 "of next finer level\n";

    showTemporalStatus_<withTemporal>();
    std::cout << "create algorithm : " << E.name() << " to electric\n";

    auto destID = resources_.getIDs(E);
    showVec_(destID);

    auto refineAlg = makeRefineAlgForField_(E, model_.electromag.E);
    algorithmsElectricOut_.try_emplace(E.name(), refineAlg);
}
template<typename Variable, bool withTemporal>
void registerIonsOut(Variable &ions, BooleanSelector<withTemporal>)
{
    std::cout << "register a transaction from variable(current level) to ions"
                 "of next finer level\n";

    showTemporalStatus_<withTemporal>();

    for (auto const &pop : ions)
    {
        std::cout << "create algorithm : " << pop.name() << " to ionsPop\n";

        auto destID = resources_.getIDs(pop);
        showVec_(destID);

        algorithmsIonsOut_.try_emplace(pop.name(), std::make_shared<Algorithm>(dimension_));
    }
}

template<typename Variable>
void registerMagneticSync(Variable &B)
{
    std::cout << "register electric coarse transaction\n";
    std::cout << "create algorithm : " << B.name() << " to magnetic\n";

    auto destID = resources_.getIDs(B);
    showVec_(destID);

    auto coarseAlg = makeCoarseAlgForField_(B, model_.electromag.B);

    algorithmsMagneticSync_.try_emplace(B.name(), coarseAlg);
}
template<typename Variable>
void registerElectricSync(Variable &E)
{
    std::cout << "register magnetic coarse transaction\n";
    std::cout << "create algorithm : " << E.name() << " to electric\n";

    auto destID = resources_.getIDs(E);
    showVec_(destID);

    auto coarseAlg = makeCoarseAlgForField_(E, model_.electromag.E);

    algorithmsMagneticSync_.try_emplace(E.name(), coarseAlg);
}

void setStatus(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const *hierarchy,
               std::shared_ptr<SAMRAI::hier::PatchLevel> const *currentLevel,
               std::shared_ptr<SAMRAI::hier::PatchLevel> const *oldLevel)
{
    hierarchy_    = hierarchy;
    currentLevel_ = currentLevel;
    oldLevel_     = oldLevel;

    currentLevelNumber_ = (*currentLevel)->getLevelNumber();

    if ((*currentLevel)->inHierarchy())
    {
        relativeLevel_ = currentLevelNumber_ - startLevel_;
    }
    else
    {
        relativeLevel_ = 0;
    }
}

void updateSchedule(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const &hierarchy,
                    int const levelNumber)
{
    //
}


void init(SAMRAI::hier::Patch const &patch)
{
    resources_.allocate(electricNew_, patch);
    resources_.allocate(magneticNew_, patch);
}


private:
template<typename VariableSrc, typename VariableDest>
auto makeRefineAlgForField_(VariableSrc &src, VariableDest &dest)
{
    auto destID = resources_.getIDs(dest);

    auto srcID = resources_.getIDs(src);


    auto refineAlg = std::make_shared<Algorithm>(dimension_);

    refineAlg->refOperator = &fieldRefineOp_;


    SAMRAI::xfer::RefineAlgorithm algorithm{};
    for (auto iComponent = 0u; iComponent < NBR_COMPO; ++iComponent)
    {
        algorithm.registerRefine(destID[iComponent], srcID[iComponent], destID[iComponent],
                                 *refineAlg->refOperator);
    }

    return refineAlg;
}

template<typename VariableSrc, typename VariableDest>
auto makeCoarseAlgForField_(VariableSrc &src, VariableDest &dest)
{
    auto destID = resources_.getIDs(dest);

    auto srcID = resources_.getIDs(src);


    auto coarseAlg = std::make_shared<Algorithm>(dimension_);

    coarseAlg->coarseOperator = &fieldCoarsenOp_;

    SAMRAI::xfer::CoarsenAlgorithm algorithm{SAMRAI::tbox::Dimension{VariableSrc::dimension}};
    for (auto iComponent = 0u; iComponent < NBR_COMPO; ++iComponent)
    {
        algorithm.registerCoarsen(destID[iComponent], srcID[iComponent],
                                  *coarseAlg->coarseOperator);
    }

    return coarseAlg;
}

void showVec_(std::vector<int> const &IDs)
{
    for (auto const &value : IDs)
    {
        std::cout << value << "\n";
    }
}
template<bool withTemporal>
void showTemporalStatus_()
{
    if constexpr (withTemporal)
    {
        std::cout << "Perform temporal transaction\n";
    }
    else if constexpr (!withTemporal)
    {
        std::cout << "Perform non temporal transaction\n";
    }
}
template<bool withTemporal, typename Enum, Enum fillType>
void showAction_()
{
    std::cout << "***************************************\n";
    showTemporalStatus_<withTemporal>();

    if constexpr (fillType == Enum::SameLevel)
    {
        std::cout << "Perform same level transaction\n";
    }
    if constexpr (fillType == Enum::GhostRegionOnSameLevel)
    {
        std::cout << "Perform ghost region on same level transaction\n";
    }
    if constexpr (fillType == Enum::LevelBorderOnly)
    {
        std::cout << "Perform coarse to fine level boundary  only transaction\n";
    }
    if constexpr (fillType == Enum::GhostRegion)
    {
        std::cout << "Perform ghost region transaction\n";
    }
    if constexpr (fillType == Enum::EraseDestination)
    {
        std::cout << "Perform ghost and interior transaction (interpolation + initial split) "
                     "transaction\n";
    }

    std::cout << "***************************************\n";
}

HybridState const &model_;

std::remove_reference_t<decltype(model_.electromag.E)> electricNew_;
std::remove_reference_t<decltype(model_.electromag.B)> magneticNew_;

std::map<std::string, std::shared_ptr<Algorithm>> algorithmsMagneticIn_;
std::map<std::string, std::shared_ptr<Algorithm>> algorithmsMagneticOut_;

std::map<std::string, std::shared_ptr<Algorithm>> algorithmsElectricIn_;
std::map<std::string, std::shared_ptr<Algorithm>> algorithmsElectricOut_;

std::map<std::string, std::shared_ptr<Algorithm>> algorithmsIonsIn_;
std::map<std::string, std::shared_ptr<Algorithm>> algorithmsIonsOut_;

std::map<std::string, std::shared_ptr<Algorithm>> algorithmsMagneticSync_;
std::map<std::string, std::shared_ptr<Algorithm>> algorithmsElectricSync_;

ResourcesManager &resources_;

std::shared_ptr<SAMRAI::hier::PatchHierarchy> const *hierarchy_;
std::shared_ptr<SAMRAI::hier::PatchLevel> const *currentLevel_;
std::shared_ptr<SAMRAI::hier::PatchLevel> const *oldLevel_;

std::shared_ptr<SAMRAI::hier::RefineOperator> const &fieldRefineOp_;
std::shared_ptr<SAMRAI::hier::CoarsenOperator> const &fieldCoarsenOp_;

std::shared_ptr<SAMRAI::hier::RefineOperator> const &particlesRefineOp_;

int const startLevel_;

int relativeLevel_;
int currentLevelNumber_;

SAMRAI::tbox::Dimension dimension_;

bool alreadyRegisterElectricNew_{false};
bool alreadyRegisterMagneticNew_{false};
#endif
};

template<typename HybridModel>
const std::string HybridHybridTransactionStrategy<HybridModel>::stratName
    = "HybridModel-HybridModel";

/*
template<typename HybridState, typename ResourcesManager>
class HybridTransactionStrategyFactory
{
private:
using VecFieldT   = decltype(std::declval<HybridState>().electromag.E);
using IonsT       = decltype(std::declval<HybridState>().ions);
using HybridStrat = HybridTransactionStrategy<HybridState>;
using HybHybStrat = HybridHybridTransactionStrategy<HybridState, ResourcesManager>;
using MHDHybStrat = MHDHybridTransactionStrategy<HybridState>;
// using mhdmhdStrat = MHDMHDTransactionStrategy;


public:
static std::unique_ptr<HybridStrat> createStrategy(HybridTransactionStrategyType type)
{
    switch (type)
    {
        case HybridTransactionStrategyType::HybridHybrid:
            return std::make_unique<HybHybStrat>();
            break;

        case HybridTransactionStrategyType::MHDHybrid:
            return std::make_unique<MHDHybStrat>();
            break;
    }
}
};

*/



} // namespace PHARE

#endif
