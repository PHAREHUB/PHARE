
#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_H
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_H

#include "data/field/coarsening/field_coarsen_operator.h"
#include "data/field/refine/field_refine_operator.h"
#include "data/field/time_interpolate/field_linear_time_interpolate.h"
#include "data/particles/refine/particles_data_split.h"
#include "evolution/messengers/hybrid_messenger_info.h"
#include "evolution/messengers/hybrid_messenger_strategy.h"
#include "physical_models/physical_model.h"
#include "quantity_refiner.h"
#include "tools/resources_manager_utilities.h"

#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>


#include <optional>
#include <utility>


namespace PHARE
{
/** \brief An HybridMessenger purpose is to manage messenger from an to : the electric field,
 * magnetic field, and the ions.
 *
 *
 * TODO: we would also want to add rho soon
 * TODO: we have a lot of get* send* maybe we should use some enum instead,so that we can make
 * combinaison with other parameter easier
 *
 */
template<typename HybridModel>
class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
{
    using IonsT     = decltype(std::declval<HybridModel>().state.ions);
    using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

public:
    static const std::string stratName;


    HybridHybridMessengerStrategy(
        std::shared_ptr<typename HybridModel::resources_manager_type> manager, int const firstLevel)
        : HybridMessengerStrategy<HybridModel>{stratName}
        , resourcesManager_{std::move(manager)}
        , firstLevel_{firstLevel}
    {
        resourcesManager_->registerResources(EM_old_.E);
        resourcesManager_->registerResources(EM_old_.B);
    }

    virtual ~HybridHybridMessengerStrategy() = default;



    /* ------------------------------------------------------------------------
                    methods used for the IMessenger interface
       ------------------------------------------------------------------------ */


    /**
     * @brief allocate the messenger strategy internal variables to the model resourceManager
     */
    virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
    {
        // hybModel.resourcesManager->allocate(EM_old_.E, patch, allocateTime);
        // hybModel.resourcesManager->allocate(EM_old_.B, patch, allocateTime);
        // hybModel.resourcesManager->allocate(EM_old_, patch, allocateTime);
        resourcesManager_->allocate(EM_old_, patch, allocateTime);
    }



    /**
     * @brief setup creates all SAMRAI algorithms to communicate data involved in a messenger
     * between the coarse and fine levels.
     *
     * This method creates the SAMRAI algorithms for communications associated between pairs of
     * variables. The function does not create the SAMRAI schedules since they depend on the levels
     */
    virtual void
    registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                       [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo) override
    {
        std::unique_ptr<HybridMessengerInfo> hybridInfo{
            dynamic_cast<HybridMessengerInfo*>(fromCoarserInfo.release())};


        setupEandBGhostMessengers_(hybridInfo);
        setupEandBInitMessengers_(hybridInfo);

        // on a tous les algos pour remplir les ghosts de
        // - model E, B vers model E, B
        // - model E, B vers solver list of E_internal and list of B_internal

        // TODO setup particle messengers
    }



    /**
     * @brief setLevel creates SAMRAI schedules for all variables relevant for hybrid to hybrid
     * communications, whether they are from model, solver or messenger internals
     *
     * The method takes maps associating quantityName to the pair (algo, map of schedules)
     * for each algo, it creates a bunch of associated schedules and store them in the vector
     */
    virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               int const levelNumber) override
    {
        // update the schedules for ghost communications of the electric and magnetic field

        auto level = hierarchy->getPatchLevel(levelNumber);
        createGhostSchedules_(magneticGhostsRefiners_, hierarchy, level);
        createGhostSchedules_(electricGhostsRefiners_, hierarchy, level);

        // now update electric and magnetic initialization
        // root level is not initialized with a schedule using coarser level data
        // so we don't create these schedules if root level
        if (levelNumber != 0)
        {
            createInitSchedules_(magneticInitRefiners_, hierarchy, level);
            createInitSchedules_(electricInitRefiners_, hierarchy, level);
        }
    }



    /**
     * @brief regrid performs the regriding communications for Hybrid to Hybrid messengers
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
            for (auto& [key, refiner] : refiners.qtyRefiners)
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

        doRegrid(magneticInitRefiners_);
        doRegrid(electricInitRefiners_);
    }




    virtual std::string fineModelName() const override { return HybridModel::model_name; }




    virtual std::string coarseModelName() const override { return HybridModel::model_name; }




    virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
    {
        return std::make_unique<HybridMessengerInfo>();
    }




    virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
    {
        return std::make_unique<HybridMessengerInfo>();
    }




    virtual void initLevel(int const levelNumber, double const initDataTime) const override
    {
        applyInitSchedules_(levelNumber, initDataTime, magneticInitRefiners_);
        applyInitSchedules_(levelNumber, initDataTime, electricInitRefiners_);
    }



    /* ------------------------------------------------------------------------
                 methods used for the HybridMessenger interface
       ------------------------------------------------------------------------ */


    /**
     * @brief see IMessenger::fillMagneticGhosts for documentation

     * Note on the HybridHybrid version:
     * The function throws if the given magnetic field B has not been registered
     * in the ghostMagnetic field of the HybridMessengerInfo
     *
     * The method finds if the name of the VecField is
     *
     */
    virtual void fillMagneticGhosts(VecFieldT& B, int const levelNumber,
                                    double const fillTime) override
    {
        std::cout << "perform the magnetic ghost fill\n";


        auto schedule = magneticGhostsRefiners_.find(B.name(), levelNumber);
        if (*schedule)
        {
            (*schedule)->fillData(fillTime);
        }
        else
        {
            throw std::runtime_error("no schedule for the magnetic field " + B.name());
        }
    }



    virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber,
                                    double const fillTime) override
    {
        std::cout << "perform the electric ghost fill\n";

        auto schedule = electricGhostsRefiners_.find(E.name(), levelNumber);
        if (*schedule)
        {
            (*schedule)->fillData(fillTime);
        }
        else
        {
            throw std::runtime_error("no schedule for the electric field " + E.name());
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
    // using Refiners =




    /**
     * @brief setupEandBGhostMessengers_ creates the SAMRAI algorithms to transfer
     * electromagnetic fields on ghost nodes of level patches and register all refine operations.
     *
     * For each variable to register, we use the registerRefine overload that allows time
     * interpolation.
     *
     * The operation must use 'src_told' and 'src_tnew' versions of each field and
     * time/space interpolate that onto the ghost nodes of the same quantity onto the destination
     * level patches.
     *
     * The schedules originating from these algorithms will be executed in the solver
     * to fill the ghosts on a given level at time t_coarse + n*dt_fine.At this time, the ghosts are
     * obtained from time interpolation of the model fields at times t_coarse (src_told) and
     * t_coarse+dt_coarse (src_tnew).
     *
     * src_told IDs are obtained from the messenger Eold and Bold, which are copies of the model
     * fields at t_coarse (they are initialized with the level from the coarse level and swapped
     * with the solution at lastStep().)
     *
     * src_tnew IDs are obtained from the model IDs since on the coarser level, the model fields
     * already are at t_coarse + dt_coarse when the schedule is executed.
     *
     *
     * Each operation thus takes:
     *  - src       : variable IDs in model (not sure what it is used for here)
     *  - dest      : variable IDs in the MessengerInfo.ghostXXXXnames (could be from model, could
     * be from solver)
     *  - src_told  : variable IDs in Messenger
     *  - src_tnew  : variable IDs in model
     *
     *
     * The SAMRAI algorithms created by this function are stored in RefineAlgosAndSchedule
     * objects, themselves stored in an appropriate map with a key provided by the
     * HybridMessengerInfo
     *
     * Note: each time one of the 6 schedules created here will be executed, i.e. at each fine time
     * step, the fields at src_tnew will be communicated from the coarser level even though they
     * have not changed. This may be optimized later, for instance by storing the coarser model
     * fields at firstStep() and doing the time interpolation ourselves.
     *
     */
    void setupEandBGhostMessengers_(std::unique_ptr<HybridMessengerInfo> const& info)
    {
        auto const& Eold = EM_old_.E;
        auto const& Bold = EM_old_.B;


        // register a refine operation for all components in src/dest componentNames
        // using the components of the vecField_old as the old source
        auto registerRefine = [this](auto const& srcComponentNames, auto const& destComponentNames,
                                     auto const& names_old) //
        {
            QuantityRefiner refiner;
            refiner.algo = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();


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
                    refiner.algo->registerRefine(*dest_id, // dest
                                                 *src_id,  // source at same time
                                                 *old_id,  // source at past time (for time interp)
                                                 *src_id, // source at future time (for time interp)
                                                 *dest_id, // scratch
                                                 fieldRefineOp_, fieldTimeOp_);
                }
            }
            return refiner;
        };



        auto addToPool = [&registerRefine](VecFieldNames const& modelVec,
                                           std::vector<VecFieldNames> const& ghostVec,
                                           auto const& oldVecNames, auto& refinerPool) //
        {
            auto nbrVariables = ghostVec.size();
            for (auto i = 0u; i < nbrVariables; ++i)
            {
                std::array<std::string, 3> destComponentNames{
                    {ghostVec[i].xName, ghostVec[i].yName, ghostVec[i].zName}};

                std::array<std::string, 3> srcComponentNames{
                    {modelVec.xName, modelVec.yName, modelVec.zName}};

                refinerPool.add(registerRefine(srcComponentNames, destComponentNames, oldVecNames),
                                ghostVec[i].vecName);
            }
        };


        addToPool(info->modelMagnetic, info->ghostMagnetic, extractNames(Bold),
                  magneticGhostsRefiners_);
        addToPool(info->modelElectric, info->ghostElectric, extractNames(Eold),
                  electricGhostsRefiners_);
    }




    void setupEandBInitMessengers_(std::unique_ptr<HybridMessengerInfo> const& info)
    {
        auto registerRefine = [this](auto& srcComponentNames, auto& destComponentNames) //
        {
            QuantityRefiner refiner;
            refiner.algo = std::make_unique<SAMRAI::xfer::RefineAlgorithm>();

            for (auto componentIndex = 0u; componentIndex < 3; ++componentIndex)
            {
                auto src_id  = resourcesManager_->getID(srcComponentNames[componentIndex]);
                auto dest_id = resourcesManager_->getID(destComponentNames[componentIndex]);

                auto id1 = *src_id;
                auto id2 = *dest_id;
                if (src_id && dest_id)
                {
                    // dest, src, old, new, scratch
                    refiner.algo->registerRefine(*dest_id, // dest
                                                 *src_id,  // source at same time
                                                 *dest_id, // scratch
                                                 fieldRefineOp_);

                    return refiner;
                }
                else
                {
                    throw std::runtime_error("invalid IDs");
                }
            }
        };


        auto addToPool
            = [&registerRefine](std::vector<VecFieldNames> const& names, auto& refinerPool) //
        {
            auto nbrVariables = names.size();
            for (auto i = 0u; i < nbrVariables; ++i)
            {
                std::array<std::string, 3> componentNames{
                    {names[i].xName, names[i].yName, names[i].zName}};

                refinerPool.add(registerRefine(componentNames, componentNames), names[i].vecName);
            }
        };


        addToPool(info->initMagnetic, magneticInitRefiners_);
        addToPool(info->initElectric, electricInitRefiners_);
    }




    void createGhostSchedules_(RefinerPool& refiners,
                               std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
    {
        // clang-format off
        for (auto& [key, refiner] : refiners.qtyRefiners)
        // clang-format on
        {
            auto& algo    = refiner.algo;
            auto schedule = algo->createSchedule(level, level->getNextCoarserHierarchyLevelNumber(),
                                                 hierarchy);
            refiner.add(schedule, level->getLevelNumber());
        }
    }




    void createInitSchedules_(RefinerPool& refiners,
                              std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                              std::shared_ptr<SAMRAI::hier::PatchLevel> const& level)
    {
        // clang-format off
        for (auto& [key, refiner] : refiners.qtyRefiners)
        // clang-format on
        {
            auto& algo       = refiner.algo;
            auto levelNumber = level->getLevelNumber();

            // note that here we must take that createsSchedule() overload and put nullptr as src
            // since we want to take from coarser level everywhere.
            // using the createSchedule overload that takes level, next_coarser_level only
            // would result in interior ghost nodes to be filled with interior of neighbor patches
            // but there is nothing there.
            refiner.add(algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy),
                        levelNumber);
        }
    }



    void applyInitSchedules_(int levelNumber, double initDataTime,
                             RefinerPool const& refiners) const
    {
        // clang-format off
        for (auto& [key, refiner] : refiners.qtyRefiners)
        // clang-format on
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
    Electromag EM_old_{stratName + "_EM_old"}; // TODO needs to be allocated somewhere and
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
    RefinerPool magneticGhostsRefiners_;

    //! stores the algo and its associated schedules for getting magnetic field at initialization
    RefinerPool magneticInitRefiners_;

    // same as for the magnetic field
    RefinerPool electricGhostsRefiners_;
    RefinerPool electricInitRefiners_;


    // algo and schedule used to initialize domain particles
    // from coarser level using particleRefineOp<domain>
    QuantityRefiner particleInitRefine_;

    // keys : model particles (initialization and 2nd push), temporaryParticles (firstPush)
    RefinerPool particleGhostExchange_;

    // at first step of advance:
    // from temporaryParticle of coarseLevel to model PRA1 ( + PRA1 copy into PRA)
    // and from modelParticles of coarserLevel to model PRA2
    // the copy of PRA1 vector to PRA is done after the schedule in a PatchStrategy post truc
    // these are ran before solver->advanceLevel() in the MultiPhysics::advanceLevel()
    // with a method : messenger.firstStepOperation() or somthg like that...
    // solver->advanceLevel() when starting, has all PRAs set correctly

    // keys: PRA1, PRA2 , chosen by messenger
    RefinerPool particlePRA_;
};

template<typename HybridModel>
const std::string HybridHybridMessengerStrategy<HybridModel>::stratName = "HybridModel-HybridModel";



} // namespace PHARE

#endif
