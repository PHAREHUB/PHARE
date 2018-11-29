
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
        magneticGhostsRefiners_.createGhostSchedules(hierarchy, level);
        electricGhostsRefiners_.createGhostSchedules(hierarchy, level);

        // now update electric and magnetic initialization
        // root level is not initialized with a schedule using coarser level data
        // so we don't create these schedules if root level
        if (levelNumber != 0)
        {
            magneticInitRefiners_.createInitSchedules(hierarchy, level);
            electricInitRefiners_.createInitSchedules(hierarchy, level);
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


        auto schedule = magneticGhostsRefiners_.findSchedule(B.name(), levelNumber);
        if (schedule)
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


        auto schedule = electricGhostsRefiners_.findSchedule(E.name(), levelNumber);
        if (schedule)
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

    using Electromag = decltype(std::declval<HybridModel>().state.electromag);

    using GridLayoutT = typename HybridModel::gridLayout_type;

    using field_type =
        typename decltype(std::declval<HybridModel>().state.electromag)::vecfield_type::field_type;


    /**
     * @brief setupEandBGhostMessengers_ adds to the ghost refiner pool all electric and magnetic
     * VecFieldDescriptors that are present in the given info.
     *
     */
    void setupEandBGhostMessengers_(std::unique_ptr<HybridMessengerInfo> const& info)
    {
        auto const& Eold = EM_old_.E;
        auto const& Bold = EM_old_.B;

        addToGhostRefinerPool_(info->ghostElectric, info->modelElectric, VecFieldDescriptor{Eold},
                               electricGhostsRefiners_);
        addToGhostRefinerPool_(info->ghostMagnetic, info->modelMagnetic, VecFieldDescriptor{Bold},
                               magneticGhostsRefiners_);
    }



    /**
     * @brief addToGhostRefinerPool_ adds to the ghost refiner pool all VecFieldDescriptor of the
     * given vector field.
     *
     * Each of the ghost VecFieldDescriptor will have an entry in the ghost refiner pool
     *
     * @param ghostVec is the collection of VecFieldDescriptor. Each VecFieldDescriptor corresponds
     * to a VecField for which ghosts will be needed.
     * @param modelVec is VecFieldDescriptor for the model VecField associated with the VecField for
     * which ghosts are needed. When ghosts are filled, this quantity is taken on the coarser level
     * and is definer at t_coarse+dt_coarse
     * @param oldModelVec is the VecFieldDescriptor for the VecField for which ghosts are needed, at
     * t_coarse. These are typically internal variables of the messenger, like Eold or Bold.
     * @param refinerPool is the RefinerPool to which we add the refiner to.
     */
    void addToGhostRefinerPool_(std::vector<VecFieldDescriptor> const& ghostVec,
                                VecFieldDescriptor const& modelVec,
                                VecFieldDescriptor const& oldModelVec, RefinerPool& refinerPool)
    {
        auto nbrVectors = ghostVec.size();
        for (auto i = 0u; i < nbrVectors; ++i)
        {
            refinerPool.add(makeGhostRefiner(ghostVec[i], modelVec, oldModelVec, resourcesManager_,
                                             fieldRefineOp_, fieldTimeOp_),
                            ghostVec[i].vecName);
        }
    }




    void setupEandBInitMessengers_(std::unique_ptr<HybridMessengerInfo> const& info)
    {
        addToInitRefinerPool_(info->initMagnetic, magneticInitRefiners_);
        addToInitRefinerPool_(info->initElectric, electricInitRefiners_);
    }




    void addToInitRefinerPool_(std::vector<VecFieldDescriptor> const& descriptors,
                               RefinerPool& refinerPool)
    {
        for (auto const& descriptor : descriptors)
        {
            refinerPool.add(makeInitRefiner(descriptor, resourcesManager_, fieldRefineOp_),
                            descriptor.vecName);
        }
    }


    /*

        void createGhostSchedules_(RefinerPool& refiners,
                                   std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   std::shared_ptr<SAMRAI::hier::PatchLevel>& level)
        {
            // clang-format off
            for (auto& [key, refiner] : refiners.qtyRefiners)
            // clang-format on
            {
                auto& algo    = refiner.algo;
                auto schedule = algo->createSchedule(level,
       level->getNextCoarserHierarchyLevelNumber(), hierarchy); refiner.add(schedule,
       level->getLevelNumber());
            }
        }
    */


    /*
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

                // note that here we must take that createsSchedule() overload and put nullptr as
       src
                // since we want to take from coarser level everywhere.
                // using the createSchedule overload that takes level, next_coarser_level only
                // would result in interior ghost nodes to be filled with interior of neighbor
       patches
                // but there is nothing there.
                refiner.add(algo->createSchedule(level, nullptr, levelNumber - 1, hierarchy),
                            levelNumber);
            }
        }
    */


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




    // field data refine op
    /*std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_{
        std::make_shared<FieldRefineOperator<gridlayout_type, field_type>>()};

    // field data time op
    std::shared_ptr<FieldLinearTimeInterpolate<gridlayout_type, field_type>> fieldTimeOp_{
        std::make_shared<FieldLinearTimeInterpolate<gridlayout_type, field_type>>()};
*/


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



    using FieldT = typename VecFieldT::field_type;


    std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_{
        std::make_shared<FieldRefineOperator<GridLayoutT, FieldT>>()};

    // field data time op
    std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> fieldTimeOp_{
        std::make_shared<FieldLinearTimeInterpolate<GridLayoutT, FieldT>>()};
};

template<typename HybridModel>
const std::string HybridHybridMessengerStrategy<HybridModel>::stratName = "HybridModel-HybridModel";



} // namespace PHARE

#endif
