
#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_H
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_H

#include "data/field/coarsening/field_coarsen_operator.h"
#include "data/field/refine/field_refine_operator.h"
#include "data/field/time_interpolate/field_linear_time_interpolate.h"
#include "data/particles/refine/particles_data_split.h"
#include "data/particles/refine/split.h"
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
/** \brief An HybridMessenger is the specialization of a HybridMessengerStrategy for hybrid to
 * hybrid data communications.
 */
template<typename HybridModel>
class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
{
    using IonsT                                = typename HybridModel::ions_type;
    using ElectromagT                          = typename HybridModel::electromag_type;
    using VecFieldT                            = typename HybridModel::vecfield_type;
    using GridLayoutT                          = typename HybridModel::gridLayout_type;
    using FieldT                               = typename VecFieldT::field_type;
    using ResourcesManagerT                    = typename HybridModel::resources_manager_type;
    static constexpr std::size_t dimension     = GridLayoutT::dimension;
    static constexpr std::size_t interpOrder   = GridLayoutT::interp_order;
    using SplitT                               = Split<dimension, interpOrder>;
    static constexpr std::size_t nbRefinedPart = 2; // TODO stop hard-coding this
    using InteriorParticleRefineOp
        = ParticlesRefineOperator<dimension, interpOrder, ParticlesDataSplitType::interior,
                                  nbRefinedPart, SplitT>;

    using CoarseToFineRefineOp
        = ParticlesRefineOperator<dimension, interpOrder, ParticlesDataSplitType::coarseBoundaryOld,
                                  nbRefinedPart, SplitT>;




public:
    static const std::string stratName;
    static constexpr std::size_t rootLevelNumber = 0;


    HybridHybridMessengerStrategy(std::shared_ptr<ResourcesManagerT> manager, int const firstLevel)
        : HybridMessengerStrategy<HybridModel>{stratName}
        , resourcesManager_{std::move(manager)}
        , firstLevel_{firstLevel}
    {
        resourcesManager_->registerResources(EM_old_);
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

        registerGhosts_(hybridInfo);
        registerInit_(hybridInfo);
    }



    /**
     * @brief setLevel creates SAMRAI schedules for all resources in the ghost and init RefinerPools
     *
     * Needing ghosts to be filled:
     *  - magnetic fields
     *  - electric fields
     *  - particles
     *
     *  ion moments do not need to be filled on ghost node by SAMRAI schedules
     *  since they will be filled with coarseToFine particles on level ghost nodes
     *  and computed by ghost particles on interior patch ghost nodes
     *
     * Needing to be initialized
     *  - magnetic fields
     *  - electric fields
     *  - ion bulk velocity (total)
     *  - ion density (total)
     *  - ion interior particle arrays
     *  - ion coarseToFineOld particle arrays
     */
    virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                               int const levelNumber) override
    {
        auto level = hierarchy->getPatchLevel(levelNumber);

        magneticGhostsRefiners_.createGhostSchedules(hierarchy, level);
        electricGhostsRefiners_.createGhostSchedules(hierarchy, level);
        ghostParticleRefiners_.createGhostSchedules(hierarchy, level);

        // root level is not initialized with a schedule using coarser level data
        // so we don't create these schedules if root level
        if (levelNumber != rootLevelNumber)
        {
            magneticInitRefiners_.createInitSchedules(hierarchy, level);
            electricInitRefiners_.createInitSchedules(hierarchy, level);
            ionBulkInitRefiners_.createInitSchedules(hierarchy, level);
            ionDensityInitRefiners_.createInitSchedules(hierarchy, level);
            interiorParticleRefiners_.createInitSchedules(hierarchy, level);
            coarseToFineParticleRefiners_.createInitSchedules(hierarchy, level);
        }
    }



    /**
     * @brief regrid performs the regriding communications for Hybrid to Hybrid messengers
     *
     * basically, all quantities that are in initialization refiners need to be regridded
     */
    virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                        const int levelNumber,
                        std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                        double const initDataTime) override
    {
        magneticInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
        electricInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
        ionBulkInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
        ionDensityInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
        // TODO regrid particle arrays
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
        magneticInitRefiners_.initialize(levelNumber, initDataTime);
        electricInitRefiners_.initialize(levelNumber, initDataTime);
        ionBulkInitRefiners_.initialize(levelNumber, initDataTime);
        ionDensityInitRefiners_.initialize(levelNumber, initDataTime);
        interiorParticleRefiners_.initialize(levelNumber, initDataTime);
        coarseToFineParticleRefiners_.initialize(levelNumber, initDataTime);
        // ghostParticleRefiners_.initialize(levelNumber, initDataTime);
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
        magneticGhostsRefiners_.fillVecFieldGhosts(B, levelNumber, fillTime);
    }




    virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber,
                                    double const fillTime) override
    {
        std::cout << "perform the electric ghost fill\n";
        electricGhostsRefiners_.fillVecFieldGhosts(E, levelNumber, fillTime);
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


    virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) override
    {
        auto& hybridModel = static_cast<HybridModel&>(model);
        for (auto& patch : level)
        {
            auto dataOnPatch
                = resourcesManager_->setOnPatch(*patch, hybridModel.state.electromag, EM_old_);

            auto& EM = hybridModel.state.electromag;
            EM_old_.copyData(EM);
        }
    }


private:
    void registerGhosts_(std::unique_ptr<HybridMessengerInfo> const& info)
    {
        auto const& Eold = EM_old_.E;
        auto const& Bold = EM_old_.B;


        addToGhostRefinerPool_(info->ghostElectric, info->modelElectric, VecFieldDescriptor{Eold},
                               electricGhostsRefiners_);

        addToGhostRefinerPool_(info->ghostMagnetic, info->modelMagnetic, VecFieldDescriptor{Bold},
                               magneticGhostsRefiners_);


        // TODO add particle infos to particle ghost refiner pool
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




    void registerInit_(std::unique_ptr<HybridMessengerInfo> const& info)
    {
        auto makeKeys = [](auto const& descriptor) {
            std::vector<std::string> keys;
            std::transform(std::begin(descriptor), std::end(descriptor), std::back_inserter(keys),
                           [](auto const& d) { return d.vecName; });
            return keys;
        };

        addToInitRefinerPool_(info->initMagnetic, fieldRefineOp_, magneticInitRefiners_,
                              makeKeys(info->initMagnetic));

        addToInitRefinerPool_(info->initElectric, fieldRefineOp_, electricInitRefiners_,
                              makeKeys(info->initElectric));

        addToInitRefinerPool_(info->initIonBulk, fieldRefineOp_, ionBulkInitRefiners_,
                              makeKeys(info->initIonBulk));

        addToInitRefinerPool_(info->initIonDensity, fieldRefineOp_, ionDensityInitRefiners_,
                              info->initIonDensity);

        addToInitRefinerPool_(info->interiorParticles, interiorParticleRefineOp_,
                              interiorParticleRefiners_, info->interiorParticles);


        addToInitRefinerPool_(info->coarseToFineParticles, coarseToFineRefineOp_,
                              coarseToFineParticleRefiners_, info->coarseToFineParticles);


        addToInitRefinerPool_(info->ghostParticles, nullptr, ghostParticleRefiners_,
                              info->ghostParticles);
    }



    template<typename Descriptors>
    void addToInitRefinerPool_(Descriptors const& descriptors,
                               std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                               RefinerPool& refinerPool, std::vector<std::string> keys)
    {
        auto key = std::begin(keys);
        for (auto const& descriptor : descriptors)
        {
            refinerPool.add(makeInitRefiner(descriptor, resourcesManager_, refineOp), *key++);
        }
    }

    //! keeps a copy of the model electromagnetic field at t=n
    ElectromagT EM_old_{stratName + "_EM_old"}; // TODO needs to be allocated somewhere and
                                                // updated to t=n before advanceLevel()


    //! ResourceManager shared with other objects (like the HybridModel)
    std::shared_ptr<typename HybridModel::resources_manager_type> resourcesManager_;


    int const firstLevel_;


    //! store refiners for magnetic fields that need ghosts to be filled
    RefinerPool magneticGhostsRefiners_;

    //! store refiners for magnetic fields that need to be initialized
    RefinerPool magneticInitRefiners_;

    //! store refiners for electric fields that need ghosts to be filled
    RefinerPool electricGhostsRefiners_;

    //! store refiners for electric fields that need to be initializes
    RefinerPool electricInitRefiners_;

    //! store refiners for ion bulk velocity resources that need to be initialized
    RefinerPool ionBulkInitRefiners_;

    //! store refiners for total ion density resources that need to be initialized
    RefinerPool ionDensityInitRefiners_;

    // algo and schedule used to initialize domain particles
    // from coarser level using particleRefineOp<domain>
    RefinerPool interiorParticleRefiners_;

    //! store refiners for coarse to fine particles
    RefinerPool coarseToFineParticleRefiners_;

    // keys : model particles (initialization and 2nd push), temporaryParticles (firstPush)
    RefinerPool ghostParticleRefiners_;

    // at first step of advance:
    // from temporaryParticle of coarseLevel to model PRA1 ( + PRA1 copy into PRA)
    // and from modelParticles of coarserLevel to model PRA2
    // the copy of PRA1 vector to PRA is done after the schedule in a PatchStrategy post truc
    // these are ran before solver->advanceLevel() in the MultiPhysics::advanceLevel()
    // with a method : messenger.firstStepOperation() or somthg like that...
    // solver->advanceLevel() when starting, has all PRAs set correctly

    // keys: PRA1, PRA2 , chosen by messenger
    RefinerPool particlePRA_;



    std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_{
        std::make_shared<FieldRefineOperator<GridLayoutT, FieldT>>()};

    // field data time op
    std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> fieldTimeOp_{
        std::make_shared<FieldLinearTimeInterpolate<GridLayoutT, FieldT>>()};


    std::shared_ptr<SAMRAI::hier::RefineOperator> interiorParticleRefineOp_{
        std::make_shared<InteriorParticleRefineOp>()};

    std::shared_ptr<SAMRAI::hier::RefineOperator> coarseToFineRefineOp_{
        std::make_shared<CoarseToFineRefineOp>()};
};

template<typename HybridModel>
const std::string HybridHybridMessengerStrategy<HybridModel>::stratName = "HybridModel-HybridModel";



} // namespace PHARE

#endif
