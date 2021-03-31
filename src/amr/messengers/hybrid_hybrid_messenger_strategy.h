
#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_H
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_H

#include "communicators.h"
#include "amr/data/field/coarsening/field_coarsen_operator.h"
#include "amr/data/field/refine/field_refine_operator.h"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.h"
#include "amr/data/particles/refine/particles_data_split.h"
#include "amr/data/particles/refine/split.h"
#include "amr/messengers/hybrid_messenger_info.h"
#include "amr/messengers/hybrid_messenger_strategy.h"
#include "amr/resources_manager/amr_utils.h"
#include "amr/resources_manager/resources_manager_utilities.h"

#include "core/numerics/interpolator/interpolator.h"
#include "core/numerics/moments/moments.h"
#include "core/hybrid/hybrid_quantities.h"



#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>


#include <iterator>
#include <optional>
#include <utility>
#include <iomanip>


namespace PHARE
{
namespace amr
{
    /** \brief An HybridMessenger is the specialization of a HybridMessengerStrategy for hybrid to
     * hybrid data communications.
     */
    template<typename HybridModel, typename RefinementParams>
    class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
    {
        using IonsT                              = typename HybridModel::ions_type;
        using ElectromagT                        = typename HybridModel::electromag_type;
        using VecFieldT                          = typename HybridModel::vecfield_type;
        using GridLayoutT                        = typename HybridModel::gridlayout_type;
        using FieldT                             = typename VecFieldT::field_type;
        using ResourcesManagerT                  = typename HybridModel::resources_manager_type;
        static constexpr std::size_t dimension   = GridLayoutT::dimension;
        static constexpr std::size_t interpOrder = GridLayoutT::interp_order;
        using IPhysicalModel                     = typename HybridModel::Interface;

        using InteriorParticleRefineOp = typename RefinementParams::InteriorParticleRefineOp;
        using CoarseToFineRefineOpOld  = typename RefinementParams::CoarseToFineRefineOpOld;
        using CoarseToFineRefineOpNew  = typename RefinementParams::CoarseToFineRefineOpNew;

        template<typename HybridModel_1, typename IPhysicalModel_1>
        friend class StaticHybridHybridMessengerStrategy;


    public:
        static const std::string stratName;
        static constexpr std::size_t rootLevelNumber = 0;


        HybridHybridMessengerStrategy(std::shared_ptr<ResourcesManagerT> manager,
                                      int const firstLevel)
            : HybridMessengerStrategy<HybridModel>{stratName}
            , resourcesManager_{std::move(manager)}
            , firstLevel_{firstLevel}
        {
            resourcesManager_->registerResources(EM_old_);
            resourcesManager_->registerResources(Jold_);
        }

        virtual ~HybridHybridMessengerStrategy() = default;



        /* ------------------------------------------------------------------------
                        methods used for the IMessenger interface
           ------------------------------------------------------------------------ */


        /**
         * @brief allocate the messenger strategy internal variables to the model resourceManager
         */
        void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            resourcesManager_->allocate(EM_old_, patch, allocateTime);
            resourcesManager_->allocate(Jold_, patch, allocateTime);
        }



        /**
         * @brief setup creates all SAMRAI algorithms to communicate data involved in a messenger
         * between the coarse and fine levels.
         *
         * This method creates the SAMRAI algorithms for communications associated between pairs of
         * variables. The function does not create the SAMRAI schedules since they depend on the
         * levels
         */
        void registerQuantities([[maybe_unused]] std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<HybridMessengerInfo> hybridInfo{
                dynamic_cast<HybridMessengerInfo*>(fromFinerInfo.release())};

            registerGhostComms_(hybridInfo);
            registerInitComms(hybridInfo);
            registerSyncComms(hybridInfo);
        }



        /**
         * @brief registerLevel registers the level for all Communicators
         *
         * The level must always be registered to ghost Communicators but only to init Communicators
         * if it is not the root level since root is not initialized by coarser data.
         *
         * Ghost communicators that need to know this level :
         *
         *  - magnetic fields
         *  - electric fields
         *  - patch ghost particles
         *
         *  ion moments : do not need to be filled on ghost node by SAMRAI schedules
         *  since they will be filled with levelGhostParticles[old,new] on level ghost nodes
         *  and computed by ghost particles on interior patch ghost nodes
         *
         *
         * Init communicators that need to know this level :
         *
         *  - magnetic fields
         *  - electric fields
         *  - ion interior particle arrays
         *  - ion levelGhostParticlesOld particle arrays
         */
        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magneticGhosts_.registerLevel(hierarchy, level);
            electricGhosts_.registerLevel(hierarchy, level);
            currentGhosts_.registerLevel(hierarchy, level);
            patchGhostParticles_.registerLevel(hierarchy, level);
            densityGhosts_.registerLevel(hierarchy, level);

            // root level is not initialized with a schedule using coarser level data
            // so we don't create these schedules if root level
            if (levelNumber != rootLevelNumber)
            {
                // those are for refinement
                magneticInit_.registerLevel(hierarchy, level);
                electricInit_.registerLevel(hierarchy, level);
                interiorParticles_.registerLevel(hierarchy, level);
                levelGhostParticlesOld_.registerLevel(hierarchy, level);
                levelGhostParticlesNew_.registerLevel(hierarchy, level);

                // and these for coarsening
                magnetoSynchronizers_.registerLevel(hierarchy, level);
                electroSynchronizers_.registerLevel(hierarchy, level);
                densitySynchronizers_.registerLevel(hierarchy, level);
                ionBulkVelSynchronizers_.registerLevel(hierarchy, level);
            }
        }



        /**
         * @brief regrid performs the regriding communications for Hybrid to Hybrid messengers
         *
         * basically, all quantities that are in initialization refiners need to be regridded
         */
        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                    const int levelNumber,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                    IPhysicalModel& model, double const initDataTime) override
        {
            auto level = hierarchy->getPatchLevel(levelNumber);
            magneticInit_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            electricInit_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            interiorParticles_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            levelGhostParticlesOld_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            copyLevelGhostOldToPushable_(*level, model);

            // computeIonMoments_(*level, model);
            // levelGhostNew will be refined in next firstStep
        }




        std::string fineModelName() const override { return HybridModel::model_name; }



        std::string coarseModelName() const override { return HybridModel::model_name; }




        std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<HybridMessengerInfo>();
        }




        std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<HybridMessengerInfo>();
        }



        /**
         * @brief initLevel is used to initialize hybrid data on the level levelNumer at time
         * initDataTime from hybrid coarser data.
         */
        void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                       double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();

            magneticInit_.fill(levelNumber, initDataTime);
            electricInit_.fill(levelNumber, initDataTime);

            // no need to call these :
            // magneticGhosts_.fill(levelNumber, initDataTime);
            // electricGhosts_.fill(levelNumber, initDataTime);
            // because the SAMRAI schedules in the 'init' communicators
            // already fill the patch ghost box from the neighbor interior box.
            // so ghost nodes are already filled .

            interiorParticles_.fill(levelNumber, initDataTime);
            // however we need to call the ghost communicator for patch ghost particles
            // since the interior schedules have a restriction to the interior of the patch.
            patchGhostParticles_.fill(levelNumber, initDataTime);


            levelGhostParticlesOld_.fill(levelNumber, initDataTime);


            // levelGhostParticles will be pushed during the advance phase
            // they need to be identical to levelGhostParticlesOld before advance
            copyLevelGhostOldToPushable_(level, model);

            // computeIonMoments_(level, model);
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
        void fillMagneticGhosts(VecFieldT& B, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillMagneticGhosts");
            magneticGhosts_.fill(B, levelNumber, fillTime);
        }




        void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillElectricGhosts");
            electricGhosts_.fill(E, levelNumber, fillTime);
        }




        void fillCurrentGhosts(VecFieldT& J, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillCurrentGhosts");
            currentGhosts_.fill(J, levelNumber, fillTime);
        }


        void fillDensityGhosts(int const levelNumber, double const fillTime) override
        {
            densityGhosts_.fill(levelNumber, fillTime);
        }



        /**
         * @brief fillIonGhostParticles will fill the interior ghost particle array from neighbor
         * patches of the same level. Before doing that, it empties the array for all populations
         */
        void fillIonGhostParticles(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                   double const fillTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillIonGhostParticles");

            for (auto patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                for (auto& pop : ions)
                {
                    empty(pop.patchGhostParticles());
                }
            }

            patchGhostParticles_.fill(level.getLevelNumber(), fillTime);
        }




        /**
         * @brief fillIonMomentGhosts will compute the ion moments for all species on ghost nodes.
         *
         * For patch ghost nodes, patch ghost particles are used.
         * For level ghost nodes, levelGhostParticlesOld and new are both projected with a time
         * interpolation coef.
         */
        void fillIonMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                 double const beforePushTime, double const afterPushTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillIonMomentGhosts");

            auto alpha = timeInterpCoef_(beforePushTime, afterPushTime, level.getLevelNumber());
            if (level.getLevelNumber() > 0 and (alpha < 0 or alpha > 1))
            {
                std::cout << std::setprecision(12) << alpha << "\n";
                throw std::runtime_error("ion moment ghost time interp coef invalid : alpha: "
                                         + std::to_string(alpha) + " beforePushTime "
                                         + std::to_string(beforePushTime) + " afterPushTime "
                                         + std::to_string(afterPushTime) + " on level "
                                         + std::to_string(level.getLevelNumber()));
            }
            for (auto patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                auto layout      = layoutFromPatch<GridLayoutT>(*patch);

                for (auto& pop : ions)
                {
                    // first thing to do is to project patchGhostParitcles moments
                    auto& patchGhosts = pop.patchGhostParticles();
                    auto& density     = pop.density();
                    auto& flux        = pop.flux();

                    interpolate_(std::begin(patchGhosts), std::end(patchGhosts), density, flux,
                                 layout);

                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // then grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        auto& levelGhostOld = pop.levelGhostParticlesOld();
                        interpolate_(std::begin(levelGhostOld), std::end(levelGhostOld), density,
                                     flux, layout, 1. - alpha);

                        auto& levelGhostNew = pop.levelGhostParticlesNew();
                        interpolate_(std::begin(levelGhostNew), std::end(levelGhostNew), density,
                                     flux, layout, alpha);
                    }
                }
            }
        }



        /**
         * @brief firstStep : in the HybridHybridMessengerStrategy, the firstStep method is used to
         * get level border ghost particles from the next coarser level. These particles are defined
         * in the future at the time the method is called because the coarser level is ahead in
         * time. These particles are communicated only at first step of a substepping cycle. They
         * will be used with the levelGhostParticlesOld particles to get the moments on level border
         * nodes.
         * The method is does nothing if the level is the root level because the root level
         * cannot get levelGhost from next coarser (it has none).
         */
        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& level,
                       std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                       double const currentTime, double const prevCoarserTime,
                       double const newCoarserTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::firstStep");

            auto levelNumber = level.getLevelNumber();
            if (newCoarserTime < prevCoarserTime)
                throw std::runtime_error(
                    "Error : prevCoarserTime (" + std::to_string(prevCoarserTime)
                    + ") should be < newCoarserTime (" + std::to_string(prevCoarserTime) + ")");

            // root level has no levelghost particles
            if (levelNumber != 0)
            {
                PHARE_LOG_START("HybridHybridMessengerStrategy::firstStep.fill");
                levelGhostParticlesNew_.fill(levelNumber, currentTime);
                PHARE_LOG_STOP("HybridHybridMessengerStrategy::firstStep.fill");

                // during firstStep() coarser level and current level are at the same time
                // so 'time' is also the beforePushCoarseTime_
                beforePushCoarseTime_[levelNumber] = prevCoarserTime;
                afterPushCoarseTime_[levelNumber]  = newCoarserTime;
            }
        }


        /**
         * @brief lastStep is used to perform operations at the last step of a substepping cycle.
         * It is called after the level is advanced. Here for hybrid-hybrid messages, the method
         * moves levelGhostParticlesNew particles into levelGhostParticlesOld ones. Then
         * levelGhostParticlesNew are emptied since it will be filled again at firstStep of the next
         * substepping cycle. the new CoarseToFineOld content is then copied to levelGhostParticles
         * so that they can be pushed during the next subcycle
         */
        void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) override
        {
            if (level.getLevelNumber() > 0)
            {
                PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::lastStep");

                auto& hybridModel = static_cast<HybridModel&>(model);
                for (auto& patch : level)
                {
                    auto& ions       = hybridModel.state.ions;
                    auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                    for (auto& pop : ions)
                    {
                        auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                        auto& levelGhostParticlesNew = pop.levelGhostParticlesNew();
                        auto& levelGhostParticles    = pop.levelGhostParticles();

                        core::swap(levelGhostParticlesNew, levelGhostParticlesOld);
                        core::empty(levelGhostParticlesNew);
                        core::empty(levelGhostParticles);
                        std::copy(std::begin(levelGhostParticlesOld),
                                  std::end(levelGhostParticlesOld),
                                  std::back_inserter(levelGhostParticles));

                        if (level.getLevelNumber() == 0)
                        {
                            if (levelGhostParticlesNew.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesNew detected in root level : "
                                    + std::to_string(levelGhostParticlesNew.size()));
                            if (levelGhostParticles.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticles detected in root level : "
                                    + std::to_string(levelGhostParticles.size()));
                            if (levelGhostParticlesOld.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesOld detected in root level : "
                                    + std::to_string(levelGhostParticlesOld.size()));
                        }
                    }
                }
            }
        }



        /**
         * @brief prepareStep is the concrete implementation of the
         * HybridMessengerStrategy::prepareStep method For hybrid-Hybrid communications.
         * This method copies the current model electromagnetic field and current, defined at
         * t=n. Since prepareStep() is called just before advancing the level, this operation
         * actually saves the t=n electromagnetic field and current into the messenger. When the
         * time comes that the next finer level needs to time interpolate the electromagnetic
         * field and current at its ghost nodes, this level will have its model EM field  and
         * current at t=n+1 and thanks to this methods, the t=n field will be in the messenger.
         */
        void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                         double currentTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::prepareStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, hybridModel.state.electromag, hybridModel.state.J, EM_old_, Jold_);
                resourcesManager_->setTime(EM_old_, *patch, currentTime);
                resourcesManager_->setTime(Jold_, *patch, currentTime);

                auto& EM = hybridModel.state.electromag;
                auto& J  = hybridModel.state.J;
                EM_old_.copyData(EM);
                Jold_.copyData(J);
            }
        }




        void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                            double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();
            assert(levelNumber == 0);

            auto& hybridModel = static_cast<HybridModel&>(model);

            magneticGhosts_.fill(hybridModel.state.electromag.B, levelNumber, initDataTime);
            electricGhosts_.fill(hybridModel.state.electromag.E, levelNumber, initDataTime);
            patchGhostParticles_.fill(levelNumber, initDataTime);

            // at some point in the future levelGhostParticles could be filled with injected
            // particles depending on the domain boundary condition.
        }



        void synchronize(SAMRAI::hier::PatchLevel& level) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::synchronize");

            auto levelNumber = level.getLevelNumber();

            // call coarsning schedules...
            magnetoSynchronizers_.sync(levelNumber);
            electroSynchronizers_.sync(levelNumber);
            // densitySynchronizers_.sync(levelNumber);
            // ionBulkVelSynchronizers_.sync(levelNumber);
        }

        void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                             double const time) override
        {
            auto levelNumber  = level.getLevelNumber();
            auto& hybridModel = static_cast<HybridModel&>(model);
            magneticGhosts_.fill(hybridModel.state.electromag.B, levelNumber, time);
            electricGhosts_.fill(hybridModel.state.electromag.E, levelNumber, time);
        }

    private:
        void registerGhostComms_(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto const& Eold = EM_old_.E;
            auto const& Bold = EM_old_.B;


            fillRefiners_(info->ghostElectric, info->modelElectric, VecFieldDescriptor{Eold},
                          electricGhosts_);

            fillRefiners_(info->ghostMagnetic, info->modelMagnetic, VecFieldDescriptor{Bold},
                          magneticGhosts_);

            fillRefiners_(info->ghostCurrent, info->modelCurrent, VecFieldDescriptor{Jold_},
                          currentGhosts_);

            fillRefiners_(info->ghostIonDensity, densityGhosts_);
        }




        void registerInitComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& descriptor) {
                std::vector<std::string> keys;
                std::transform(std::begin(descriptor), std::end(descriptor),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            fillRefiners_(info->initMagnetic, fieldRefineOp_, magneticInit_,
                          makeKeys(info->initMagnetic));

            fillRefiners_(info->initElectric, fieldRefineOp_, electricInit_,
                          makeKeys(info->initElectric));


            fillRefiners_(info->interiorParticles, interiorParticleRefineOp_, interiorParticles_,
                          info->interiorParticles);


            fillRefiners_(info->levelGhostParticlesOld, levelGhostParticlesOldOp_,
                          levelGhostParticlesOld_, info->levelGhostParticlesOld);


            fillRefiners_(info->levelGhostParticlesNew, levelGhostParticlesNewOp_,
                          levelGhostParticlesNew_, info->levelGhostParticlesNew);


            fillRefiners_(info->patchGhostParticles, nullptr, patchGhostParticles_,
                          info->patchGhostParticles);
        }




        void registerSyncComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            magnetoSynchronizers_.add(info->modelMagnetic, resourcesManager_, fieldCoarseningOp_,
                                      info->modelMagnetic.vecName);

            electroSynchronizers_.add(info->modelElectric, resourcesManager_, fieldCoarseningOp_,
                                      info->modelElectric.vecName);

            ionBulkVelSynchronizers_.add(info->modelIonBulkVelocity, resourcesManager_,
                                         fieldCoarseningOp_, info->modelIonBulkVelocity.vecName);


            densitySynchronizers_.add(info->modelIonDensity, fieldCoarseningOp_,
                                      info->modelIonDensity, resourcesManager_);
        }


        /**
         * @brief makeCommunicators_ adds to the ghost communicators all VecFieldDescriptor of
         * the given vector field.
         *
         * Each of the ghost VecFieldDescriptor will have an entry in the ghost communicators
         *
         * @param ghostVec is the collection of VecFieldDescriptor. Each VecFieldDescriptor
         * corresponds to a VecField for which ghosts will be needed.
         * @param modelVec is VecFieldDescriptor for the model VecField associated with the
         * VecField for which ghosts are needed. When ghosts are filled, this quantity is taken
         * on the coarser level and is definer at t_coarse+dt_coarse
         * @param oldModelVec is the VecFieldDescriptor for the VecField for which ghosts are
         * needed, at t_coarse. These are typically internal variables of the messenger, like
         * Eold or Bold.
         */
        void fillRefiners_(std::vector<VecFieldDescriptor> const& ghostVecs,
                           VecFieldDescriptor const& modelVec,
                           VecFieldDescriptor const& oldModelVec,
                           RefinerPool<RefinerType::GhostField>& refiners)
        {
            for (auto const& ghostVec : ghostVecs)
            {
                refiners.add(ghostVec, modelVec, oldModelVec, resourcesManager_, fieldRefineOp_,
                             fieldTimeOp_, ghostVec.vecName);
            }
        }




        template<typename Descriptors, typename RefinerPool>
        void fillRefiners_(Descriptors const& descriptors,
                           std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                           RefinerPool& refiners, std::vector<std::string> keys)
        {
            auto key = std::begin(keys);
            for (auto const& descriptor : descriptors)
            {
                refiners.add(descriptor, refineOp, *key++, resourcesManager_);
            }
        }

        template<typename RefinerPool>
        void fillRefiners_(FieldDescriptor const& descriptor, RefinerPool& refiner)
        {
            refiner.add(descriptor, nullptr, descriptor, resourcesManager_);
        }


        void copyLevelGhostOldToPushable_(SAMRAI::hier::PatchLevel& level, IPhysicalModel& model)
        {
            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level)
            {
                auto& ions       = hybridModel.state.ions;
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                for (auto& pop : ions)
                {
                    auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                    auto& levelGhostParticles    = pop.levelGhostParticles();

                    core::empty(levelGhostParticles);
                    std::copy(std::begin(levelGhostParticlesOld), std::end(levelGhostParticlesOld),
                              std::back_inserter(levelGhostParticles));
                }
            }
        }




        double timeInterpCoef_(double const beforePushTime, double const afterPushTime,
                               std::size_t levelNumber)
        {
            return (afterPushTime - beforePushTime)
                   / (afterPushCoarseTime_[levelNumber] - beforePushCoarseTime_[levelNumber]);
        }




        //! keeps a copy of the model electromagnetic field at t=n
        ElectromagT EM_old_{stratName + "_EM_old"}; // TODO needs to be allocated somewhere and
                                                    // updated to t=n before advanceLevel()

        VecFieldT Jold_{stratName + "_Jold", core::HybridQuantity::Vector::J};


        //! ResourceManager shared with other objects (like the HybridModel)
        std::shared_ptr<ResourcesManagerT> resourcesManager_;


        int const firstLevel_;
        std::unordered_map<std::size_t, double> beforePushCoarseTime_;
        std::unordered_map<std::size_t, double> afterPushCoarseTime_;

        core::Interpolator<dimension, interpOrder> interpolate_;


        //! store communicators for magnetic fields that need ghosts to be filled
        RefinerPool<RefinerType::GhostField> magneticGhosts_;

        //! store communicators for magnetic fields that need to be initialized
        RefinerPool<RefinerType::InitField> magneticInit_;

        //! store refiners for electric fields that need ghosts to be filled
        RefinerPool<RefinerType::GhostField> electricGhosts_;

        //! store communicators for electric fields that need to be initializes
        RefinerPool<RefinerType::InitField> electricInit_;


        RefinerPool<RefinerType::GhostField> currentGhosts_;

        RefinerPool<RefinerType::GhostField> densityGhosts_;

        // algo and schedule used to initialize domain particles
        // from coarser level using particleRefineOp<domain>
        RefinerPool<RefinerType::InitInteriorPart> interiorParticles_;

        //! store communicators for coarse to fine particles old
        RefinerPool<RefinerType::LevelBorderParticles> levelGhostParticlesOld_;

        //! store communicators for coarse to fine particles new
        RefinerPool<RefinerType::LevelBorderParticles> levelGhostParticlesNew_;

        // keys : model particles (initialization and 2nd push), temporaryParticles (firstPush)
        RefinerPool<RefinerType::InteriorGhostParticles> patchGhostParticles_;

        SynchronizerPool<dimension> densitySynchronizers_;

        SynchronizerPool<dimension> ionBulkVelSynchronizers_;

        SynchronizerPool<dimension> electroSynchronizers_;

        SynchronizerPool<dimension> magnetoSynchronizers_;


        std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_{
            std::make_shared<FieldRefineOperator<GridLayoutT, FieldT>>()};

        // field data time op
        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> fieldTimeOp_{
            std::make_shared<FieldLinearTimeInterpolate<GridLayoutT, FieldT>>()};


        std::shared_ptr<SAMRAI::hier::RefineOperator> interiorParticleRefineOp_{
            std::make_shared<InteriorParticleRefineOp>()};

        std::shared_ptr<SAMRAI::hier::RefineOperator> levelGhostParticlesOldOp_{
            std::make_shared<CoarseToFineRefineOpOld>()};

        std::shared_ptr<SAMRAI::hier::RefineOperator> levelGhostParticlesNewOp_{
            std::make_shared<CoarseToFineRefineOpNew>()};


        std::shared_ptr<SAMRAI::hier::CoarsenOperator> fieldCoarseningOp_{
            std::make_shared<FieldCoarsenOperator<GridLayoutT, FieldT>>()};
    };

    template<typename HybridModel, typename RefinementParams>
    const std::string HybridHybridMessengerStrategy<HybridModel, RefinementParams>::stratName
        = "HybridModel-HybridModel";

} // namespace amr

} // namespace PHARE

#endif
