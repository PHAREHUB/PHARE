#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP

#include "core/def.hpp" // IWYU pragma: keep
#include "core/logger.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "refiner_pool.hpp"
#include "synchronizer_pool.hpp"

#include "amr/data/field/coarsening/moments_coarsener.hpp"
#include "amr/types/amr_types.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"
#include "amr/data/field/coarsening/electric_field_coarsener.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_composite_refiner.hpp"
#include "amr/messengers/refinement_config.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/particles/particles_variable_fill_pattern.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"

#include <deque>

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/hier/CoarseFineBoundary.h>
#include <SAMRAI/xfer/BoxGeometryVariableFillPattern.h>

#include <memory>
#include <string>
#include <utility>
#include <iomanip>
#include <iostream>




namespace PHARE
{
namespace amr
{
    /** \brief An HybridMessenger is the specialization of a HybridMessengerStrategy for hybrid
     * to hybrid data communications.
     */
    template<typename HybridModel, typename RefinementParams>
    class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
    {
        using amr_types   = PHARE::amr::SAMRAI_Types;
        using level_t     = amr_types::level_t;
        using patch_t     = amr_types::patch_t;
        using hierarchy_t = amr_types::hierarchy_t;

        using Super = HybridMessengerStrategy<HybridModel>;

        using GridT             = HybridModel::grid_type;
        using IonsT             = HybridModel::ions_type;
        using ElectromagT       = HybridModel::electromag_type;
        using VecFieldT         = HybridModel::vecfield_type;
        using TensorFieldT      = IonsT::tensorfield_type;
        using GridLayoutT       = HybridModel::gridlayout_type;
        using FieldT            = VecFieldT::field_type;
        using VectorFieldDataT  = TensorFieldData<1, GridLayoutT, GridT, core::HybridQuantity>;
        using ResourcesManagerT = HybridModel::resources_manager_type;
        using IPhysicalModel    = HybridModel::Interface;

        static constexpr std::size_t dimension   = GridLayoutT::dimension;
        static constexpr std::size_t interpOrder = GridLayoutT::options.interp_order;

        using InteriorParticleRefineOp = RefinementParams::InteriorParticleRefineOp;
        using CoarseToFineRefineOpOld  = RefinementParams::CoarseToFineRefineOpOld;
        using CoarseToFineRefineOpNew  = RefinementParams::CoarseToFineRefineOpNew;

        using FieldTimeInterp = FieldLinearTimeInterpolate<GridLayoutT, GridT>;

        using VecFieldTimeInterp
            = VecFieldLinearTimeInterpolate<GridLayoutT, GridT, core::HybridQuantity>;

        template<typename Policy>
        using FieldCoarsenOp = FieldCoarsenOperator<GridLayoutT, GridT, Policy>;

        template<typename Policy>
        using VecFieldCoarsenOp
            = VecFieldCoarsenOperator<GridLayoutT, GridT, Policy, core::HybridQuantity>;

        using DefaultFieldCoarsenOp    = FieldCoarsenOp<DefaultFieldCoarsener<dimension>>;
        using DefaultVecFieldCoarsenOp = VecFieldCoarsenOp<DefaultFieldCoarsener<dimension>>;
        using MomentsFieldCoarsenOp    = FieldCoarsenOp<MomentsCoarsener<dimension>>;
        using MomentsVecFieldCoarsenOp = VecFieldCoarsenOp<MomentsCoarsener<dimension>>;
        using ElectricFieldCoarsenOp   = VecFieldCoarsenOp<ElectricFieldCoarsener<dimension>>;

    public:
        static inline std::string const stratName    = "HybridModel-HybridModel";
        static constexpr std::size_t rootLevelNumber = 0;


        HybridHybridMessengerStrategy(std::shared_ptr<ResourcesManagerT> const& manager,
                                      int const firstLevel,
                                      RefinementConfig const& refinementConfig = {})
            : HybridMessengerStrategy<HybridModel>{stratName}
            , resourcesManager_{manager}
            , firstLevel_{firstLevel}
        {
            makeRefineOperators_(refinementConfig);

            resourcesManager_->registerResources(Bold_);
            resourcesManager_->registerResources(Eold_);
            resourcesManager_->registerResources(NiOld_);
            resourcesManager_->registerResources(ViOld_);
            resourcesManager_->registerResources(sumTensor_);
        }

        virtual ~HybridHybridMessengerStrategy() = default;



        /* ------------------------------------------------------------------------
                        methods used for the IMessenger interface
           ------------------------------------------------------------------------ */


        /**
         * @brief allocate the messenger strategy internal variables to the model
         * resourceManager
         */
        void allocate(patch_t& patch, double const allocateTime) const override
        {
            resourcesManager_->allocate(Bold_, patch, allocateTime);
            resourcesManager_->allocate(Eold_, patch, allocateTime);
            resourcesManager_->allocate(NiOld_, patch, allocateTime);
            resourcesManager_->allocate(ViOld_, patch, allocateTime);
            resourcesManager_->allocate(sumTensor_, patch, allocateTime);
        }



        /**
         * @brief setup creates all SAMRAI algorithms to communicate data involved in a
         * messenger between the coarse and fine levels.
         *
         * This method creates the SAMRAI algorithms for communications associated between pairs
         * of variables. The function does not create the SAMRAI schedules since they depend on
         * the levels
         */
        void registerQuantities([[maybe_unused]] std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<HybridMessengerInfo> hybridInfo{
                dynamic_cast<HybridMessengerInfo*>(fromFinerInfo.release())};

            auto&& [b_id] = resourcesManager_->getIDsList(hybridInfo->modelMagnetic);

            magneticRefinePatchStrategy_.registerIDs(b_id);

            // we do not overwrite interior on patch ghost filling. In theory this doesn't matter
            // much since the only interior values are the outermost layer of faces of the domain,
            // and should be near equal from one patch to the other.
            BalgoPatchGhost.registerRefine(b_id, b_id, b_id, nullptr,
                                           nonOverwriteInteriorTFfillPattern);


            // for regrid, we need to overwrite the interior or else only the new ghosts would be
            // filled. We also need to use the regrid operator, which checks for nans before filling
            // the new values, as we do not want to overwrite the copy that was already done for the
            // faces that were already there before regrid.
            BregridAlgo.registerRefine(b_id, b_id, b_id, BRefineOp_,
                                       overwriteInteriorTFfillPattern);

            auto&& [e_id] = resourcesManager_->getIDsList(hybridInfo->modelElectric);


            EalgoPatchGhost.registerRefine(e_id, e_id, e_id, EfieldRefineOp_,
                                           nonOverwriteInteriorTFfillPattern);

            auto&& [e_reflux_id]  = resourcesManager_->getIDsList(hybridInfo->refluxElectric);
            auto&& [e_fluxsum_id] = resourcesManager_->getIDsList(hybridInfo->fluxSumElectric);


            RefluxAlgo.registerCoarsen(e_reflux_id, e_fluxsum_id, electricFieldCoarseningOp_);

            // we then need to refill the ghosts so that they agree with the newly refluxed cells

            PatchGhostRefluxedAlgo.registerRefine(e_reflux_id, e_reflux_id, e_reflux_id,
                                                  EfieldRefineOp_,
                                                  nonOverwriteInteriorTFfillPattern);

            registerGhostComms_(hybridInfo);
            registerInitComms(hybridInfo);
            registerSyncComms(hybridInfo);
        }

        void registerQuantities(IPhysicalModel& /*coarseModel*/, IPhysicalModel& fineModel) override
        {
            auto& hybridModel = dynamic_cast<HybridModel&>(fineModel);
            hybridModel.scheduler(Super::scheduler());
            ions_ = &hybridModel.state.ions; // for scheduler fill functions
        }


        /**
         * @brief all RefinerPool must be notified the level levelNumber now exist.
         * not doing so will result in communication to/from that level being impossible
         */
        void registerLevel(std::shared_ptr<hierarchy_t> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);


            magPatchGhostsRefineSchedules[levelNumber]
                = BalgoPatchGhost.createSchedule(level, &magneticRefinePatchStrategy_);

            elecPatchGhostsRefineSchedules[levelNumber] = EalgoPatchGhost.createSchedule(level);

            // technically not needed for finest as refluxing is not done onto it.
            patchGhostRefluxedSchedules[levelNumber] = PatchGhostRefluxedAlgo.createSchedule(level);

            auto& scheduler = Super::scheduler();

            for (auto const& view : electricGhostViews_)
                scheduler.template add<RefinerType::GhostField, ResourcesManagerT>(
                    view.name(), hierarchy, level);
            for (auto const& view : magneticGhostViews_)
                scheduler.template add<RefinerType::GhostField, ResourcesManagerT>(
                    view.name(), hierarchy, level);
            // chargeDensityLevelGhostsRefiners_.registerLevel(hierarchy, level);
            // velLevelGhostsRefiners_.registerLevel(hierarchy, level);
            scheduler.template add<RefinerType::ExteriorGhostParticles, ResourcesManagerT>(
                HybridScheduleTraits<HybridSchedule::IonGhostParticles>::key, hierarchy, level);

            chargeDensityPatchGhostsRefiners_.registerLevel(hierarchy, level);
            velPatchGhostsRefiners_.registerLevel(hierarchy, level);

            scheduler.template add<RefinerType::PatchIonPopBorderSum, ResourcesManagerT>(
                HybridScheduleTraits<HybridSchedule::IonPopBorderSum>::key, hierarchy, level);

            scheduler.template add<RefinerType::PatchIonPopBorderMax, ResourcesManagerT>(
                HybridScheduleTraits<HybridSchedule::IonBorderMax>::key, hierarchy, level);

            // root level is not initialized with a schedule using coarser level data
            // so we don't create these schedules if root level
            // TODO this 'if' may not be OK if L0 is regrided
            if (levelNumber != rootLevelNumber)
            {
                // refluxing
                auto const& coarseLevel      = hierarchy->getPatchLevel(levelNumber - 1);
                refluxSchedules[levelNumber] = RefluxAlgo.createSchedule(coarseLevel, level);

                // those are for refinement
                magInitRefineSchedules[levelNumber] = BalgoInit.createSchedule(
                    level, nullptr, levelNumber - 1, hierarchy, &magneticRefinePatchStrategy_);

                electricInitRefiners_.registerLevel(hierarchy, level);
                domainParticlesRefiners_.registerLevel(hierarchy, level);
                lvlGhostPartOldRefiners_.registerLevel(hierarchy, level);
                lvlGhostPartNewRefiners_.registerLevel(hierarchy, level);

                // and these for coarsening
                electroSynchronizers_.registerLevel(hierarchy, level);
                chargeDensitySynchronizers_.registerLevel(hierarchy, level);
                ionBulkVelSynchronizers_.registerLevel(hierarchy, level);
            }
        }



        /**
         * @brief regrid performs the regriding communications for Hybrid to Hybrid messengers
         , all quantities that are in initialization refiners need to be regridded
         */
        void regrid(std::shared_ptr<hierarchy_t> const& hierarchy, int const levelNumber,
                    std::shared_ptr<level_t> const& oldLevel, IPhysicalModel& model,
                    double const initDataTime) override
        {
            auto& hybridModel = dynamic_cast<HybridModel&>(model);
            auto level        = hierarchy->getPatchLevel(levelNumber);

            bool const isRegriddingL0 = levelNumber == 0 and oldLevel;

            magneticRegriding_(hierarchy, level, oldLevel, hybridModel, initDataTime);
            electricInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            domainParticlesRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);


            // we now call only levelGhostParticlesOld.fill() and not .regrid()
            // regrid() would refine from next coarser in regions of level not overlaping
            // oldLevel, but copy from domain particles of oldLevel where there is an
            // overlap while we do not a priori see why this could be wrong,but this led to
            // occasional failures of the SAMRAI MPI module. See
            // https://github.com/PHAREHUB/PHARE/issues/604 calling .fill() ensures that
            // levelGhostParticlesOld particles are filled exclusively from spliting next
            // coarser domain ones like when a new finest level is created.


            if (levelNumber != rootLevelNumber)
            {
                lvlGhostPartOldRefiners_.fill(levelNumber, initDataTime);
                copyLevelGhostOldToPushable_(*level, model);
            }

            // levelGhostNew will be refined in next firstStep

            // after filling the new level with the regrid schedule, some
            // nodes may not have been copied correctly, due to a bug in SAMRAI
            // it seems these nodes are only on ghost box border if that border
            // overlaps an old level patch border. See https://github.com/LLNL/SAMRAI/pull/293
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
         * @brief initLevel is used to initialize hybrid data on the level levelNumer at
         * time initDataTime from hybrid coarser data.
         */
        void initLevel(IPhysicalModel& model, level_t& level, double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();

            auto& hybridModel = static_cast<HybridModel&>(model);

            magInitRefineSchedules[levelNumber]->fillData(initDataTime);
            electricInitRefiners_.fill(levelNumber, initDataTime);

            // no need to call these :
            // magGhostsRefiners_.fill(levelNumber, initDataTime);
            // elecGhostsRefiners_.fill(levelNumber, initDataTime);
            // because the SAMRAI schedules in the 'init' communicators
            // already fill the patch ghost box from the neighbor interior box.
            // so ghost nodes are already filled .

            PHARE_LOG_START(3, "hybhybmessengerStrat::initLevel : interior part fill schedule");
            domainParticlesRefiners_.fill(levelNumber, initDataTime);
            PHARE_LOG_STOP(3, "hybhybmessengerStrat::initLevel : interior part fill schedule");

            lvlGhostPartOldRefiners_.fill(levelNumber, initDataTime);


            // levelGhostParticles will be pushed during the advance phase
            // they need to be identical to levelGhostParticlesOld before advance
            copyLevelGhostOldToPushable_(level, model);
        }



        /* ------------------------------------------------------------------------
                     methods used for the HybridMessenger interface
           ------------------------------------------------------------------------ */


        void fillMagneticGhosts(VecFieldT& B, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillMagneticGhosts");

            Super::scheduler().template fill<RefinerType::GhostField>(B.name(), level, fillTime);
        }

        void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillElectricGhosts");

            Super::scheduler().template fill<RefinerType::GhostField>(E.name(), level, fillTime);
        }




        /**
         * @brief fillIonGhostParticles will fill the interior ghost particle array from
         * neighbor patches of the same level. Before doing that, it empties the array for
         * all populations
         */
        void fillIonGhostParticles(IonsT& /*ions*/, level_t& level, double const fillTime) override
        {
            Super::scheduler().template fill<RefinerType::ExteriorGhostParticles>(
                HybridScheduleTraits<HybridSchedule::IonGhostParticles>::key, level, fillTime);
        }



        void fillIonPopBorders(IonsT& /*ions*/, level_t& level, double const fillTime) override
        {
            Super::scheduler().template fill<RefinerType::PatchIonPopBorderSum>(
                HybridScheduleTraits<HybridSchedule::IonPopBorderSum>::key, level, fillTime);
        }


        void fillIonBorders(IonsT& /*ions*/, level_t& level, double const fillTime) override
        {
            Super::scheduler().template fill<RefinerType::PatchIonPopBorderMax>(
                HybridScheduleTraits<HybridSchedule::IonBorderMax>::key, level, fillTime);
        }



        /**
         * @brief fillIonPopMomentGhosts works on moment ghost nodes
         *
         * level border nodes are completed by the deposition
         * of level ghost [old,new] particles for all populations, linear time interpolation
         * is used to get the contribution of old/new particles
         */
        void fillIonPopMomentGhosts(IonsT& ions, level_t& level,
                                    double const afterPushTime) override
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillIonPopMomentGhosts");

            auto alpha = timeInterpCoef_(afterPushTime, level.getLevelNumber());
            if (level.getLevelNumber() > 0 and (alpha < 0 or alpha > 1))
            {
                std::cout << std::setprecision(12) << alpha << "\n";
                throw std::runtime_error("ion moment ghost time interp coef invalid : alpha: "
                                         + std::to_string(alpha) + " beforePushTime "
                                         + std::to_string(afterPushTime) + " on level "
                                         + std::to_string(level.getLevelNumber()));
            }
            for (auto const& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                auto layout      = layoutFromPatch<GridLayoutT>(*patch);

                for (auto& pop : ions)
                {
                    auto& particleDensity = pop.particleDensity();
                    auto& chargeDensity   = pop.chargeDensity();
                    auto& flux            = pop.flux();

                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        auto& levelGhostOld = pop.levelGhostParticlesOld();
                        interpolate_(makeRange(levelGhostOld), particleDensity, chargeDensity, flux,
                                     layout, 1. - alpha);

                        auto& levelGhostNew = pop.levelGhostParticlesNew();
                        interpolate_(makeRange(levelGhostNew), particleDensity, chargeDensity, flux,
                                     layout, alpha);
                    }
                }
            }
        }


        /* pure (patch and level) ghost nodes are filled by applying a regular ghost
         * schedule i.e. that does not overwrite the border patch node previously well
         * calculated from particles Note : the ghost schedule only fills the total density
         * and bulk velocity and NOT population densities and fluxes. These partial moments
         * are already completed by the "sum" schedules (+= on incomplete nodes)*/
        // virtual void fillIonMomentGhosts(IonsT& ions, level_t& level,
        //                                  double const afterPushTime) override
        // {
        //     PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillIonMomentGhosts");
        //     auto& chargeDensity = ions.chargeDensity();
        //     auto& velocity      = ions.velocity();
        //
        //
        //     if (level.getLevelNumber() != 0)
        //     {
        //         setNaNsOnFieldGhosts(chargeDensity, level);
        //         setNaNsOnVecfieldGhosts(velocity, level);
        //     }
        //     chargeDensityLevelGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
        //     velLevelGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
        // }

        /**
         * @brief firstStep : in the HybridHybridMessengerStrategy, the firstStep method is
         * used to get level border ghost particles from the next coarser level. These
         * particles are defined in the future at the time the method is called because the
         * coarser level is ahead in time. These particles are communicated only at first
         * step of a substepping cycle. They will be used with the levelGhostParticlesOld
         * particles to get the moments on level border nodes. The method is does nothing if
         * the level is the root level because the root level cannot get levelGhost from
         * next coarser (it has none).
         */
        void firstStep(IPhysicalModel& /*model*/, level_t& level,
                       std::shared_ptr<hierarchy_t> const& /*hierarchy*/, double const currentTime,
                       double const prevCoarserTime, double const newCoarserTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::firstStep");

            auto levelNumber = level.getLevelNumber();
            if (newCoarserTime < prevCoarserTime)
                throw std::runtime_error(
                    "Error : prevCoarserTime (" + std::to_string(prevCoarserTime)
                    + ") should be < newCoarserTime (" + std::to_string(prevCoarserTime) + ")");

            // root level has no levelghost particles
            if (levelNumber != 0)
            {
                PHARE_LOG_START(3, "HybridHybridMessengerStrategy::firstStep.fill");
                lvlGhostPartNewRefiners_.fill(levelNumber, currentTime);
                PHARE_LOG_STOP(3, "HybridHybridMessengerStrategy::firstStep.fill");

                // during firstStep() coarser level and current level are at the same time
                // so 'time' is also the beforePushCoarseTime_
                beforePushCoarseTime_[levelNumber] = prevCoarserTime;
                afterPushCoarseTime_[levelNumber]  = newCoarserTime;
            }
        }


        /**
         * @brief lastStep is used to perform operations at the last step of a substepping
         * cycle. It is called after the level is advanced. Here for hybrid-hybrid messages,
         * the method moves levelGhostParticlesNew particles into levelGhostParticlesOld
         * ones. Then levelGhostParticlesNew are emptied since it will be filled again at
         * firstStep of the next substepping cycle. the new CoarseToFineOld content is then
         * copied to levelGhostParticles so that they can be pushed during the next subcycle
         */
        void lastStep(IPhysicalModel& model, level_t& level) override
        {
            if (level.getLevelNumber() == 0)
                return;

            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::lastStep");

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

                    std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
                    levelGhostParticlesNew.clear();
                    levelGhostParticles = levelGhostParticlesOld;
                }
            }
        }




        /**
         * @brief prepareStep is the concrete implementation of the
         * HybridMessengerStrategy::prepareStep method For hybrid-Hybrid communications.
         * This method copies the electromagnetic field, the density and bulk velocity,
         * defined at t=n. Since prepareStep() is called just before advancing the level,
         * this operation actually saves the t=n versions of B, E, Ni, Vi into the
         * messenger. When the time comes that the next finer level needs to time
         * interpolate these quantities at its ghost nodes, this level will be able to
         * interpolate at required time because the t=n fields of previous next coarser
         * step will be in the messenger.
         */
        void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::prepareStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, hybridModel.state.electromag, hybridModel.state.J,
                    hybridModel.state.ions, Bold_, Eold_, NiOld_, ViOld_);

                resourcesManager_->setTime(Bold_, *patch, currentTime);
                resourcesManager_->setTime(Eold_, *patch, currentTime);
                resourcesManager_->setTime(NiOld_, *patch, currentTime);
                resourcesManager_->setTime(ViOld_, *patch, currentTime);

                auto& Vi = hybridModel.state.ions.velocity();
                auto& Ni = hybridModel.state.ions.chargeDensity();

                Bold_.copyData(hybridModel.state.electromag.B);
                Eold_.copyData(hybridModel.state.electromag.E);
                ViOld_.copyData(Vi);
                NiOld_.copyData(Ni);
            }
        }




        void fillRootGhosts(IPhysicalModel& model, level_t& level,
                            double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();
            assert(levelNumber == 0);

            auto& hybridModel = static_cast<HybridModel&>(model);

            // TODO : still useful?
            Super::scheduler().template call<RefinerType::GhostField>(
                hybridModel.state.electromag.E.name(), levelNumber, initDataTime);

            // at some point in the future levelGhostParticles could be filled with injected
            // particles depending on the domain boundary condition.
            //
            // J ghosts need no filling here: J is computed locally by Ampere on
            // physical+1, exactly the layer the Ohm laplacian reads.
            //
            // do we need moment ghosts filled here?
            // a priori no because those are at this time only needed for coarsening, which
            // will not happen before the first advance
        }



        void synchronize(level_t& level) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::synchronize");

            auto levelNumber = level.getLevelNumber();
            PHARE_LOG_LINE_STR("synchronizing level " + std::to_string(levelNumber));

            // call coarsning schedules...
            electroSynchronizers_.sync(levelNumber);
            chargeDensitySynchronizers_.sync(levelNumber);
            ionBulkVelSynchronizers_.sync(levelNumber);
        }

        // this function coarsens the fluxSum onto the corresponding coarser fluxes (E in hybrid),
        // and fills the patch ghosts, making it ready for the faraday in the solver.reflux()
        void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                    double const syncTime) override
        {
            refluxSchedules[fineLevelNumber]->coarsenData();
            patchGhostRefluxedSchedules[coarserLevelNumber]->fillData(syncTime);
        }

        // after coarsening, domain nodes have been updated and therefore patch ghost nodes
        // will probably stop having the exact same value as their overlapped neighbor
        // note that NONE of the refiners used here overwrite interior nodes
        // because the point is to re-sync patch ghost nodes with overlaped domain nodes
        void postSynchronize(IPhysicalModel& model, level_t& level, double const time) override
        {
            auto levelNumber  = level.getLevelNumber();
            auto& hybridModel = static_cast<HybridModel&>(model);

            PHARE_LOG_LINE_STR("postSynchronize level " + std::to_string(levelNumber));

            // this electric schedule should probably only be a patch ghost one
            // since levelghost nodes are not affected by the coarsening
            Super::scheduler().template call<RefinerType::GhostField>(
                hybridModel.state.electromag.E.name(), levelNumber, time);

            chargeDensityPatchGhostsRefiners_.fill(levelNumber, time);
            velPatchGhostsRefiners_.fill(hybridModel.state.ions.velocity(), levelNumber, time);
        }

    private:
        // Select the field-refinement operators once at construction. The composite runtime
        // kernels are built from the configured order; B uses the stage-1 magnetic kernel (fills
        // every fine face; the magnetic patch strategy runs the stage-2 divB touch-up afterward).
        // Particle refine operators (interior / level-ghost) are NOT touched.
        void makeRefineOperators_(RefinementConfig const& config)
        {
            auto fieldKernel = [&] {
                return std::make_shared<KernelFieldRefineOperator<GridLayoutT, GridT>>(
                    makeRefineKernel<GridLayoutT, GridT>(config.order));
            };
            auto vecKernel = [&] {
                return std::make_shared<KernelVecFieldRefineOperator<VectorFieldDataT>>(
                    makeRefineKernel<GridLayoutT, GridT>(config.order));
            };
            auto magKernel = [&] {
                return std::make_shared<KernelVecFieldRefineOperator<VectorFieldDataT>>(
                    makeMagneticRefineKernel<GridLayoutT, GridT>(config.order));
            };

            fieldRefineOp_    = fieldKernel();
            vecFieldRefineOp_ = vecKernel();
            BInitRefineOp_    = magKernel();
            BRefineOp_        = magKernel();
            EfieldRefineOp_   = vecKernel();
        }

        void registerGhostComms_(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            // *********************************************************************
            // all of the ghost refiners take the nonOverwriteInteriorTFfillPattern as
            // domain nodes are already computed
            // *********************************************************************


            // model-E and Eavg C-F ghosts are time interpolated between the coarse E
            // snapshots Eold_ (t=n) and model E (t=n+1). Eavg is filled at the half-time
            // abscissa by the solver; model-E at newTime.
            for (auto const& ghost : info->ghostElectric)
                addGhostVecFieldFill_(electricGhostViews_, ghost, core::HybridQuantity::Vector::E,
                                      info->modelElectric, Eold_.name(), EfieldRefineOp_);


            // we need a separate patch strategy for each refiner so that each one can register
            // their required ids
            magneticPatchStratPerGhostRefiner_ = [&]() {
                std::vector<std::shared_ptr<MagneticRefinePatchStrategy<VectorFieldDataT>>>
                    result;

                result.reserve(info->ghostMagnetic.size());

                for (auto const& key : info->ghostMagnetic)
                {
                    auto&& [id] = resourcesManager_->getIDsList(key);

                    auto patch_strat
                        = std::make_shared<MagneticRefinePatchStrategy<VectorFieldDataT>>();

                    patch_strat->registerIDs(id);

                    result.push_back(patch_strat);
                }
                return result;
            }();

            // B (and Bpred) C-F ghosts are time interpolated between Bold_ (t=n) and model
            // B (t=n+1). Per-refiner magnetic patch strategy is carried through so the
            // divB touch-up still runs after spatial refine.
            for (size_t i = 0; i < info->ghostMagnetic.size(); ++i)
                addGhostVecFieldFill_(magneticGhostViews_, info->ghostMagnetic[i],
                                      core::HybridQuantity::Vector::B, info->modelMagnetic,
                                      Bold_.name(), BRefineOp_,
                                      magneticPatchStratPerGhostRefiner_[i].get());


            // J has no refiner: it is never communicated, each level computes it
            // locally via Ampere on physical+1 from the time-refined B ghosts.

            // chargeDensityLevelGhostsRefiners_.addTimeRefiner(
            //     info->modelIonDensity, info->modelIonDensity, NiOld_.name(), fieldRefineOp_,
            //     fieldTimeOp_, info->modelIonDensity, nonOverwriteInteriorFieldFillPattern);
            //
            //
            // velLevelGhostsRefiners_.addTimeRefiners(
            //     info->ghostBulkVelocity, info->modelIonBulkVelocity, ViOld_.name(),
            //     vecFieldRefineOp_, vecFieldTimeOp_, nonOverwriteInteriorTFfillPattern);

            chargeDensityPatchGhostsRefiners_.addTimeRefiner(
                info->modelIonDensity, info->modelIonDensity, NiOld_.name(), fieldRefineOp_,
                fieldTimeOp_, info->modelIonDensity, nonOverwriteInteriorFieldFillPattern);

            velPatchGhostsRefiners_.addTimeRefiners(
                info->ghostBulkVelocity, info->modelIonBulkVelocity, ViOld_.name(),
                vecFieldRefineOp_, vecFieldTimeOp_, nonOverwriteInteriorTFfillPattern);
        }




        void registerInitComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto& scheduler = Super::scheduler();

            auto b_id = resourcesManager_->getID(info->modelMagnetic);
            BalgoInit.registerRefine(*b_id, *b_id, *b_id, BInitRefineOp_,
                                     overwriteInteriorTFfillPattern);

            // no fill pattern given for this init
            // will use boxgeometryvariable fillpattern, itself using the
            // field geometry with overwrite_interior true from SAMRAI
            // we could set the overwriteInteriorTFfillPattern it would be the same
            electricInitRefiners_.addStaticRefiners(info->initElectric, EfieldRefineOp_,
                                                    info->initElectric);


            domainParticlesRefiners_.addStaticRefiners(
                info->interiorParticles, interiorParticleRefineOp_, info->interiorParticles);


            lvlGhostPartOldRefiners_.addStaticRefiners(info->levelGhostParticlesOld,
                                                       levelGhostParticlesOldOp_,
                                                       info->levelGhostParticlesOld);


            lvlGhostPartNewRefiners_.addStaticRefiners(info->levelGhostParticlesNew,
                                                       levelGhostParticlesNewOp_,
                                                       info->levelGhostParticlesNew);


            // one algorithm per population, the population patch ghost particles are moved into
            // the domain of the patch they now belong to, then emptied
            for (std::size_t i = 0; i < info->patchGhostParticles.size(); ++i)
            {
                using Schedule = HybridScheduleTraits<HybridSchedule::IonGhostParticles>;
                auto const idx
                    = scheduler.template n_algorithms<Schedule::type>(Schedule::key);
                scheduler
                    .template add_algorithm<Schedule::type>(
                        Schedule::key,
                        [this, i, idx](auto const& level, double const fillTime) {
                            PHARE_LOG_SCOPE(
                                1, "HybridHybridMessengerStrategy::fillIonGhostParticles");
                            Super::scheduler().template call<Schedule::type>(
                                Schedule::key, level, fillTime, idx);

                            for (auto patch : resourcesManager_->enumerate(level, *ions_))
                                (*ions_)[i].patchGhostParticles().clear();
                        })
                    .template register_resource<Schedule::type>(
                        Schedule::key, resourcesManager_, info->patchGhostParticles[i],
                        info->patchGhostParticles[i], info->patchGhostParticles[i], nullptr,
                        std::make_shared<ParticleDomainFromGhostFillPattern<GridLayoutT>>());
            }


            // flux (VecField/TensorFieldData) and both density fields (FieldData) are registered
            // onto the SAME algorithm/schedule per population so they are communicated together
            // in a single set of MPI messages (see MultiFieldBorderOpTransactionFactory).
            assert(info->sumBorderFields.size() == 2 * info->ghostFlux.size());
            for (std::size_t i = 0; i < info->ghostFlux.size(); ++i)
            {
                // we cannot have the schedule doing the += in place in the flux/density arrays
                // because some overlaps could be counted several times.
                // we therefore first copy flux and both density fields into scratch buffers and
                // then execute the (fused) schedule onto them before copying the result back
                // onto the real arrays.
                using Schedule = HybridScheduleTraits<HybridSchedule::IonPopBorderSum>;
                auto const idx
                    = scheduler.template n_algorithms<Schedule::type>(Schedule::key);
                scheduler
                    .template add_algorithm<Schedule::type>(
                        Schedule::key,
                        [this, i, idx](auto const& level, double const fillTime) {
                            auto& ions = *ions_;
                            for (auto patch : resourcesManager_->enumerate(
                                     level, ions, tmpVecField0_, tmpField0_, tmpField1_))
                                copyIonPopBorderScratch_(ions[i], /*toScratch=*/true);

                            Super::scheduler().template call<Schedule::type>(
                                Schedule::key, level, fillTime, idx);

                            for (auto patch : resourcesManager_->enumerate(
                                     level, ions, tmpVecField0_, tmpField0_, tmpField1_))
                                copyIonPopBorderScratch_(ions[i], /*toScratch=*/false);
                        })
                    .template register_resource<Schedule::type>(
                        Schedule::key, resourcesManager_, tmpVecField0_.name(), info->ghostFlux[i],
                        tmpVecField0_.name(), nullptr,
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>())
                    .template register_resource<Schedule::type>(
                        Schedule::key, resourcesManager_, tmpField0_.name(),
                        info->sumBorderFields[2 * i], tmpField0_.name(), nullptr,
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>())
                    .template register_resource<Schedule::type>(
                        Schedule::key, resourcesManager_, tmpField1_.name(),
                        info->sumBorderFields[2 * i + 1], tmpField1_.name(), nullptr,
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());
            }


            // mass density, charge density (FieldData) and bulk velocity (VecFieldData) are
            // registered onto the SAME algorithm/schedule so they are communicated together in
            // a single set of MPI messages (see MultiFieldBorderOpTransactionFactory). Max is
            // idempotent (max(x,x) == x) so, unlike the sum refiners above, these operate
            // directly on the real fields with no scratch buffer needed.
            using MaxSchedule = HybridScheduleTraits<HybridSchedule::IonBorderMax>;
            scheduler.template add_algorithm<MaxSchedule::type>(MaxSchedule::key);

            assert(info->maxBorderFields.size() == 2); // mass & charge densities
            for (auto const& field : info->maxBorderFields)
                scheduler.template register_resource<MaxSchedule::type>(
                    MaxSchedule::key, resourcesManager_, field, field, field, nullptr,
                    std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());

            assert(info->maxBorderVecFields.size() == 1);
            for (auto const& vecfield : info->maxBorderVecFields)
                scheduler.template register_resource<MaxSchedule::type>(
                    MaxSchedule::key, resourcesManager_, vecfield, vecfield, vecfield, nullptr,
                    std::make_shared<
                        TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
        }



        void registerSyncComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            electroSynchronizers_.add(info->modelElectric, electricFieldCoarseningOp_,
                                      info->modelElectric);

            ionBulkVelSynchronizers_.add(info->modelIonBulkVelocity, vecFieldMomentsCoarseningOp_,
                                         info->modelIonBulkVelocity);

            chargeDensitySynchronizers_.add(info->modelIonDensity, fieldMomentsCoarseningOp_,
                                            info->modelIonDensity);
        }




        //! one time interpolated GhostField algorithm for `ghost`, filled by its name, its ghosts
        //! are NaN-ed before each fill through a view kept in `views`
        void addGhostVecFieldFill_(std::deque<VecFieldT>& views, std::string const& ghost,
                                   auto const quantity, std::string const& model,
                                   std::string const& old, auto const& refineOp,
                                   SAMRAI::xfer::RefinePatchStrategy* patchStrat = nullptr)
        {
            auto& view = views.emplace_back(ghost, quantity);
            Super::scheduler()
                .template add_algorithm<RefinerType::GhostField>(
                    ghost,
                    [this, &view](auto const& level, double const fillTime) {
                        PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillGhostVecField");
                        setNaNsOnVecfieldGhosts(view, level);
                        Super::scheduler().template call<RefinerType::GhostField>(
                            view.name(), level, fillTime, 0);
                    },
                    patchStrat)
                .template register_time_interpolated_resource<RefinerType::GhostField>(
                    ghost, resourcesManager_, ghost, ghost, old, model, refineOp, vecFieldTimeOp_,
                    nonOverwriteInteriorTFfillPattern);
        }


        //! copies a population's flux and densities to (or back from) the scratch fields
        //! the fused border sum schedule operates on, see fillIonPopBorders()
        void copyIonPopBorderScratch_(auto& pop, bool const toScratch)
        {
            auto const copy = [&](auto& field, auto& scratch) {
                auto const bytes = field.size() * sizeof(*field.data());
                if (toScratch)
                    std::memcpy(scratch.data(), field.data(), bytes);
                else
                    std::memcpy(field.data(), scratch.data(), bytes);
            };

            auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();
            for (std::uint8_t c = 0; c < N; ++c)
                copy(pop.flux()[c], tmpVecField0_[c]);
            copy(pop.particleDensity(), tmpField0_);
            copy(pop.chargeDensity(), tmpField1_);
        }


        void copyLevelGhostOldToPushable_(level_t& level, IPhysicalModel& model)
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

                    levelGhostParticles = levelGhostParticlesOld;
                }
            }
        }




        double timeInterpCoef_(double const afterPushTime, std::size_t levelNumber)
        {
            return (afterPushTime - beforePushCoarseTime_[levelNumber])
                   / (afterPushCoarseTime_[levelNumber] - beforePushCoarseTime_[levelNumber]);
        }




        void magneticRegriding_(std::shared_ptr<hierarchy_t> const& hierarchy,
                                std::shared_ptr<level_t> const& level,
                                std::shared_ptr<level_t> const& oldLevel, HybridModel& hybridModel,
                                double const initDataTime)
        {
            auto magSchedule = BregridAlgo.createSchedule(
                level, oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy,
                &magneticRefinePatchStrategy_);
            magSchedule->fillData(initDataTime);
        }


        /** * @brief setNaNsOnFieldGhosts sets NaNs on all the ghost nodes of the field on this
         * patch.
         *
         * The magnetic/electric ghost fills (see addGhostVecFieldFill_) run one
         * SAMRAI RefineSchedule per fill. It first copies from the interiors of the other
         * patches of the same level (patch ghosts, periodic images included), then refines from
         * the coarser level only where no patch of this level overlaps (level ghosts). The refine
         * operators only write nodes that are still NaN, so the NaNs left after the copy are
         * exactly the level ghost nodes, and the nodes the copy already filled are not
         * overwritten. Setting NaNs on the patch ghosts too is therefore harmless: the copy
         * fills them before the refinement reads anything.
         *
         * This is needed because the schedule copy is done before refinement
         * as a result of FieldVariable::fineBoundaryRepresentsVariable=false
         */
        void setNaNsOnFieldGhosts(FieldT& field)
        {
            auto const nGhosts = GridLayoutT::nDNbrGhosts(core::QtyCentering::primal)[0];
            field[core::NdArrayMask{0, nGhosts - 1}]
                = std::numeric_limits<typename FieldT::type>::quiet_NaN();
        }

        void setNaNsOnFieldGhosts(FieldT& field, level_t const& level)
        {
            for ([[maybe_unused]] auto& patch : resourcesManager_->enumerate(level, field))
                setNaNsOnFieldGhosts(field);
        }

        void setNaNsOnVecfieldGhosts(VecFieldT& vf, level_t const& level)
        {
            for ([[maybe_unused]] auto& patch : resourcesManager_->enumerate(level, vf))
                for (auto& field : vf)
                    setNaNsOnFieldGhosts(field);
        }

        VecFieldT Bold_{stratName + "_Bold", core::HybridQuantity::Vector::B};
        VecFieldT Eold_{stratName + "_Eold", core::HybridQuantity::Vector::E};
        VecFieldT ViOld_{stratName + "_VBulkOld", core::HybridQuantity::Vector::V};
        FieldT NiOld_{stratName + "_NiOld", core::HybridQuantity::Scalar::rho};

        TensorFieldT sumTensor_{"PHARE_sumTensor", core::HybridQuantity::Tensor::M};

        //! shared temporaries declared by HybridModel; registered/allocated there, used here
        //! as scratch space for the fused ion population border-sum schedule (flux + both
        //! density fields communicated together, see fillIonPopBorders())
        VecFieldT tmpVecField0_{"tmpVecField0", core::HybridQuantity::Vector::V};
        FieldT tmpField0_{"tmpField0", core::HybridQuantity::all_primal_field};
        FieldT tmpField1_{"tmpField1", core::HybridQuantity::all_primal_field};



        //! ResourceManager shared with other objects (like the HybridModel)
        std::shared_ptr<ResourcesManagerT> resourcesManager_;


        int const firstLevel_;
        std::unordered_map<std::size_t, double> beforePushCoarseTime_;
        std::unordered_map<std::size_t, double> afterPushCoarseTime_;

        core::Interpolator<dimension, interpOrder> interpolate_;

        using rm_t                    = ResourcesManagerT;
        using RefineOperator          = SAMRAI::hier::RefineOperator;
        using TimeInterpolateOperator = SAMRAI::hier::TimeInterpolateOperator;

        // these refiners are used to initialize electromagnetic fields when creating
        // a new level (initLevel) or regridding (regrid)
        using InitRefinerPool             = RefinerPool<rm_t, RefinerType::InitField>;
        using GhostRefinerPool            = RefinerPool<rm_t, RefinerType::GhostField>;
        using InitDomPartRefinerPool      = RefinerPool<rm_t, RefinerType::InitInteriorPart>;
        using LevelBorderFieldRefinerPool = RefinerPool<rm_t, RefinerType::LevelBorderField>;
        using DomainGhostPartRefinerPool  = RefinerPool<rm_t, RefinerType::ExteriorGhostParticles>;
        using PatchGhostRefinerPool       = RefinerPool<rm_t, RefinerType::PatchGhostField>;
        using FieldFillPattern_t          = FieldFillPattern<dimension>;
        using TensorFieldFillPattern_t    = TensorFieldFillPattern<dimension /*, rank=1*/>;


        InitRefinerPool electricInitRefiners_{resourcesManager_};


        SAMRAI::xfer::RefineAlgorithm BalgoPatchGhost;
        SAMRAI::xfer::RefineAlgorithm BalgoInit;
        SAMRAI::xfer::RefineAlgorithm BregridAlgo;
        SAMRAI::xfer::RefineAlgorithm EalgoPatchGhost;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magInitRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> BpredGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> elecPatchGhostsRefineSchedules;

        SAMRAI::xfer::CoarsenAlgorithm RefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::RefineAlgorithm PatchGhostRefluxedAlgo;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> refluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> patchGhostRefluxedSchedules;

        //! views on the electric/magnetic fields that need ghosts to be filled, see
        //! addGhostVecFieldFill_
        std::deque<VecFieldT> electricGhostViews_;
        std::deque<VecFieldT> magneticGhostViews_;

        //! the fine model ions, used by the scheduler fill functions
        IonsT* ions_ = nullptr;


        // moment ghosts
        // The border node is already complete by the deposit of ghost particles
        // these refiners are used to fill ghost nodes, and therefore, owing to
        // the GhostField tag, will only assign pure ghost nodes. Border nodes will
        // be overwritten only on level borders, which does not seem to be an issue.
        // ******
        // NOTE :
        // *****
        // these and all the code that use them is commented
        // the reason for not deleting is that in its current state the code
        // only deposits levelghost particles which therefore leaves some of the level
        // ghost nodes incomplete (missing the outside contribution).
        // We thought about replacing levelghost particle deposit by filling a level ghost schedule
        // but at interp order >=2, levelghost particles will contribute to inner domain nodes
        // which a schedule will not do so we need them.
        // Keeping this code here is a way to ease the filling of pure level ghost nodes
        // if we decide to do so one day. This would overwrite what level ghost particles
        // have deposited on level ghost nodes, but since it is incomplete it does not matter
        // on the other hand this would be necessary if we wanted to have a multiple point
        // coarsening operator for moments.
        // LevelBorderFieldRefinerPool chargeDensityLevelGhostsRefiners_{resourcesManager_};
        // LevelBorderFieldRefinerPool velLevelGhostsRefiners_{resourcesManager_};

        PatchGhostRefinerPool chargeDensityPatchGhostsRefiners_{resourcesManager_};
        PatchGhostRefinerPool velPatchGhostsRefiners_{resourcesManager_};

        // pool of refiners for interior particles of each population
        // and the associated refinement operator
        InitDomPartRefinerPool domainParticlesRefiners_{resourcesManager_};

        using RefOp_ptr = std::shared_ptr<RefineOperator>;

        RefOp_ptr interiorParticleRefineOp_{std::make_shared<InteriorParticleRefineOp>()};

        //! store communicators for coarse to fine particles old
        // pools of refiners to fill level ghost particles, old and new ones
        // and their associated refinement operator
        static auto constexpr LGRefT = RefinerType::LevelBorderParticles;
        RefinerPool<rm_t, LGRefT> lvlGhostPartOldRefiners_{resourcesManager_};
        RefinerPool<rm_t, LGRefT> lvlGhostPartNewRefiners_{resourcesManager_};
        RefOp_ptr levelGhostParticlesOldOp_{std::make_shared<CoarseToFineRefineOpOld>()};
        RefOp_ptr levelGhostParticlesNewOp_{std::make_shared<CoarseToFineRefineOpNew>()};


        //! to grab particle leaving neighboring patches and inject into domain

        SynchronizerPool<rm_t> chargeDensitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> ionBulkVelSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> electroSynchronizers_{resourcesManager_};


        // built in the ctor body (makeRefineOperators_) from the composite runtime kernels.
        RefOp_ptr fieldRefineOp_;
        RefOp_ptr vecFieldRefineOp_;


        RefOp_ptr BInitRefineOp_;
        RefOp_ptr BRefineOp_;
        RefOp_ptr EfieldRefineOp_;
        std::shared_ptr<FieldFillPattern_t> nonOverwriteInteriorFieldFillPattern
            = std::make_shared<FieldFillPattern<dimension>>(); // stateless (mostly)

        std::shared_ptr<FieldFillPattern_t> overwriteInteriorFieldFillPattern
            = std::make_shared<FieldFillPattern<dimension>>(
                /*overwrite_interior=*/true); // stateless (mostly)

        std::shared_ptr<TensorFieldFillPattern_t> nonOverwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>();

        std::shared_ptr<TensorFieldFillPattern_t> overwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>(
                /*overwrite_interior=*/true);

        std::shared_ptr<TimeInterpolateOperator> fieldTimeOp_{std::make_shared<FieldTimeInterp>()};
        std::shared_ptr<TimeInterpolateOperator> vecFieldTimeOp_{
            std::make_shared<VecFieldTimeInterp>()};

        using CoarsenOperator_ptr = std::shared_ptr<SAMRAI::hier::CoarsenOperator>;

        CoarsenOperator_ptr fieldMomentsCoarseningOp_{std::make_shared<MomentsFieldCoarsenOp>()};
        CoarsenOperator_ptr vecFieldMomentsCoarseningOp_{
            std::make_shared<MomentsVecFieldCoarsenOp>()};
        CoarsenOperator_ptr electricFieldCoarseningOp_{std::make_shared<ElectricFieldCoarsenOp>()};

        MagneticRefinePatchStrategy<VectorFieldDataT> magneticRefinePatchStrategy_{};

        std::vector<std::shared_ptr<MagneticRefinePatchStrategy<VectorFieldDataT>>>
            magneticPatchStratPerGhostRefiner_;
    };


} // namespace amr

} // namespace PHARE

#endif
