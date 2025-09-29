#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP

#include "core/def.hpp" // IWYU pragma: keep
#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep


#include "core/hybrid/hybrid_quantities.hpp"
#include "core/numerics/interpolator/interpolator.hpp"

#include "refiner_pool.hpp"
#include "synchronizer_pool.hpp"
#include "amr/types/amr_types.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/refine/field_refiner.hpp"
#include "amr/data/field/refine/field_moments_refiner.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"

#include "amr/data/field/coarsening/electric_field_coarsener.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/electric_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_regrider.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/particles/particles_variable_fill_pattern.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/resources_manager/amr_utils.hpp"

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
        static constexpr std::size_t interpOrder = GridLayoutT::interp_order;

        using InteriorParticleRefineOp = RefinementParams::InteriorParticleRefineOp;
        using CoarseToFineRefineOpOld  = RefinementParams::CoarseToFineRefineOpOld;
        using CoarseToFineRefineOpNew  = RefinementParams::CoarseToFineRefineOpNew;

        template<typename Policy>
        using FieldRefineOp = FieldRefineOperator<GridLayoutT, GridT, Policy>;

        template<typename Policy>
        using VecFieldRefineOp = VecFieldRefineOperator<GridLayoutT, GridT, Policy>;

        using DefaultFieldRefineOp    = FieldRefineOp<DefaultFieldRefiner<dimension>>;
        using DefaultVecFieldRefineOp = VecFieldRefineOp<DefaultFieldRefiner<dimension>>;
        using FieldMomentsRefineOp    = FieldRefineOp<FieldMomentsRefiner<dimension>>;
        using VecFieldMomentsRefineOp = VecFieldRefineOp<FieldMomentsRefiner<dimension>>;
        using MagneticFieldRefineOp   = VecFieldRefineOp<MagneticFieldRefiner<dimension>>;
        using MagneticFieldRegridOp   = VecFieldRefineOp<MagneticFieldRegrider<dimension>>;
        using ElectricFieldRefineOp   = VecFieldRefineOp<ElectricFieldRefiner<dimension>>;
        using FieldTimeInterp         = FieldLinearTimeInterpolate<GridLayoutT, GridT>;

        using VecFieldTimeInterp
            = VecFieldLinearTimeInterpolate<GridLayoutT, GridT, core::HybridQuantity>;

        template<typename Policy>
        using FieldCoarsenOp = FieldCoarsenOperator<GridLayoutT, GridT, Policy>;

        template<typename Policy>
        using VecFieldCoarsenOp
            = VecFieldCoarsenOperator<GridLayoutT, GridT, Policy, core::HybridQuantity>;

        using DefaultFieldCoarsenOp    = FieldCoarsenOp<DefaultFieldCoarsener<dimension>>;
        using DefaultVecFieldCoarsenOp = VecFieldCoarsenOp<DefaultFieldCoarsener<dimension>>;
        using ElectricFieldCoarsenOp   = VecFieldCoarsenOp<ElectricFieldCoarsener<dimension>>;

    public:
        static inline std::string const stratName    = "HybridModel-HybridModel";
        static constexpr std::size_t rootLevelNumber = 0;


        HybridHybridMessengerStrategy(std::shared_ptr<ResourcesManagerT> manager,
                                      int const firstLevel)
            : HybridMessengerStrategy<HybridModel>{stratName}
            , resourcesManager_{std::move(manager)}
            , firstLevel_{firstLevel}
        {
            resourcesManager_->registerResources(Jold_);
            resourcesManager_->registerResources(NiOld_);
            resourcesManager_->registerResources(ViOld_);
            resourcesManager_->registerResources(sumVec_);
            resourcesManager_->registerResources(sumField_);
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
            resourcesManager_->allocate(Jold_, patch, allocateTime);
            resourcesManager_->allocate(NiOld_, patch, allocateTime);
            resourcesManager_->allocate(ViOld_, patch, allocateTime);
            resourcesManager_->allocate(sumVec_, patch, allocateTime);
            resourcesManager_->allocate(sumField_, patch, allocateTime);
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

            // auto&& [b_model, b_pred] =
            // resourcesManager_->getIDsList(hybridInfo->ghostMagnetic[0],
            //                                                          hybridInfo->ghostMagnetic[1]);


            magneticRefinePatchStrategy_.registerIDs(b_id);
            // BpredRefinePatchStrategy_.registerIDs(b_pred);

            // we do not overwrite interior on patch ghost filling. In theory this doesn't matter
            // much since the only interior values are the outermost layer of faces of the domain,
            // and should be near equal from one patch to the other.
            BalgoPatchGhost.registerRefine(b_id, b_id, b_id, BfieldRefineOp_,
                                           nonOverwriteInteriorTFfillPattern);


            // for regrid, we need to overwrite the interior or else only the new ghosts would be
            // filled. We also need to use the regrid operator, which checks for nans before filling
            // the new values, as we do not want to overwrite the copy that was already done for the
            // faces that were already there before regrid.
            BregridAlgo.registerRefine(b_id, b_id, b_id, BfieldRegridOp_,
                                       overwriteInteriorTFfillPattern);

            // this is a bit ugly, should be refactored asap also fills ghosts for both each time
            // which is not great
            // BghostAlgo.registerRefine(b_model, b_model, b_model, BfieldRegridOp_,
            //                           overwriteInteriorTFfillPattern);
            // BPredGhostAlgo.registerRefine(b_pred, b_pred, b_pred, BfieldRegridOp_,
            //                               overwriteInteriorTFfillPattern);

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

            // magGhostsRefineSchedules[levelNumber]
            //     = BghostAlgo.createSchedule(level, level->getNextCoarserHierarchyLevelNumber(),
            //                                 hierarchy, &magneticRefinePatchStrategy_);
            //
            // BpredGhostsRefineSchedules[levelNumber]
            //     = BPredGhostAlgo.createSchedule(level,
            //     level->getNextCoarserHierarchyLevelNumber(),
            //                                     hierarchy, &BpredRefinePatchStrategy_);

            elecPatchGhostsRefineSchedules[levelNumber] = EalgoPatchGhost.createSchedule(level);

            // technically not needed for finest as refluxing is not done onto it.
            patchGhostRefluxedSchedules[levelNumber] = PatchGhostRefluxedAlgo.createSchedule(level);

            elecGhostsRefiners_.registerLevel(hierarchy, level);
            currentGhostsRefiners_.registerLevel(hierarchy, level);
            chargeDensityGhostsRefiners_.registerLevel(hierarchy, level);
            velGhostsRefiners_.registerLevel(hierarchy, level);
            domainGhostPartRefiners_.registerLevel(hierarchy, level);

            chargeDensityPatchGhostsRefiners_.registerLevel(hierarchy, level);
            velPatchGhostsRefiners_.registerLevel(hierarchy, level);

            for (auto& refiner : popFluxBorderSumRefiners_)
                refiner.registerLevel(hierarchy, level);

            for (auto& refiner : popDensityBorderSumRefiners_)
                refiner.registerLevel(hierarchy, level);

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

            // computeIonMoments_(*level, model);
            // levelGhostNew will be refined in next firstStep

            // after filling the new level with the regrid schedule, some
            // nodes may not have been copied correctly, due to a bug in SAMRAI
            // it seems these nodes are only on ghost box border if that border
            // overlaps an old level patch border. See https://github.com/LLNL/SAMRAI/pull/293

            // magPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);
            // elecPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);
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
            // computeIonMoments_(level, model);
        }



        /* ------------------------------------------------------------------------
                     methods used for the HybridMessenger interface
           ------------------------------------------------------------------------ */


        void fillMagneticGhosts(VecFieldT& B, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillMagneticGhosts");

            // setNaNsOnVecfieldGhosts(B, level);
            // if (B.name() == "EM_B")
            //     magGhostsRefineSchedules[level.getLevelNumber()]->fillData(fillTime);
            // else if (B.name() == "EMPred_B")
            //     BpredGhostsRefineSchedules[level.getLevelNumber()]->fillData(fillTime);
            // else
            //     throw std::runtime_error("unknown magnetic field name : " + B.name());
        }

        void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillElectricGhosts");

            setNaNsOnVecfieldGhosts(E, level);
            elecGhostsRefiners_.fill(E, level.getLevelNumber(), fillTime);
        }




        void fillCurrentGhosts(VecFieldT& J, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillCurrentGhosts");
            setNaNsOnVecfieldGhosts(J, level);
            currentGhostsRefiners_.fill(J, level.getLevelNumber(), fillTime);
        }




        /**
         * @brief fillIonGhostParticles will fill the interior ghost particle array from
         * neighbor patches of the same level. Before doing that, it empties the array for
         * all populations
         */
        void fillIonGhostParticles(IonsT& ions, level_t& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillIonGhostParticles");

            domainGhostPartRefiners_.fill(level.getLevelNumber(), fillTime);

            for (auto patch : resourcesManager_->enumerate(level, ions))
                for (auto& pop : ions)
                    pop.patchGhostParticles().clear();
        }



        void fillFluxBorders(IonsT& ions, level_t& level, double const fillTime) override
        {
            auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();
            using value_type = FieldT::value_type;


            // we cannot have the schedule doign the += in place in the flux array
            // because some overlaps could be counted several times.
            // we therefore first copy flux into a sumVec buffer and then
            // execute the schedule onto that before copying it back onto the flux array
            for (std::size_t i = 0; i < ions.size(); ++i)
            {
                for (auto patch : resourcesManager_->enumerate(level, ions, sumVec_))
                    for (std::uint8_t c = 0; c < N; ++c)
                        std::memcpy(sumVec_[c].data(), ions[i].flux()[c].data(),
                                    ions[i].flux()[c].size() * sizeof(value_type));


                popFluxBorderSumRefiners_[i].fill(level.getLevelNumber(), fillTime);

                for (auto patch : resourcesManager_->enumerate(level, ions, sumVec_))
                    for (std::uint8_t c = 0; c < N; ++c)
                        std::memcpy(ions[i].flux()[c].data(), sumVec_[c].data(),
                                    ions[i].flux()[c].size() * sizeof(value_type));
            }
        }

        void fillDensityBorders(IonsT& ions, level_t& level, double const fillTime) override
        {
            using value_type = FieldT::value_type;

            std::size_t const fieldsPerPop = popDensityBorderSumRefiners_.size() / ions.size();

            for (std::size_t i = 0; i < ions.size(); ++i)
            {
                for (auto patch : resourcesManager_->enumerate(level, ions, sumField_))
                    std::memcpy(sumField_.data(), ions[i].particleDensity().data(),
                                ions[i].particleDensity().size() * sizeof(value_type));


                popDensityBorderSumRefiners_[i * fieldsPerPop].fill(level.getLevelNumber(),
                                                                    fillTime);

                for (auto patch : resourcesManager_->enumerate(level, ions, sumField_))
                    std::memcpy(ions[i].particleDensity().data(), sumField_.data(),
                                ions[i].particleDensity().size() * sizeof(value_type));

                //

                for (auto patch : resourcesManager_->enumerate(level, ions, sumField_))
                    std::memcpy(sumField_.data(), ions[i].chargeDensity().data(),
                                ions[i].chargeDensity().size() * sizeof(value_type));

                popDensityBorderSumRefiners_[i * fieldsPerPop + 1].fill(level.getLevelNumber(),
                                                                        fillTime);

                for (auto patch : resourcesManager_->enumerate(level, ions, sumField_))
                    std::memcpy(ions[i].chargeDensity().data(), sumField_.data(),
                                ions[i].chargeDensity().size() * sizeof(value_type));
            }
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
                    // first thing to do is to project patchGhostParitcles moments


                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // then grab levelGhostParticlesOld and levelGhostParticlesNew
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
        virtual void fillIonMomentGhosts(IonsT& ions, level_t& level,
                                         double const afterPushTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillIonMomentGhosts");
            auto& chargeDensity = ions.chargeDensity();
            auto& velocity      = ions.velocity();
            chargeDensityGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
            velGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
        }

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
         * This method copies the  density J, and the density and bulk velocity, defined at t=n.
         * Since prepareStep() is called just before advancing the level, this operation
         * actually saves the t=n versions of J, Ni, Vi into the messenger. When the time comes
         * that the next finer level needs to time interpolate the electromagnetic field and
         * current at its ghost nodes, this level will be able to interpolate at required time
         * because the t=n Vi,Ni,J fields of previous next coarser step will be in the
         * messenger.
         */
        void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::prepareStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, hybridModel.state.electromag, hybridModel.state.J,
                    hybridModel.state.ions, Jold_, NiOld_, ViOld_);

                resourcesManager_->setTime(Jold_, *patch, currentTime);
                resourcesManager_->setTime(NiOld_, *patch, currentTime);
                resourcesManager_->setTime(ViOld_, *patch, currentTime);

                auto& J  = hybridModel.state.J;
                auto& Vi = hybridModel.state.ions.velocity();
                auto& Ni = hybridModel.state.ions.chargeDensity();

                Jold_.copyData(J);
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

            elecGhostsRefiners_.fill(hybridModel.state.electromag.E, levelNumber, initDataTime);

            // at some point in the future levelGhostParticles could be filled with injected
            // particles depending on the domain boundary condition.
            //
            // Do we need J ghosts filled here?
            // This method is only called when root level is initialized
            // but J ghosts are needed a priori for the laplacian when the first Ohm is
            // calculated so I think we do, not having them here is just having the
            // laplacian wrong on L0 borders for the initial E, which is not the end of the
            // world...
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
        // domain node we thus fill ghost nodes. note that we first fill shared border nodes
        // with the sharedNode refiners so that these shared nodes agree on their value at
        // MPI process boundaries. then regular refiner fill are called, which fill only
        // pure ghost nodes. note also that moments are not filled on border nodes since
        // already OK from particle deposition
        void postSynchronize(IPhysicalModel& model, level_t& level, double const time) override
        {
            auto levelNumber  = level.getLevelNumber();
            auto& hybridModel = static_cast<HybridModel&>(model);

            PHARE_LOG_LINE_STR("postSynchronize level " + std::to_string(levelNumber))

            // should we keep the filling on electrif ghosts if done in reflux?
            elecGhostsRefiners_.fill(hybridModel.state.electromag.E, levelNumber, time);
            chargeDensityPatchGhostsRefiners_.fill(levelNumber, time);
            velPatchGhostsRefiners_.fill(hybridModel.state.ions.velocity(), levelNumber, time);
        }

    private:
        void registerGhostComms_(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            // all of the ghost refiners take the nonOverwriteInteriorTFfillPattern as they should
            // only ever modify the ghost and never the interior domain
            elecGhostsRefiners_.addStaticRefiners(info->ghostElectric, EfieldRefineOp_,
                                                  info->ghostElectric,
                                                  nonOverwriteInteriorTFfillPattern);

            // static refinement for J  because it is a temporary, so keeping its
            // state updated after each regrid is not a priority. However if we do not correctly
            // refine on regrid, the post regrid state is not up to date (in our case it will be nan
            // since we nan-initialise) and thus is is better to rely on static refinement, which
            // uses the state after computation of ampere.
            currentGhostsRefiners_.addStaticRefiners(info->ghostCurrent, EfieldRefineOp_,
                                                     info->ghostCurrent,
                                                     nonOverwriteInteriorTFfillPattern);

            chargeDensityGhostsRefiners_.addTimeRefiner(
                info->modelIonDensity, info->modelIonDensity, NiOld_.name(), fieldMomentsRefineOp_,
                fieldTimeOp_, info->modelIonDensity, overwriteInteriorFieldFillPattern);


            velGhostsRefiners_.addTimeRefiners(info->ghostBulkVelocity, info->modelIonBulkVelocity,
                                               ViOld_.name(), vecFieldMomentsRefineOp_,
                                               vecFieldTimeOp_, overwriteInteriorTFfillPattern);

            chargeDensityPatchGhostsRefiners_.addTimeRefiner(
                info->modelIonDensity, info->modelIonDensity, NiOld_.name(), fieldMomentsRefineOp_,
                fieldTimeOp_, info->modelIonDensity, defaultFieldFillPattern);

            velPatchGhostsRefiners_.addTimeRefiners(
                info->ghostBulkVelocity, info->modelIonBulkVelocity, ViOld_.name(),
                vecFieldMomentsRefineOp_, vecFieldTimeOp_, nonOverwriteInteriorTFfillPattern);
        }




        void registerInitComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto b_id = resourcesManager_->getID(info->modelMagnetic);
            BalgoInit.registerRefine(*b_id, *b_id, *b_id, BfieldRefineOp_,
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


            domainGhostPartRefiners_.addStaticRefiners(
                info->patchGhostParticles, nullptr, info->patchGhostParticles,
                std::make_shared<ParticleDomainFromGhostFillPattern<GridLayoutT>>());


            for (auto const& vecfield : info->ghostFlux)
            {
                popFluxBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        sumVec_.name(), vecfield, nullptr, sumVec_.name(),
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
            }

            for (auto const& field : info->sumBorderFields)
                popDensityBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        sumField_.name(), field, nullptr, sumField_.name(),
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());
        }



        void registerSyncComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            electroSynchronizers_.add(info->modelElectric, electricFieldCoarseningOp_,
                                      info->modelElectric);

            ionBulkVelSynchronizers_.add(info->modelIonBulkVelocity, vecFieldCoarseningOp_,
                                         info->modelIonBulkVelocity);

            chargeDensitySynchronizers_.add(info->modelIonDensity, fieldCoarseningOp_,
                                            info->modelIonDensity);
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


        /** * @brief setNaNsFieldOnGhosts sets NaNs on the ghost nodes of the field
         *
         * NaNs are set on all ghost nodes, patch ghost or level ghost nodes
         * so that the refinement operators can know nodes at NaN have not been
         * touched by schedule copy.
         *
         * This is needed when the schedule copy is done before refinement
         * as a result of FieldVariable::fineBoundaryRepresentsVariable=false
         */
        void setNaNsOnFieldGhosts(FieldT& field, patch_t const& patch)
        {
            auto const qty         = field.physicalQuantity();
            using qty_t            = std::decay_t<decltype(qty)>;
            using field_geometry_t = FieldGeometry<GridLayoutT, qty_t>;

            auto const box    = patch.getBox();
            auto const layout = layoutFromPatch<GridLayoutT>(patch);

            // we need to remove the box from the ghost box
            // to use SAMRAI::removeIntersections we do some conversions to
            // samrai box.
            // note gbox is a fieldBox (thanks to the layout)

            auto const gbox  = layout.AMRGhostBoxFor(field.physicalQuantity());
            auto const sgbox = samrai_box_from(gbox);
            auto const fbox  = field_geometry_t::toFieldBox(box, qty, layout);

            // we have field samrai boxes so we can now remove one from the other
            SAMRAI::hier::BoxContainer ghostLayerBoxes{};
            ghostLayerBoxes.removeIntersections(sgbox, fbox);

            // and now finally set the NaNs on the ghost boxes
            for (auto const& gb : ghostLayerBoxes)
                for (auto const& index : layout.AMRToLocal(phare_box_from<dimension>(gb)))
                    field(index) = std::numeric_limits<typename VecFieldT::value_type>::quiet_NaN();
        }

        void setNaNsOnFieldGhosts(FieldT& field, level_t const& level)
        {
            for (auto& patch : resourcesManager_->enumerate(level, field))
                setNaNsOnFieldGhosts(field, *patch);
        }

        void setNaNsOnVecfieldGhosts(VecFieldT& vf, level_t const& level)
        {
            for (auto& patch : resourcesManager_->enumerate(level, vf))
                for (auto& component : vf)
                    setNaNsOnFieldGhosts(component, *patch);
        }


        VecFieldT Jold_{stratName + "_Jold", core::HybridQuantity::Vector::J};
        VecFieldT ViOld_{stratName + "_VBulkOld", core::HybridQuantity::Vector::V};
        FieldT NiOld_{stratName + "_NiOld", core::HybridQuantity::Scalar::rho};

        TensorFieldT sumTensor_{"PHARE_sumTensor", core::HybridQuantity::Tensor::M};
        VecFieldT sumVec_{"PHARE_sumVec", core::HybridQuantity::Vector::V};
        FieldT sumField_{"PHARE_sumField", core::HybridQuantity::Scalar::rho};



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
        using FieldGhostSumRefinerPool    = RefinerPool<rm_t, RefinerType::PatchFieldBorderSum>;
        using VecFieldGhostSumRefinerPool = RefinerPool<rm_t, RefinerType::PatchVecFieldBorderSum>;
        using FieldFillPattern_t          = FieldFillPattern<dimension>;
        using TensorFieldFillPattern_t    = TensorFieldFillPattern<dimension /*, rank=1*/>;

        //! += flux on ghost box overlap incomplete population moment nodes
        std::vector<VecFieldGhostSumRefinerPool> popFluxBorderSumRefiners_;
        //! += density on ghost box overlap incomplete population moment nodes
        std::vector<FieldGhostSumRefinerPool> popDensityBorderSumRefiners_;

        InitRefinerPool electricInitRefiners_{resourcesManager_};


        SAMRAI::xfer::RefineAlgorithm BalgoPatchGhost;
        SAMRAI::xfer::RefineAlgorithm BghostAlgo;
        SAMRAI::xfer::RefineAlgorithm BPredGhostAlgo;
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

        //! store refiners for electric fields that need ghosts to be filled
        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};

        GhostRefinerPool currentGhostsRefiners_{resourcesManager_};

        // moment ghosts
        // The border node is already complete by the deposit of ghost particles
        // these refiners are used to fill ghost nodes, and therefore, owing to
        // the GhostField tag, will only assign pure ghost nodes. Border nodes will
        // be overwritten only on level borders, which does not seem to be an issue.
        LevelBorderFieldRefinerPool chargeDensityGhostsRefiners_{resourcesManager_};
        LevelBorderFieldRefinerPool velGhostsRefiners_{resourcesManager_};

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
        DomainGhostPartRefinerPool domainGhostPartRefiners_{resourcesManager_};

        SynchronizerPool<rm_t> chargeDensitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> ionBulkVelSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> electroSynchronizers_{resourcesManager_};


        RefOp_ptr fieldRefineOp_{std::make_shared<DefaultFieldRefineOp>()};
        RefOp_ptr vecFieldRefineOp_{std::make_shared<DefaultVecFieldRefineOp>()};

        RefOp_ptr fieldMomentsRefineOp_{std::make_shared<FieldMomentsRefineOp>()};
        RefOp_ptr vecFieldMomentsRefineOp_{std::make_shared<VecFieldMomentsRefineOp>()};

        RefOp_ptr BfieldRefineOp_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr BfieldRegridOp_{std::make_shared<MagneticFieldRegridOp>()};
        RefOp_ptr EfieldRefineOp_{std::make_shared<ElectricFieldRefineOp>()};
        std::shared_ptr<FieldFillPattern_t> defaultFieldFillPattern
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
        CoarsenOperator_ptr fieldCoarseningOp_{std::make_shared<DefaultFieldCoarsenOp>()};
        CoarsenOperator_ptr vecFieldCoarseningOp_{std::make_shared<DefaultVecFieldCoarsenOp>()};
        CoarsenOperator_ptr electricFieldCoarseningOp_{std::make_shared<ElectricFieldCoarsenOp>()};

        MagneticRefinePatchStrategy<ResourcesManagerT, VectorFieldDataT>
            magneticRefinePatchStrategy_{*resourcesManager_};

        // MagneticRefinePatchStrategy<ResourcesManagerT, VectorFieldDataT>
        // BpredRefinePatchStrategy_{
        //     *resourcesManager_};
    };


} // namespace amr

} // namespace PHARE

#endif
