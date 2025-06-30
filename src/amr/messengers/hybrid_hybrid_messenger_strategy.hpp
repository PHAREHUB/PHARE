#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP

#include "core/def.hpp"
#include "core/debug.hpp"
#include "core/logger.hpp"
#include "core/def/phare_mpi.hpp"

#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/IntVector.h"
#include "core/utilities/index/index.hpp"
#include "refiner_pool.hpp"
#include "synchronizer_pool.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/field/coarsening/magnetic_field_coarsener.hpp"
#include "amr/data/field/refine/field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_refiner.hpp"
#include "amr/data/field/refine/electric_field_refiner.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/data/vecfield/vecfield.hpp"
#include "core/utilities/point/point.hpp"



#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/xfer/BoxGeometryVariableFillPattern.h"


#include <iterator>
#include <optional>
#include <utility>
#include <iomanip>
#include <iostream>
#include <string>


namespace PHARE
{
namespace amr
{
    // when registering different components to the same algorithm in SAMRAI, as we want to do for
    // vecfields, we need those components not to be considered as equivalent_classes by SAMRAI.
    // Without this precaution SAMRAI will assume the same geometry for all.
    class XVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
    {
    };

    class YVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
    {
    };

    class ZVariableFillPattern : public SAMRAI::xfer::BoxGeometryVariableFillPattern
    {
    };

    /** \brief An HybridMessenger is the specialization of a HybridMessengerStrategy for hybrid
     * to hybrid data communications.
     */
    template<typename HybridModel, typename RefinementParams>
    class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
    {
        using GridT             = typename HybridModel::grid_type;
        using IonsT             = typename HybridModel::ions_type;
        using ElectromagT       = typename HybridModel::electromag_type;
        using VecFieldT         = typename HybridModel::vecfield_type;
        using GridLayoutT       = typename HybridModel::gridlayout_type;
        using FieldT            = typename VecFieldT::field_type;
        using FieldDataT        = FieldData<GridLayoutT, GridT>;
        using ResourcesManagerT = typename HybridModel::resources_manager_type;
        using IPhysicalModel    = typename HybridModel::Interface;

        static constexpr std::size_t dimension   = GridLayoutT::dimension;
        static constexpr std::size_t interpOrder = GridLayoutT::interp_order;

        using InteriorParticleRefineOp = typename RefinementParams::InteriorParticleRefineOp;
        using CoarseToFineRefineOpOld  = typename RefinementParams::CoarseToFineRefineOpOld;
        using CoarseToFineRefineOpNew  = typename RefinementParams::CoarseToFineRefineOpNew;

        template<typename Policy>
        using BaseRefineOp          = FieldRefineOperator<GridLayoutT, GridT, Policy>;
        using DefaultFieldRefineOp  = BaseRefineOp<DefaultFieldRefiner<dimension>>;
        using MagneticFieldRefineOp = BaseRefineOp<MagneticFieldRefiner<dimension>>;
        using ElectricFieldRefineOp = BaseRefineOp<ElectricFieldRefiner<dimension>>;
        using FieldTimeInterp       = FieldLinearTimeInterpolate<GridLayoutT, GridT>;

        template<typename Policy>
        using BaseCoarsenOp     = FieldCoarsenOperator<GridLayoutT, GridT, Policy>;
        using MagneticCoarsenOp = BaseCoarsenOp<MagneticFieldCoarsener<dimension>>;
        using DefaultCoarsenOp  = BaseCoarsenOp<DefaultFieldCoarsener<dimension>>;

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
        }

        virtual ~HybridHybridMessengerStrategy() = default;



        /* ------------------------------------------------------------------------
                        methods used for the IMessenger interface
           ------------------------------------------------------------------------ */


        /**
         * @brief allocate the messenger strategy internal variables to the model
         * resourceManager
         */
        void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            resourcesManager_->allocate(Jold_, patch, allocateTime);
            resourcesManager_->allocate(NiOld_, patch, allocateTime);
            resourcesManager_->allocate(ViOld_, patch, allocateTime);
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


            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> xVariableFillPattern
                = std::make_shared<XVariableFillPattern>();

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> yVariableFillPattern
                = std::make_shared<YVariableFillPattern>();

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> zVariableFillPattern
                = std::make_shared<ZVariableFillPattern>();

            auto bx_id = resourcesManager_->getID(hybridInfo->modelMagnetic.xName);
            auto by_id = resourcesManager_->getID(hybridInfo->modelMagnetic.yName);
            auto bz_id = resourcesManager_->getID(hybridInfo->modelMagnetic.zName);

            if (!bx_id or !by_id or !bz_id)
            {
                throw std::runtime_error(
                    "HybridHybridMessengerStrategy: missing magnetic field variable IDs");
            }

            magneticRefinePatchStrategy_.registerIDs(*bx_id, *by_id, *bz_id);

            Balgo.registerRefine(*bx_id, *bx_id, *bx_id, BfieldRefineOp_, xVariableFillPattern);
            Balgo.registerRefine(*by_id, *by_id, *by_id, BfieldRefineOp_, yVariableFillPattern);
            Balgo.registerRefine(*bz_id, *bz_id, *bz_id, BfieldRefineOp_, zVariableFillPattern);

            BalgoNode.registerRefine(*bx_id, *bx_id, *bx_id, BfieldNodeRefineOp_,
                                     xVariableFillPattern);
            BalgoNode.registerRefine(*by_id, *by_id, *by_id, BfieldNodeRefineOp_,
                                     yVariableFillPattern);
            BalgoNode.registerRefine(*bz_id, *bz_id, *bz_id, BfieldNodeRefineOp_,
                                     zVariableFillPattern);

            registerGhostComms_(hybridInfo);
            registerInitComms(hybridInfo);
            registerSyncComms(hybridInfo);
        }



        /**
         * @brief all RefinerPool must be notified the level levelNumber now exist.
         * not doing so will result in communication to/from that level being impossible
         */
        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magSharedNodeRefineSchedules[levelNumber]
                = BalgoNode.createSchedule(level, &magneticRefinePatchStrategy_);

            magPatchGhostsRefineSchedules[levelNumber]
                = Balgo.createSchedule(level, &magneticRefinePatchStrategy_);

            magGhostsRefineSchedules[levelNumber] = Balgo.createSchedule(
                level, levelNumber - 1, hierarchy, &magneticRefinePatchStrategy_);

            elecSharedNodesRefiners_.registerLevel(hierarchy, level);
            currentSharedNodesRefiners_.registerLevel(hierarchy, level);

            elecGhostsRefiners_.registerLevel(hierarchy, level);
            currentGhostsRefiners_.registerLevel(hierarchy, level);

            rhoGhostsRefiners_.registerLevel(hierarchy, level);
            velGhostsRefiners_.registerLevel(hierarchy, level);

            patchGhostPartRefiners_.registerLevel(hierarchy, level);


            // root level is not initialized with a schedule using coarser level data
            // so we don't create these schedules if root level
            // TODO this 'if' may not be OK if L0 is regrided
            if (levelNumber != rootLevelNumber)
            {
                // those are for refinement
                magInitRefineSchedules[levelNumber] = Balgo.createSchedule(
                    level, nullptr, levelNumber - 1, hierarchy, &magneticRefinePatchStrategy_);
                electricInitRefiners_.registerLevel(hierarchy, level);
                domainParticlesRefiners_.registerLevel(hierarchy, level);
                lvlGhostPartOldRefiners_.registerLevel(hierarchy, level);
                lvlGhostPartNewRefiners_.registerLevel(hierarchy, level);

                // and these for coarsening
                magnetoSynchronizers_.registerLevel(hierarchy, level);
                electroSynchronizers_.registerLevel(hierarchy, level);
                densitySynchronizers_.registerLevel(hierarchy, level);
                ionBulkVelSynchronizers_.registerLevel(hierarchy, level);
            }
        }



        /**
         * @brief regrid performs the regriding communications for Hybrid to Hybrid messengers
         , all quantities that are in initialization refiners need to be regridded
         */
        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                    int const levelNumber,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                    IPhysicalModel& model, double const initDataTime) override
        {
            auto& hybridModel = dynamic_cast<HybridModel&>(model);
            auto level        = hierarchy->getPatchLevel(levelNumber);

            bool isRegriddingL0 = levelNumber == 0 and oldLevel;

            magneticRegriding_(hierarchy, level, oldLevel, hybridModel, initDataTime);
            electricInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            domainParticlesRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            patchGhostPartRefiners_.fill(levelNumber, initDataTime);


            // regriding will fill the new level wherever it has points that overlap
            // old level. This will include its level border points.
            // These new level border points will thus take values that where previous
            // domain values. Magnetic flux is thus not necessarily consistent with
            // the Loring et al. method to sync the induction between coarse and fine faces.
            // Specifically, we need all fine faces to have equal magnetic field and also
            // equal to that of the shared coarse face.
            // This means that we now need to fill ghosts and border included

            if (!isRegriddingL0)
            {
                auto& E = hybridModel.state.electromag.E;
                elecGhostsRefiners_.fill(E, levelNumber, initDataTime);
            }

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

            magPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);

            PHARE_DEBUG_SCOPE("HyHyMessStrat/regrid/");
            PHARE_DEBUG_CHECK_LEVEL(GridLayoutT, *resourcesManager_, *level);
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
        void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                       double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();


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
            // however we need to call the ghost communicator for patch ghost particles
            // since the interior schedules have a restriction to the interior of the patch.
            PHARE_LOG_START(3, "hybhybmessengerStrat::initLevel : patch ghost part fill schedule");
            patchGhostPartRefiners_.fill(levelNumber, initDataTime);
            PHARE_LOG_STOP(3, "hybhybmessengerStrat::initLevel : patch ghost part fill schedule");


            lvlGhostPartOldRefiners_.fill(levelNumber, initDataTime);


            // levelGhostParticles will be pushed during the advance phase
            // they need to be identical to levelGhostParticlesOld before advance
            copyLevelGhostOldToPushable_(level, model);

            // computeIonMoments_(level, model);
        }



        /* ------------------------------------------------------------------------
                     methods used for the HybridMessenger interface
           ------------------------------------------------------------------------ */



        void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillElectricGhosts");
            elecSharedNodesRefiners_.fill(E, levelNumber, fillTime);
            elecGhostsRefiners_.fill(E, levelNumber, fillTime);
        }




        void fillCurrentGhosts(VecFieldT& J, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillCurrentGhosts");
            currentSharedNodesRefiners_.fill(J, levelNumber, fillTime);
            currentGhostsRefiners_.fill(J, levelNumber, fillTime);
        }




        /**
         * @brief fillIonGhostParticles will fill the interior ghost particle array from
         * neighbor patches of the same level. Before doing that, it empties the array for
         * all populations
         */
        void fillIonGhostParticles(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                   double const fillTime) override
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillIonGhostParticles");

            for (auto patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                for (auto& pop : ions)
                {
                    pop.patchGhostParticles().clear();
                }
            }
            patchGhostPartRefiners_.fill(level.getLevelNumber(), fillTime);
        }




        /**
         * @brief fillIonPopMomentGhosts works on moment ghost nodes
         *
         * patch border node moments are completed by the deposition of patch ghost
         * particles for all populations level border nodes are completed by the deposition
         * of level ghost [old,new] particles for all populations, linear time interpolation
         * is used to get the contribution of old/new particles
         */
        void fillIonPopMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                    double const afterPushTime) override
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::fillIonMomentGhosts");

            auto alpha = timeInterpCoef_(afterPushTime, level.getLevelNumber());
            if (level.getLevelNumber() > 0 and (alpha < 0 or alpha > 1))
            {
                std::cout << std::setprecision(12) << alpha << "\n";
                throw std::runtime_error("ion moment ghost time interp coef invalid : alpha: "
                                         + std::to_string(alpha) + " beforePushTime "
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

                    interpolate_(makeRange(patchGhosts), density, flux, layout);

                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // then grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        auto& levelGhostOld = pop.levelGhostParticlesOld();
                        interpolate_(makeRange(levelGhostOld), density, flux, layout, 1. - alpha);

                        auto& levelGhostNew = pop.levelGhostParticlesNew();
                        interpolate_(makeRange(levelGhostNew), density, flux, layout, alpha);
                    }
                }
            }
        }


        /* pure (patch and level) ghost nodes are filled by applying a regular ghost
         * schedule i.e. that does not overwrite the border patch node previously well
         * calculated from particles Note : the ghost schedule only fills the total density
         * and bulk velocity and NOT population densities and fluxes. These partial
         * densities and fluxes are thus not available on ANY ghost node.*/
        virtual void fillIonMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                         double const afterPushTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillIonMomentGhosts");
            rhoGhostsRefiners_.fill(level.getLevelNumber(), afterPushTime);
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
        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& level,
                       std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                       double const currentTime, double const prevCoarserTime,
                       double const newCoarserTime) override
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
        void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) override
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
        void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                         double currentTime) override
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
                auto& Ni = hybridModel.state.ions.density();

                Jold_.copyData(J);
                ViOld_.copyData(Vi);
                NiOld_.copyData(Ni);
            }
        }




        void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                            double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();
            assert(levelNumber == 0);

            auto& hybridModel = static_cast<HybridModel&>(model);

            elecSharedNodesRefiners_.fill(hybridModel.state.electromag.E, levelNumber,
                                          initDataTime);

            elecGhostsRefiners_.fill(hybridModel.state.electromag.E, levelNumber, initDataTime);
            patchGhostPartRefiners_.fill(levelNumber, initDataTime);

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



        void synchronize(SAMRAI::hier::PatchLevel& level) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::synchronize");

            auto levelNumber = level.getLevelNumber();
            PHARE_LOG_LINE_STR("synchronizing level " + std::to_string(levelNumber));

            // call coarsning schedules...
            magnetoSynchronizers_.sync(levelNumber);
            electroSynchronizers_.sync(levelNumber);
            densitySynchronizers_.sync(levelNumber);
            ionBulkVelSynchronizers_.sync(levelNumber);

            PHARE_DEBUG_SCOPE("HyHyMessStrat/after_sync/");
            PHARE_DEBUG_CHECK_LEVEL(GridLayoutT, *resourcesManager_, level);
        }

        // after coarsening, domain nodes have been updated and therefore patch ghost nodes
        // will probably stop having the exact same value as their overlapped neighbor
        // domain node we thus fill ghost nodes. note that we first fill shared border nodes
        // with the sharedNode refiners so that these shared nodes agree on their value at
        // MPI process boundaries. then regular refiner fill are called, which fill only
        // pure ghost nodes. note also that moments are not filled on border nodes since
        // already OK from particle deposition
        void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                             double const time) override
        {
            auto levelNumber  = level.getLevelNumber();
            auto& hybridModel = static_cast<HybridModel&>(model);

            PHARE_LOG_LINE_STR("postSynchronize level " + std::to_string(levelNumber))

            magSharedNodeRefineSchedules[levelNumber]->fillData(time);
            elecSharedNodesRefiners_.fill(hybridModel.state.electromag.E, levelNumber, time);

            // we fill magnetic field ghosts only on patch ghost nodes and not on level
            // ghosts the reason is that 1/ filling ghosts is necessary to prevent mismatch
            // between ghost and overlaped neighboring patch domain nodes resulting from
            // former coarsening which does not occur for level ghosts and 2/ overwriting
            // level border with next coarser model B would invalidate divB on the first
            // fine domain cell since its border face only received a fraction of the
            // induction that has occured on the shared coarse face.
            magPatchGhostsRefineSchedules[levelNumber]->fillData(time);
            elecGhostsRefiners_.fill(hybridModel.state.electromag.E, levelNumber, time);
            rhoGhostsRefiners_.fill(levelNumber, time);
            velGhostsRefiners_.fill(hybridModel.state.ions.velocity(), levelNumber, time);

            PHARE_DEBUG_SCOPE("HyHyMessStrat/post_sync/after/");
            PHARE_DEBUG_CHECK_LEVEL(GridLayoutT, *resourcesManager_, level);
        }

    private:
        void registerGhostComms_(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& vecFieldNames) {
                std::vector<std::string> keys;
                std::transform(std::begin(vecFieldNames), std::end(vecFieldNames),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            elecSharedNodesRefiners_.addStaticRefiners(info->ghostElectric, EfieldNodeRefineOp_,
                                                       makeKeys(info->ghostElectric));

            elecGhostsRefiners_.addStaticRefiners(info->ghostElectric, EfieldRefineOp_,
                                                  makeKeys(info->ghostElectric));

            currentSharedNodesRefiners_.addTimeRefiners(info->ghostCurrent, info->modelCurrent,
                                                        core::VecFieldNames{Jold_},
                                                        EfieldNodeRefineOp_, fieldTimeOp_);

            currentGhostsRefiners_.addTimeRefiners(info->ghostCurrent, info->modelCurrent,
                                                   core::VecFieldNames{Jold_}, EfieldRefineOp_,
                                                   fieldTimeOp_);

            rhoGhostsRefiners_.addTimeRefiner(info->modelIonDensity, info->modelIonDensity,
                                              NiOld_.name(), fieldRefineOp_, fieldTimeOp_,
                                              info->modelIonDensity);


            velGhostsRefiners_.addTimeRefiners(info->ghostBulkVelocity, info->modelIonBulkVelocity,
                                               core::VecFieldNames{ViOld_}, fieldRefineOp_,
                                               fieldTimeOp_);
        }




        void registerInitComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& descriptor) {
                std::vector<std::string> keys;
                std::transform(std::begin(descriptor), std::end(descriptor),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            electricInitRefiners_.addStaticRefiners(info->initElectric, EfieldRefineOp_,
                                                    makeKeys(info->initElectric));


            domainParticlesRefiners_.addStaticRefiners(
                info->interiorParticles, interiorParticleRefineOp_, info->interiorParticles);


            lvlGhostPartOldRefiners_.addStaticRefiners(info->levelGhostParticlesOld,
                                                       levelGhostParticlesOldOp_,
                                                       info->levelGhostParticlesOld);


            lvlGhostPartNewRefiners_.addStaticRefiners(info->levelGhostParticlesNew,
                                                       levelGhostParticlesNewOp_,
                                                       info->levelGhostParticlesNew);


            patchGhostPartRefiners_.addStaticRefiners(info->patchGhostParticles, nullptr,
                                                      info->patchGhostParticles);
        }




        void registerSyncComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            magnetoSynchronizers_.add(info->modelMagnetic, magneticCoarseningOp_,
                                      info->modelMagnetic.vecName);

            electroSynchronizers_.add(info->modelElectric, fieldCoarseningOp_,
                                      info->modelElectric.vecName);

            ionBulkVelSynchronizers_.add(info->modelIonBulkVelocity, fieldCoarseningOp_,
                                         info->modelIonBulkVelocity.vecName);

            densitySynchronizers_.add(info->modelIonDensity, fieldCoarseningOp_,
                                      info->modelIonDensity);
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

                    levelGhostParticles = levelGhostParticlesOld;
                }
            }
        }




        double timeInterpCoef_(double const afterPushTime, std::size_t levelNumber)
        {
            return (afterPushTime - beforePushCoarseTime_[levelNumber])
                   / (afterPushCoarseTime_[levelNumber] - beforePushCoarseTime_[levelNumber]);
        }




        void magneticRegriding_(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                                std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                                HybridModel& hybridModel, double const initDataTime)
        {
            // first we set all B ghost nodes to NaN so that we can later
            // postprocess them and fill them with the correct value
            for (auto& patch : *level)
            {
                auto const& layout = layoutFromPatch<GridLayoutT>(*patch);
                auto _  = resourcesManager_->setOnPatch(*patch, hybridModel.state.electromag.B);
                auto& B = hybridModel.state.electromag.B;

                auto setToNaN = [&](auto& B, core::MeshIndex<dimension> idx) {
                    B(idx) = std::numeric_limits<double>::quiet_NaN();
                };

                layout.evalOnGhostBox(B(core::Component::X), [&](auto&... args) mutable {
                    setToNaN(B(core::Component::X), {args...});
                });
                layout.evalOnGhostBox(B(core::Component::Y), [&](auto&... args) mutable {
                    setToNaN(B(core::Component::Y), {args...});
                });
                layout.evalOnGhostBox(B(core::Component::Z), [&](auto&... args) mutable {
                    setToNaN(B(core::Component::Z), {args...});
                });
            }

            // here we create the schedule on the fly because it is the only moment where we
            // have both the old and current level

            auto magSchedule = Balgo.createSchedule(
                level, oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy);
            magSchedule->fillData(initDataTime);

            // we set the new fine faces using the toth and roe (2002) formulas. This requires
            // an even number of ghost cells as we set the new fine faces using the values of
            // the fine faces shared with the corresponding coarse faces of the coarse cell.
            for (auto& patch : *level)
            {
                auto const& layout = layoutFromPatch<GridLayoutT>(*patch);
                auto _   = resourcesManager_->setOnPatch(*patch, hybridModel.state.electromag.B);
                auto& B  = hybridModel.state.electromag.B;
                auto& bx = B(core::Component::X);
                auto& by = B(core::Component::Y);
                auto& bz = B(core::Component::Z);

                if constexpr (dimension == 1)
                {
                    auto postprocessBx = [&](core::MeshIndex<dimension> idx) {
                        auto ix = idx[dirX];

                        if (std::isnan(bx(ix)))
                        {
                            assert(ix % 2 == 1);
                            MagneticRefinePatchStrategy<ResourcesManagerT,
                                                        FieldDataT>::postprocessBx1d(bx, idx);
                        }
                    };

                    layout.evalOnGhostBox(B(core::Component::X),
                                          [&](auto&... args) mutable { postprocessBx({args...}); });
                }
                else if constexpr (dimension == 2)
                {
                    auto postprocessBx = [&](core::MeshIndex<dimension> idx) {
                        auto ix = idx[dirX];
                        auto iy = idx[dirY];

                        if (std::isnan(bx(ix, iy)))
                        {
                            assert(ix % 2 == 1);
                            MagneticRefinePatchStrategy<ResourcesManagerT,
                                                        FieldDataT>::postprocessBx2d(bx, by, idx);
                        }
                    };

                    auto postprocessBy = [&](core::MeshIndex<dimension> idx) {
                        auto ix = idx[dirX];
                        auto iy = idx[dirY];

                        if (std::isnan(by(ix, iy)))
                        {
                            assert(iy % 2 == 1);
                            MagneticRefinePatchStrategy<ResourcesManagerT,
                                                        FieldDataT>::postprocessBy2d(bx, by, idx);
                        }
                    };

                    layout.evalOnGhostBox(B(core::Component::X),
                                          [&](auto&... args) mutable { postprocessBx({args...}); });

                    layout.evalOnGhostBox(B(core::Component::Y),
                                          [&](auto&... args) mutable { postprocessBy({args...}); });
                }
                else if constexpr (dimension == 3)
                {
                    auto meshSize = layout.meshSize();

                    auto postprocessBx = [&](core::MeshIndex<dimension> idx) {
                        auto ix = idx[dirX];
                        auto iy = idx[dirY];
                        auto iz = idx[dirZ];

                        if (std::isnan(bx(ix, iy, iz)))
                        {
                            assert(ix % 2 == 1);
                            MagneticRefinePatchStrategy<ResourcesManagerT,
                                                        FieldDataT>::postprocessBx3d(bx, by, bz,
                                                                                     meshSize, idx);
                        }
                    };

                    auto postprocessBy = [&](core::MeshIndex<dimension> idx) {
                        auto ix = idx[dirX];
                        auto iy = idx[dirY];
                        auto iz = idx[dirZ];

                        if (std::isnan(by(ix, iy, iz)))
                        {
                            assert(iy % 2 == 1);
                            MagneticRefinePatchStrategy<ResourcesManagerT,
                                                        FieldDataT>::postprocessBy3d(bx, by, bz,
                                                                                     meshSize, idx);
                        }
                    };

                    auto postprocessBz = [&](core::MeshIndex<dimension> idx) {
                        auto ix = idx[dirX];
                        auto iy = idx[dirY];
                        auto iz = idx[dirZ];

                        if (std::isnan(bz(ix, iy, iz)))
                        {
                            assert(iz % 2 == 1);
                            MagneticRefinePatchStrategy<ResourcesManagerT,
                                                        FieldDataT>::postprocessBz3d(bx, by, bz,
                                                                                     meshSize, idx);
                        }
                    };

                    layout.evalOnGhostBox(B(core::Component::X),
                                          [&](auto&... args) mutable { postprocessBx({args...}); });

                    layout.evalOnGhostBox(B(core::Component::Y),
                                          [&](auto&... args) mutable { postprocessBy({args...}); });

                    layout.evalOnGhostBox(B(core::Component::Z),
                                          [&](auto&... args) mutable { postprocessBz({args...}); });
                }

                auto notNan = [&](auto& b, core::MeshIndex<dimension> idx) {
                    auto check = [&](auto&&... indices) {
                        if (std::isnan(b(indices...)))
                        {
                            std::string index_str;
                            ((index_str
                              += (index_str.empty() ? "" : ", ") + std::to_string(indices)),
                             ...);
                            throw std::runtime_error("NaN found in magnetic field " + b.name()
                                                     + " at index (" + index_str + ")");
                        }
                    };

                    if constexpr (dimension == 1)
                    {
                        check(idx[dirX]);
                    }
                    else if constexpr (dimension == 2)
                    {
                        check(idx[dirX], idx[dirY]);
                    }
                    else if constexpr (dimension == 3)
                    {
                        check(idx[dirX], idx[dirY], idx[dirZ]);
                    }
                };

                auto checkNoNaNsLeft = [&]() {
                    auto checkComponent = [&](auto component) {
                        layout.evalOnGhostBox(
                            B(component), [&](auto&... args) { notNan(B(component), {args...}); });
                    };

                    checkComponent(core::Component::X);
                    checkComponent(core::Component::Y);
                    checkComponent(core::Component::Z);
                };

                PHARE_DEBUG_DO(checkNoNaNsLeft());
            }
        }




        VecFieldT Jold_{stratName + "_Jold", core::HybridQuantity::Vector::J};
        VecFieldT ViOld_{stratName + "_VBulkOld", core::HybridQuantity::Vector::V};
        FieldT NiOld_{stratName + "_NiOld", core::HybridQuantity::Scalar::rho};


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
        using InitRefinerPool           = RefinerPool<rm_t, RefinerType::InitField>;
        using SharedNodeRefinerPool     = RefinerPool<rm_t, RefinerType::SharedBorder>;
        using GhostRefinerPool          = RefinerPool<rm_t, RefinerType::GhostField>;
        using PatchGhostRefinerPool     = RefinerPool<rm_t, RefinerType::PatchGhostField>;
        using InitDomPartRefinerPool    = RefinerPool<rm_t, RefinerType::InitInteriorPart>;
        using PatchGhostPartRefinerPool = RefinerPool<rm_t, RefinerType::InteriorGhostParticles>;

        InitRefinerPool electricInitRefiners_{resourcesManager_};

        SAMRAI::xfer::RefineAlgorithm Balgo;
        SAMRAI::xfer::RefineAlgorithm BalgoNode;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magInitRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magSharedNodeRefineSchedules;


        //! store refiners for electric fields that need ghosts to be filled
        SharedNodeRefinerPool elecSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};

        GhostRefinerPool currentSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool currentGhostsRefiners_{resourcesManager_};

        // moment ghosts
        // these do not need sharedNode refiners. The reason is that
        // the border node is already complete by the deposit of ghost particles
        // these refiners are used to fill ghost nodes, and therefore, owing to
        // the GhostField tag, will only assign pur ghost nodes. Border nodes will
        // be overwritten only on level borders, which does not seem to be an issue.
        GhostRefinerPool rhoGhostsRefiners_{resourcesManager_};
        GhostRefinerPool velGhostsRefiners_{resourcesManager_};

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


        // this contains refiners for each population to exchange patch ghost particles
        PatchGhostPartRefinerPool patchGhostPartRefiners_{resourcesManager_};

        SynchronizerPool<rm_t> densitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> ionBulkVelSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> electroSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> magnetoSynchronizers_{resourcesManager_};


        RefOp_ptr fieldRefineOp_{std::make_shared<DefaultFieldRefineOp>()};
        // see field_variable_fill_pattern.hpp for explanation about this "node_only" flag
        // note that refinement operator, via the boolean argument, serve as a relay for the
        // the refinealgorithm to get the correct variablefillpattern
        RefOp_ptr BfieldNodeRefineOp_{std::make_shared<MagneticFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr BfieldRefineOp_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr EfieldNodeRefineOp_{std::make_shared<ElectricFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr EfieldRefineOp_{std::make_shared<ElectricFieldRefineOp>()};

        std::shared_ptr<TimeInterpolateOperator> fieldTimeOp_{std::make_shared<FieldTimeInterp>()};

        using CoarsenOperator_ptr = std::shared_ptr<SAMRAI::hier::CoarsenOperator>;
        CoarsenOperator_ptr fieldCoarseningOp_{std::make_shared<DefaultCoarsenOp>()};
        CoarsenOperator_ptr magneticCoarseningOp_{std::make_shared<MagneticCoarsenOp>()};

        MagneticRefinePatchStrategy<ResourcesManagerT, FieldDataT> magneticRefinePatchStrategy_{
            *resourcesManager_};
    };


} // namespace amr

} // namespace PHARE

#endif
