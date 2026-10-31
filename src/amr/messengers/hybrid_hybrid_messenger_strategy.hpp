#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP


#include "core/def.hpp" // IWYU pragma: keep
#include "core/logger.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/utilities/types.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/numerics/interpolator/interpolating.hpp"

#include "amr/data/field/field_data.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/messengers/field_operate_transaction.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"
#include "amr/data/field/coarsening/moments_coarsener.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_composite_refiner.hpp"
#include "amr/messengers/refinement_config.hpp"
#include "amr/physical_models/models.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/particles/particles_variable_fill_pattern.hpp"
#include "amr/data/field/coarsening/electric_field_coarsener.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"

#include "refiner_pool.hpp"
#include "synchronizer_pool.hpp"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/hier/CoarseFineBoundary.h>

#include <memory>
#include <string>
#include <utility>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace PHARE
{
namespace amr
{
    /** \brief An HybridMessenger is the specialization of a HybridMessengerStrategy for hybrid
     * to hybrid data communications.
     */
    // emits the moment border fills in another translation unit, see
    // hybrid_hybrid_messenger_strategy_borders.hpp
    template<typename Strategy>
    struct HybridMessengerStrategyBorders;


    template<typename HybridModel, typename RefinementParams>
    class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
    {
        using amr_types     = HybridModel::amr_types;
        using level_t       = amr_types::level_t;
        using patch_t       = amr_types::patch_t;
        using hierarchy_t   = amr_types::hierarchy_t;
        using level_ptr     = amr_types::level_ptr;
        using hierarchy_ptr = amr_types::hierarchy_ptr;

        using GridT       = HybridModel::grid_type;
        using IonsT       = HybridModel::ions_type;
        using ElectromagT = HybridModel::electromag_type;
        using GridLayoutT = HybridModel::gridlayout_type;

        using FieldT     = HybridModel::vecfield_type::field_type;
        using FieldDataT = FieldData<GridLayoutT, GridT>;

        using VecFieldT        = HybridModel::vecfield_type;
        using VectorFieldDataT = TensorFieldData<1, GridLayoutT, GridT, core::HybridQuantity>;

        using TensorFieldT     = IonsT::tensorfield_type;
        using TensorFieldDataT = TensorFieldData<2, GridLayoutT, GridT, core::HybridQuantity>;

        using ResourcesManagerT = HybridModel::resources_manager_type;
        using BoundaryManagerT  = HybridModel::boundary_manager_type;
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

        // if FieldT is tiled, Field_t is not! always not tiled!
        using Grid_t        = HybridModel::storage_t::grid_type;
        using Field_t       = HybridModel::storage_t::field_type;
        using VecField_t    = HybridModel::storage_t::vecfield_type;
        using TensorField_t = HybridModel::storage_t::tensorfield_type;

        using FieldData_t       = FieldData<GridLayoutT, Grid_t>;
        using VectorFieldData_t = TensorFieldData<1, GridLayoutT, Grid_t, core::HybridQuantity>;
        using TensorFieldData_t = TensorFieldData<2, GridLayoutT, Grid_t, core::HybridQuantity>;

        using SetMaxOp     = core::SetMax<typename FieldData_t::value_type>;
        using PlusEqualsOp = core::PlusEquals<typename FieldData_t::value_type>;
        using BorderSumOp  = core::FieldBorderSumOp<typename FieldData_t::value_type>;
        using BorderMaxOp  = core::FieldBorderMaxOp<typename FieldData_t::value_type>;

        using FieldBorderSumTransactionFactory_t
            = FieldBorderOpTransactionFactory<FieldDataT, BorderSumOp, FieldData_t>;
        using VecFieldBorderSumTransactionFactory_t
            = FieldBorderOpTransactionFactory<VectorFieldDataT, BorderSumOp, VectorFieldData_t>;
        using TensorFieldBorderSumTransactionFactory_t
            = FieldBorderOpTransactionFactory<TensorFieldDataT, BorderSumOp, TensorFieldData_t>;

        using FieldBorderMaxTransactionFactory_t
            = FieldBorderOpTransactionFactory<FieldDataT, BorderMaxOp, FieldData_t>;
        using VecFieldBorderMaxTransactionFactory_t
            = FieldBorderOpTransactionFactory<VectorFieldDataT, BorderMaxOp, VectorFieldData_t>;

        static_assert(core::is_field_v<Grid_t> and !core::is_field_tile_set_v<Grid_t>);

    public:
        static inline std::string const stratName    = "HybridModel-HybridModel";
        static constexpr std::size_t rootLevelNumber = 0;


        HybridHybridMessengerStrategy(HybridModel& model, int const firstLevel,
                                      RefinementConfig const& refinementConfig = {})
            : HybridMessengerStrategy<HybridModel>{stratName}
            , model_{model}
            , resourcesManager_{model.resourcesManager}
            , boundaryManager_{model.boundaryManager}
            , firstLevel_{firstLevel}
        {
            makeRefineOperators_(refinementConfig);

            resourcesManager_->registerResources(Bold_);
            resourcesManager_->registerResources(Eold_);
            // scratch to sum the ion moments on patch borders, see fill(Flux|Density)Borders
            model_.requestTemporaries(core::QtyCentering::primal, core::HybridQuantity::Vector::V);
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
            magneticRefinePatchStrategy_.registerScratchID(scratchID_(scratch_vecfield()));

            registerRefineScratchID(scratchID_(scratch_field()), fieldRefineOp_, BInitRefineOp_,
                                    BRefineOp_, EfieldRefineOp_);

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



        /**
         * @brief all RefinerPool must be notified the level levelNumber now exist.
         * not doing so will result in communication to/from that level being impossible
         */
        void registerLevel(hierarchy_ptr const& hierarchy, int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magPatchGhostsRefineSchedules[levelNumber]
                = BalgoPatchGhost.createSchedule(level, &magneticRefinePatchStrategy_);

            elecPatchGhostsRefineSchedules[levelNumber] = EalgoPatchGhost.createSchedule(level);

            // technically not needed for finest as refluxing is not done onto it.
            patchGhostRefluxedSchedules[levelNumber] = PatchGhostRefluxedAlgo.createSchedule(level);

            elecGhostsRefiners_.registerLevel(hierarchy, level);
            magGhostsRefiners_.registerLevel(hierarchy, level);
            domainGhostPartRefiners_.registerLevel(hierarchy, level);

            chargeDensityPatchGhostsRefiners_.registerLevel(hierarchy, level);

            {
                auto border_sum_factory = std::make_shared<VecFieldBorderSumTransactionFactory_t>();
                for (auto& refiner : popFluxBorderSumRefiners_)
                    refiner.registerLevel(hierarchy, level, border_sum_factory);
            }

            {
                auto border_sum_factory = std::make_shared<FieldBorderSumTransactionFactory_t>();
                for (auto& refiner : popDensityBorderSumRefiners_)
                    refiner.registerLevel(hierarchy, level, border_sum_factory);
            }

            {
                auto border_max_factory = std::make_shared<VecFieldBorderMaxTransactionFactory_t>();
                for (auto& refiner : ionFluxBorderMaxRefiners_)
                    refiner.registerLevel(hierarchy, level, border_max_factory);
            }
            {
                auto border_max_factory = std::make_shared<FieldBorderMaxTransactionFactory_t>();
                for (auto& refiner : ionDensityBorderMaxRefiners_)
                    refiner.registerLevel(hierarchy, level, border_max_factory);
            }
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
                ionChargeFluxSynchronizers_.registerLevel(hierarchy, level);
            }
        }



        /**
         * @brief regrid performs the regriding communications for Hybrid to Hybrid messengers
         , all quantities that are in initialization refiners need to be regridded
         */
        void regrid(hierarchy_ptr const& hierarchy, int const levelNumber,
                    level_ptr const& oldLevel, IPhysicalModel& model,
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

            setNaNsOnVecfieldGhosts(B, level);
            magGhostsRefiners_.fill(B, level.getLevelNumber(), fillTime);
            syncInnerGhostsOnVecfield(B, level);
        }

        void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillElectricGhosts");

            setNaNsOnVecfieldGhosts(E, level);
            elecGhostsRefiners_.fill(E, level.getLevelNumber(), fillTime);
            syncInnerGhostsOnVecfield(E, level);
        }




        /**
         * @brief fillIonGhostParticles will fill the interior ghost particle array from
         * neighbor patches of the same level. Before doing that, it empties the array for
         * all populations
         */
        void fillIonGhostParticles(IonsT& ions, level_t& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(2, "HybridHybridMessengerStrategy::fillIonGhostParticles");

            domainGhostPartRefiners_.fill(level.getLevelNumber(), fillTime);

            patch_exec_(
                level,
                [](auto& ions_p) {
                    for (auto& pop : ions_p)
                        pop.patchGhostParticles().clear();
                },
                ions);
        }


        // bodies in hybrid_hybrid_messenger_strategy_borders.hpp
        void fillFluxBorders(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                             double const fillTime) override
        {
            fluxBorders_(ions, level, fillTime);
        }

        void fillDensityBorders(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                double const fillTime) override
        {
            densityBorders_(ions, level, fillTime);
        }

        void fillIonBorders(IonsT& ions, level_t& level, double const fillTime) override
        {
            ionBorders_(ions, level, fillTime);
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
            if (level.getLevelNumber() == 0) // no levelGhost on root level
                return;

            auto accessor = make_level_accessor<GridLayoutT>(level, *resourcesManager_, ions);
            level_exec(
                accessor,
                [](auto& view) -> auto& { return std::get<0>(view.args)[0].particleDensity(); },
                [&](auto& view, auto const&... tile_idx) {
                    auto& [ions_p] = view.args;
                    decltype(interpolate_) interpolate; // one per task, tasks run in parallel

                    for (auto& pop : ions_p)
                    {
                        auto& particleDensity = pop.particleDensity();
                        auto& flux            = pop.flux();

                        auto const deposit = [&](auto const& particles, double const coef) {
                            using enum core::ParticleType;
                            if constexpr (sizeof...(tile_idx) > 0)
                                interpolate.template particleToMesh_tile<LevelGhost>(
                                    tile_idx..., particles, particleDensity, flux, coef);
                            else
                                interpolate.template particleToMesh<LevelGhost>(
                                    particles, view.layout, particleDensity, flux, coef);
                        };

                        // grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        deposit(pop.levelGhostParticlesOld(), 1. - alpha);
                        deposit(pop.levelGhostParticlesNew(), alpha);
                    }
                });
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
        void firstStep(IPhysicalModel& /*model*/, level_t& level, double const currentTime,
                       double const prevCoarserTime, double const newCoarserTime) override
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::firstStep");

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

            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::lastStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            patch_exec_(
                level,
                [](auto& ions) {
                    for (auto& pop : ions)
                    {
                        auto& levelGhostParticlesOld = pop.levelGhostParticlesOld();
                        auto& levelGhostParticlesNew = pop.levelGhostParticlesNew();
                        auto& levelGhostParticles    = pop.levelGhostParticles();

                        std::swap(levelGhostParticlesNew, levelGhostParticlesOld);
                        levelGhostParticlesNew.clear();
                        levelGhostParticles = levelGhostParticlesOld;
                    }
                },
                hybridModel.state.ions);
        }




        /**
         * @brief prepareStep is the concrete implementation of the
         * HybridMessengerStrategy::prepareStep method For hybrid-Hybrid communications.
         * This method copies the electromagnetic field defined at t=n. Since prepareStep()
         * is called just before advancing the level, this operation actually saves the t=n
         * versions of B and E into the messenger. When the time comes that the next finer
         * level needs to time interpolate these quantities at its ghost nodes, this level
         * will be able to interpolate at required time because the t=n fields of previous
         * next coarser step will be in the messenger.
         */
        void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) override
        {
            PHARE_LOG_SCOPE(1, "HybridHybridMessengerStrategy::prepareStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level) // SAMRAI patch data times, serial
            {
                resourcesManager_->setTime(Bold_, *patch, currentTime);
                resourcesManager_->setTime(Eold_, *patch, currentTime);
            }

            level_exec_(
                level, [](auto&, auto& Bold, auto&...) -> auto& { return Bold; },
                [](auto& electromag, auto& Bold, auto& Eold, auto const&... tile_idx) {
                    core::copy_data(Bold, electromag.B, tile_idx...);
                    core::copy_data(Eold, electromag.E, tile_idx...);
                },
                hybridModel.state.electromag, Bold_, Eold_);
        }




        void fillRootGhosts(IPhysicalModel& model, level_t& level,
                            double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();
            assert(levelNumber == 0);

            auto& hybridModel = static_cast<HybridModel&>(model);

            // TODO : still useful?
            elecGhostsRefiners_.fill(hybridModel.state.electromag.E, levelNumber, initDataTime);

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
            ionChargeFluxSynchronizers_.sync(levelNumber);
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
            elecGhostsRefiners_.fill(hybridModel.state.electromag.E, levelNumber, time);

            chargeDensityPatchGhostsRefiners_.fill(levelNumber, time);
        }

    private:
        template<typename>
        friend struct HybridMessengerStrategyBorders;

        void fluxBorders_(IonsT& ions, SAMRAI::hier::PatchLevel& level, double const fillTime);
        void densityBorders_(IonsT& ions, SAMRAI::hier::PatchLevel& level, double const fillTime);
        void ionBorders_(IonsT& ions, level_t& level, double const fillTime);

        // Select the field-refinement operators once at construction. The composite runtime
        // kernels are built from the configured order; B uses the stage-1 magnetic kernel (fills
        // every fine face; the magnetic patch strategy runs the stage-2 divB touch-up afterward).
        // Particle refine operators (interior / level-ghost) are NOT touched.
        void makeRefineOperators_(RefinementConfig const& config)
        {
            auto fieldKernel = [&] {
                return std::make_shared<KernelFieldRefineOperator<GridLayoutT, GridT>>(
                    makeRefineKernel<GridLayoutT, refine_kernel_field_t<GridT>>(config.order));
            };
            auto vecKernel = [&] {
                return std::make_shared<KernelVecFieldRefineOperator<VectorFieldDataT>>(
                    makeRefineKernel<GridLayoutT, refine_kernel_field_t<GridT>>(config.order));
            };
            auto magKernel = [&] {
                return std::make_shared<KernelVecFieldRefineOperator<VectorFieldDataT>>(
                    makeMagneticRefineKernel<GridLayoutT, refine_kernel_field_t<GridT>>(
                        config.order));
            };

            fieldRefineOp_  = fieldKernel();
            BInitRefineOp_  = magKernel();
            BRefineOp_      = magKernel();
            EfieldRefineOp_ = vecKernel();
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
            elecGhostsRefiners_.addTimeRefiners(info->ghostElectric, info->modelElectric,
                                                Eold_.name(), EfieldRefineOp_, vecFieldTimeOp_,
                                                nonOverwriteInteriorTFfillPattern);


            // we need a separate patch strategy for each refiner so that each one can register
            // their required ids
            magneticPatchStratPerGhostRefiner_ = [&]() {
                std::vector<std::shared_ptr<
                    MagneticRefinePatchStrategy<VectorFieldDataT, BoundaryManagerT>>>
                    result;

                result.reserve(info->ghostMagnetic.size());

                for (auto const& key : info->ghostMagnetic)
                {
                    auto&& [id] = resourcesManager_->getIDsList(key);

                    auto patch_strat = std::make_shared<
                        MagneticRefinePatchStrategy<VectorFieldDataT, BoundaryManagerT>>(
                        *boundaryManager_);

                    patch_strat->registerIDs(id);
                    patch_strat->registerScratchID(scratchID_(scratch_vecfield()));

                    result.push_back(patch_strat);
                }
                return result;
            }();

            // B (and Bpred) C-F ghosts are time interpolated between Bold_ (t=n) and model
            // B (t=n+1). Per-refiner magnetic patch strategy is carried through so the
            // divB touch-up still runs after spatial refine.
            for (std::size_t i = 0; i < info->ghostMagnetic.size(); ++i)
            {
                magGhostsRefiners_.addTimeRefiner(
                    info->ghostMagnetic[i], info->modelMagnetic, Bold_.name(), BRefineOp_,
                    vecFieldTimeOp_, info->ghostMagnetic[i], nonOverwriteInteriorTFfillPattern,
                    magneticPatchStratPerGhostRefiner_[i]);
            }


            // J has no refiner: it is never communicated, each level computes it
            // locally via Ampere on physical+1 from the time-refined B ghosts.

            // same level schedule: patch ghost nodes are copied from neighbor patches at the
            // same time, no time interpolation
            chargeDensityPatchGhostsRefiners_.addStaticRefiner(
                info->modelIonDensity, fieldRefineOp_, info->modelIonDensity,
                nonOverwriteInteriorFieldFillPattern);
        }




        void registerInitComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
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


            domainGhostPartRefiners_.addStaticRefiners(
                info->patchGhostParticles, nullptr, info->patchGhostParticles,
                std::make_shared<ParticleDomainFromGhostFillPattern<GridLayoutT>>());


            // the model's temporaries, see requestTemporaries
            auto const tmpVecFieldName = model_.getTmp(core::HybridQuantity::Vector::V).name();
            auto const tmpFieldName    = model_.getTmp(core::QtyCentering::primal).name();

            for (auto const& vecfield : info->ghostFlux)
                popFluxBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        tmpVecFieldName, vecfield, nullptr, tmpVecFieldName,
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());

            for (auto const& field : info->sumBorderFields)
                popDensityBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        tmpFieldName, field, nullptr, tmpFieldName,
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());


            assert(info->maxBorderFields.size() == 1); // charge density
            for (auto const& field : info->maxBorderFields)
                ionDensityBorderMaxRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        tmpFieldName, field, nullptr, tmpFieldName,
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());

            assert(info->maxBorderVecFields.size() == 1); // charge flux
            for (auto const& vecfield : info->maxBorderVecFields)
                ionFluxBorderMaxRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        tmpVecFieldName, vecfield, nullptr, tmpVecFieldName,
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
        }



        // tiled field refines borrow the untiled scratch fields of the destination patch
        int scratchID_(auto const& scratch) const
        {
            if (auto const id = resourcesManager_->getID(scratch.name()))
                return *id;
            throw std::runtime_error("HybridHybridMessengerStrategy: unregistered scratch field "
                                     + scratch.name());
        }


        void registerSyncComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            electroSynchronizers_.add(info->modelElectric, electricFieldCoarseningOp_,
                                      info->modelElectric);

            // no patch ghost refill after sync for the charge flux: the electron velocity
            // is only computed and used on the physical box
            ionChargeFluxSynchronizers_.add(info->modelIonChargeFlux, vecFieldMomentsCoarseningOp_,
                                            info->modelIonChargeFlux);

            chargeDensitySynchronizers_.add(info->modelIonDensity, fieldMomentsCoarseningOp_,
                                            info->modelIonDensity);
        }




        void copyLevelGhostOldToPushable_(level_t& level, IPhysicalModel& model)
        {
            auto& hybridModel = static_cast<HybridModel&>(model);
            patch_exec_(
                level,
                [](auto& ions) {
                    for (auto& pop : ions)
                        pop.levelGhostParticles() = pop.levelGhostParticlesOld();
                },
                hybridModel.state.ions);
        }




        double timeInterpCoef_(double const afterPushTime, std::size_t levelNumber)
        {
            return (afterPushTime - beforePushCoarseTime_[levelNumber])
                   / (afterPushCoarseTime_[levelNumber] - beforePushCoarseTime_[levelNumber]);
        }




        void magneticRegriding_(hierarchy_ptr const& hierarchy, level_ptr const& level,
                                level_ptr const& oldLevel, HybridModel& hybridModel,
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
         * The ghost refiners (magGhostsRefiners_, elecGhostsRefiners_) run one
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
        void setNaNsOnFieldGhosts(FieldT& field, patch_t const& patch)
        {
            core::fill_ghost(field, layoutFromPatch<GridLayoutT>(patch),
                             std::numeric_limits<typename FieldT::type>::quiet_NaN());
        }

        void setNaNsOnFieldGhosts(FieldT& field, level_t const& level)
        {
            for (auto& patch : resourcesManager_->enumerate(level, field))
                setNaNsOnFieldGhosts(field, *patch);
        }

        void setNaNsOnVecfieldGhosts(VecFieldT& vf, level_t const& level)
        {
            for (auto& patch : resourcesManager_->enumerate(level, vf))
                for (auto& field : vf)
                    setNaNsOnFieldGhosts(field, *patch);
        }

        void syncInnerGhostsOnVecfield(VecFieldT& vf, level_t const& level)
        {
            using field_type = typename VecFieldT::field_type;
            if constexpr (core::is_field_tile_set_v<field_type>)
                level_exec_(
                    level, [](auto& vf_p) -> auto& { return vf_p; },
                    [](auto& vf_p, auto const&... tile_idx) {
                        core::sync_inner_ghosts(vf_p, tile_idx...);
                    },
                    vf);
        }

        // fn(resources...) per patch of the level, a patch per idle pool
        void patch_exec_(auto& level, auto const& fn, auto&... resources)
        {
            auto accessor
                = make_level_accessor<GridLayoutT>(level, *resourcesManager_, resources...);
            patch_level_exec(accessor, [&](auto& view) { std::apply(fn, view.args); });
        }

        // fn(resources..., tile_idx) per tile when tiled(resources...) is tiled (a patch per idle
        // pool, its tiles on that pool's threads), else fn(resources...) per patch
        void level_exec_(auto& level, auto const& tiled, auto const& fn, auto&... resources)
        {
            auto accessor
                = make_level_accessor<GridLayoutT>(level, *resourcesManager_, resources...);
            level_exec(
                accessor, [&](auto& view) -> auto& { return std::apply(tiled, view.args); },
                [&](auto& view, auto const&... tile_idx) {
                    std::apply([&](auto&... r) { fn(r..., tile_idx...); }, view.args);
                });
        }


        VecFieldT Bold_{stratName + "_Bold", core::HybridQuantity::Vector::B};
        VecFieldT Eold_{stratName + "_Eold", core::HybridQuantity::Vector::E};

        HybridModel& model_; // owns the temporaries used here

        // the model's temporaries, never tiled, see requestTemporaries in the constructor
        //  looked up on use, the model can reallocate them on later requests
        Field_t& scratch_field() const { return model_.getTmp(core::QtyCentering::primal); }
        VecField_t& scratch_vecfield() const
        {
            return model_.getTmp(core::HybridQuantity::Vector::V);
        }


        //! ResourceManager shared with other objects (like the HybridModel)
        std::shared_ptr<ResourcesManagerT> resourcesManager_;
        std::shared_ptr<BoundaryManagerT> boundaryManager_;


        int const firstLevel_;
        std::unordered_map<std::size_t, double> beforePushCoarseTime_;
        std::unordered_map<std::size_t, double> afterPushCoarseTime_;

        core::Interpolating<dimension, interpOrder> interpolate_;

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
        using FieldGhostMaxRefinerPool    = RefinerPool<rm_t, RefinerType::PatchFieldBorderMax>;
        using VecFieldGhostMaxRefinerPool = RefinerPool<rm_t, RefinerType::PatchVecFieldBorderMax>;
        using FieldFillPattern_t          = FieldFillPattern<dimension>;
        using TensorFieldFillPattern_t    = TensorFieldFillPattern<dimension /*, rank=1*/>;

        //! += flux on ghost box overlap incomplete population moment nodes
        std::vector<VecFieldGhostSumRefinerPool> popFluxBorderSumRefiners_;
        //! += density on ghost box overlap incomplete population moment nodes
        std::vector<FieldGhostSumRefinerPool> popDensityBorderSumRefiners_;

        std::vector<FieldGhostMaxRefinerPool> ionDensityBorderMaxRefiners_;
        std::vector<VecFieldGhostMaxRefinerPool> ionFluxBorderMaxRefiners_;

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

        //! store refiners for electric fields that need ghosts to be filled
        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};

        GhostRefinerPool magGhostsRefiners_{resourcesManager_};

        // moment ghosts
        // There is no level ghost refiner for the moments: on level ghost nodes they come
        // only from the deposit of level ghost particles, which leaves the outermost level
        // ghost nodes incomplete (missing the outside contribution). A level ghost schedule
        // cannot replace that deposit: at interp order >= 2, level ghost particles also
        // contribute to domain nodes. Such a schedule would be needed to fill the pure level
        // ghost nodes of the charge density and charge flux, or to coarsen moments with a
        // multiple point operator.

        PatchGhostRefinerPool chargeDensityPatchGhostsRefiners_{resourcesManager_};

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


        DomainGhostPartRefinerPool domainGhostPartRefiners_{resourcesManager_};

        SynchronizerPool<rm_t> chargeDensitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> ionChargeFluxSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> electroSynchronizers_{resourcesManager_};


        // built in the ctor body (makeRefineOperators_) from the composite runtime kernels.
        RefOp_ptr fieldRefineOp_;


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

        std::shared_ptr<TimeInterpolateOperator> vecFieldTimeOp_{
            std::make_shared<VecFieldTimeInterp>()};

        using CoarsenOperator_ptr = std::shared_ptr<SAMRAI::hier::CoarsenOperator>;

        CoarsenOperator_ptr fieldMomentsCoarseningOp_{std::make_shared<MomentsFieldCoarsenOp>()};
        CoarsenOperator_ptr vecFieldMomentsCoarseningOp_{
            std::make_shared<MomentsVecFieldCoarsenOp>()};
        CoarsenOperator_ptr electricFieldCoarseningOp_{std::make_shared<ElectricFieldCoarsenOp>()};

        MagneticRefinePatchStrategy<VectorFieldDataT, BoundaryManagerT>
            magneticRefinePatchStrategy_{*boundaryManager_};

        std::vector<
            std::shared_ptr<MagneticRefinePatchStrategy<VectorFieldDataT, BoundaryManagerT>>>
            magneticPatchStratPerGhostRefiner_;
    };


} // namespace amr

} // namespace PHARE

#endif
