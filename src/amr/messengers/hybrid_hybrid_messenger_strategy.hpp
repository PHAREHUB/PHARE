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
#include "amr/physical_models/fully_kinetic_model.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/particles/particles_variable_fill_pattern.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "mpi/mpi_utils.hpp"

#include <SAMRAI/hier/IntVector.h>
#include <SAMRAI/hier/Patch.h>
#include <SAMRAI/xfer/RefineSchedule.h>
#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/hier/CoarseFineBoundary.h>
#include <SAMRAI/xfer/BoxGeometryVariableFillPattern.h>

#include <memory>
#include <algorithm>
#include <iterator>
#include <limits>
#include <string>
#include <vector>
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
        using amr_types     = HybridModel::amr_types;
        using level_t       = amr_types::level_t;
        using patch_t       = amr_types::patch_t;
        using hierarchy_t   = amr_types::hierarchy_t;
        using level_ptr     = amr_types::level_ptr;
        using hierarchy_ptr = amr_types::hierarchy_ptr;

        using GridT             = HybridModel::grid_type;
        using IonsT             = HybridModel::ions_type;
        using ElectromagT       = HybridModel::electromag_type;
        using VecFieldT         = HybridModel::vecfield_type;
        using TensorFieldT      = IonsT::tensorfield_type;
        using GridLayoutT       = HybridModel::gridlayout_type;
        using FieldT            = VecFieldT::field_type;
        using VectorFieldDataT  = TensorFieldData<1, GridLayoutT, GridT, core::HybridQuantity>;
        using ResourcesManagerT = HybridModel::resources_manager_type;
        using BoundaryManagerT  = HybridModel::boundary_manager_type;
        using IPhysicalModel    = HybridModel::Interface;

        static constexpr std::size_t dimension   = GridLayoutT::dimension;
        static constexpr std::size_t interpOrder = GridLayoutT::options.interp_order;

        // the fully kinetic model deposits its current
        static constexpr bool sumsParticleCurrent = solver::is_fully_kinetic_model_v<HybridModel>;

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
        static inline std::string const stratName = std::string{HybridModel::model_type_name} + "-"
                                                    + std::string{HybridModel::model_type_name};
        static constexpr std::size_t rootLevelNumber = 0;


        HybridHybridMessengerStrategy(std::shared_ptr<ResourcesManagerT> const& manager,
                                      std::shared_ptr<BoundaryManagerT> const& boundaryManager,
                                      int const firstLevel,
                                      RefinementConfig const& refinementConfig = {})

            : HybridMessengerStrategy<HybridModel>{stratName}
            , resourcesManager_{manager}
            , boundaryManager_{boundaryManager}
            , firstLevel_{firstLevel}
        {
            makeRefineOperators_(refinementConfig);

            resourcesManager_->registerResources(Bold_);
            resourcesManager_->registerResources(Eold_);
            resourcesManager_->registerResources(sumVec_);
            resourcesManager_->registerResources(sumField_);
            resourcesManager_->registerResources(sumTensor_);
            if constexpr (sumsParticleCurrent)
                resourcesManager_->registerResources(sumJ_);
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
            resourcesManager_->allocate(sumVec_, patch, allocateTime);
            resourcesManager_->allocate(sumField_, patch, allocateTime);
            resourcesManager_->allocate(sumTensor_, patch, allocateTime);
            if constexpr (sumsParticleCurrent)
                resourcesManager_->allocate(sumJ_, patch, allocateTime);
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

            // fully kinetic partial reflux (B-x): the same coarsening and patch ghost fill as the
            // reflux, from the finer level's extrapolated covered E onto a copy of the coarse E
            if constexpr (sumsParticleCurrent)
                if (!hybridInfo->partialElectric.empty())
                {
                    auto&& [e_partial_id, e_partial_flux_id] = resourcesManager_->getIDsList(
                        hybridInfo->partialElectric, hybridInfo->partialFluxElectric);
                    PartialRefluxAlgo.registerCoarsen(e_partial_id, e_partial_flux_id,
                                                      electricFieldCoarseningOp_);
                    PatchGhostPartialAlgo.registerRefine(e_partial_id, e_partial_id, e_partial_id,
                                                         EfieldRefineOp_,
                                                         nonOverwriteInteriorTFfillPattern);
                    partialReflux_ = true;
                }

            registerGhostComms_(hybridInfo);
            registerInitComms(hybridInfo);
            registerSyncComms(hybridInfo);
        }



        void registerLevel(hierarchy_ptr const& hierarchy, int const levelNumber) override
        {
            // the level is new: drop the schedules of the level it replaces, rebuild same-level
            sameLevelBuiltOn_.erase(levelNumber);
            interLevelBuiltOn_.erase(levelNumber);
            removeInterLevelSchedules_(levelNumber);
            buildSameLevelSchedules_(hierarchy, levelNumber);
            sameLevelBuiltOn_[levelNumber] = levelsOf_(hierarchy, levelNumber, levelNumber);

            // inter-level schedules read the coarser levels, which may still be replaced before
            // the hierarchy settles (e.g. a level-0 rebalance), so onHierarchyChanged builds them
        }



        // Rebuilds the schedules of levelNumber that were built from a level since replaced.
        // Building a schedule is MPI-collective: all ranks take the same branches, because
        // SAMRAI replaces a level on all ranks at once.
        void onHierarchyChanged(hierarchy_ptr const& hierarchy, int const levelNumber) override
        {
            if (!builtOnCurrentLevels_(sameLevelBuiltOn_, hierarchy, levelNumber, levelNumber))
            {
                buildSameLevelSchedules_(hierarchy, levelNumber);
                sameLevelBuiltOn_[levelNumber] = levelsOf_(hierarchy, levelNumber, levelNumber);
            }

            if (levelNumber != rootLevelNumber
                and !builtOnCurrentLevels_(interLevelBuiltOn_, hierarchy, rootLevelNumber,
                                           levelNumber))
            {
                interLevelBuiltOn_.erase(levelNumber);
                removeInterLevelSchedules_(levelNumber);
                buildInterLevelSchedules_(hierarchy, levelNumber);
                interLevelBuiltOn_[levelNumber]
                    = levelsOf_(hierarchy, rootLevelNumber, levelNumber);
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
            // the fully kinetic E is at its own half step, not at initDataTime: regridElectric
            if constexpr (!sumsParticleCurrent)
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
                lvlGhostPartOldRefiners_.fillTransient(hierarchy, level, initDataTime);
                copyLevelGhostOldToPushable_(*level, model);
            }

            // levelGhostNew will be refined in next firstStep

            // after filling the new level with the regrid schedule, some
            // nodes may not have been copied correctly, due to a bug in SAMRAI
            // it seems these nodes are only on ghost box border if that border
            // overlaps an old level patch border. See https://github.com/LLNL/SAMRAI/pull/293
        }

        // fully kinetic (P5a): the E of a level (re)created during the run, at fillTime, its own
        // half step. On a root rebalance (no coarser level) everything is copied. Above, the level
        // ghosts are filled again after the regrid: a regrid schedule can leave nodes of the ghost
        // box border uncopied where it meets an old patch border (SAMRAI PR 293), and a level
        // ghost takes the coarser E, as every later fill does
        void regridElectric(VecFieldT& E, hierarchy_ptr const& hierarchy, int const levelNumber,
                            level_ptr const& oldLevel, double const fillTime) override
        {
            if constexpr (sumsParticleCurrent)
            {
                if (elecRegridRefiners_.empty())
                    throw std::runtime_error("regridElectric: no E one step ahead registered");
                elecRegridRefiners_.regrid(hierarchy, levelNumber, oldLevel, fillTime);
                if (levelNumber == rootLevelNumber)
                    return;
                auto const level = hierarchy->getPatchLevel(levelNumber);
                setNaNsOnVecfieldGhosts(E, *level);
                elecRegridGhostsRefiners_.fillTransient(hierarchy, level, fillTime);
            }
            else
                throw std::runtime_error("regridElectric: fully kinetic only");
        }

        // fully kinetic: the model E ghosts other patches of the level own (the same-level
        // schedule of EalgoPatchGhost), the level ghosts kept
        void fillElectricPatchGhosts(VecFieldT& /*E*/, level_t const& level,
                                     double const fillTime) override
        {
            elecPatchGhostsRefineSchedules.at(level.getLevelNumber())->fillData(fillTime);
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
        void initLevel(hierarchy_ptr const& hierarchy, IPhysicalModel& model, level_t& level,
                       double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();

            auto& hybridModel = static_cast<HybridModel&>(model);

            auto const levelPtr = hierarchy->getPatchLevel(levelNumber);

            // init-only schedules are used once, so they are built at fill time, not stored
            BalgoInit
                .createSchedule(levelPtr, nullptr, levelNumber - 1, hierarchy,
                                &magneticRefinePatchStrategy_)
                ->fillData(initDataTime);
            electricInitRefiners_.fillTransient(hierarchy, levelPtr, initDataTime);

            // no need to call these :
            // magGhostsRefiners_.fill(levelNumber, initDataTime);
            // elecGhostsRefiners_.fill(levelNumber, initDataTime);
            // because the SAMRAI schedules in the 'init' communicators
            // already fill the patch ghost box from the neighbor interior box.
            // so ghost nodes are already filled .

            PHARE_LOG_START(3, "hybhybmessengerStrat::initLevel : interior part fill schedule");
            domainParticlesRefiners_.fillTransient(hierarchy, levelPtr, initDataTime);
            PHARE_LOG_STOP(3, "hybhybmessengerStrat::initLevel : interior part fill schedule");

            lvlGhostPartOldRefiners_.fillTransient(hierarchy, levelPtr, initDataTime);


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
            auto const levelNumber = level.getLevelNumber();
            if constexpr (sumsParticleCurrent)
                if (fillsMid_(B, levelNumber, fillTime))
                    return magGhostsMidRefiners_.fill(B, levelNumber, fillTime);
            magGhostsRefiners_.fill(B, levelNumber, fillTime);
        }

        void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillElectricGhosts");

            setNaNsOnVecfieldGhosts(E, level);
            auto const levelNumber = level.getLevelNumber();
            if constexpr (sumsParticleCurrent)
                if (fillsAhead_(levelNumber, fillTime))
                    return elecGhostsAheadRefiners_.fill(E, levelNumber, fillTime);
            elecGhostsRefiners_.fill(E, levelNumber, fillTime);
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

            // one particle density per population, the charge density is computed from it
            assert(popDensityBorderSumRefiners_.size() == ions.size());

            for (std::size_t i = 0; i < ions.size(); ++i)
            {
                for (auto patch : resourcesManager_->enumerate(level, ions, sumField_))
                    std::memcpy(sumField_.data(), ions[i].particleDensity().data(),
                                ions[i].particleDensity().size() * sizeof(value_type));

                popDensityBorderSumRefiners_[i].fill(level.getLevelNumber(), fillTime);

                for (auto patch : resourcesManager_->enumerate(level, ions, sumField_))
                    std::memcpy(ions[i].particleDensity().data(), sumField_.data(),
                                ions[i].particleDensity().size() * sizeof(value_type));
            }
        }


        // same copy, sum, copy back as fillFluxBorders, on the Yee centered current
        void fillCurrentBorders(VecFieldT& J, level_t& level, double const fillTime) override
        {
            if constexpr (sumsParticleCurrent)
            {
                auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();
                using value_type = FieldT::value_type;
                assert(currentBorderSumRefiners_.size() == 1);

                for (auto patch : resourcesManager_->enumerate(level, J, sumJ_))
                    for (std::uint8_t c = 0; c < N; ++c)
                        std::memcpy(sumJ_[c].data(), J[c].data(), J[c].size() * sizeof(value_type));

                currentBorderSumRefiners_[0].fill(level.getLevelNumber(), fillTime);

                for (auto patch : resourcesManager_->enumerate(level, J, sumJ_))
                    for (std::uint8_t c = 0; c < N; ++c)
                        std::memcpy(J[c].data(), sumJ_[c].data(), J[c].size() * sizeof(value_type));
            }
            else
                HybridMessengerStrategy<HybridModel>::fillCurrentBorders(J, level, fillTime);
        }

        // copies the ghost nodes of J from the patches that own them, after the border sum
        void fillCurrentGhosts(VecFieldT& J, level_t& level, double const fillTime) override
        {
            if constexpr (sumsParticleCurrent)
                currentGhostsRefiners_.fill(J, level.getLevelNumber(), fillTime);
            else
                HybridMessengerStrategy<HybridModel>::fillCurrentGhosts(J, level, fillTime);
        }

        // same as the current: copy, sum, copy back, through sumField_ (primal like rho, a
        // scratch the population densities only use inside fillDensityBorders)
        void fillChargeDensityBorders(FieldT& rho, level_t& level, double const fillTime) override
        {
            if constexpr (sumsParticleCurrent)
            {
                using value_type = FieldT::value_type;
                assert(chargeDensityBorderSumRefiners_.size() == 1);

                for (auto patch : resourcesManager_->enumerate(level, rho, sumField_))
                    std::memcpy(sumField_.data(), rho.data(), rho.size() * sizeof(value_type));

                chargeDensityBorderSumRefiners_[0].fill(level.getLevelNumber(), fillTime);

                for (auto patch : resourcesManager_->enumerate(level, rho, sumField_))
                    std::memcpy(rho.data(), sumField_.data(), rho.size() * sizeof(value_type));
            }
            else
                HybridMessengerStrategy<HybridModel>::fillChargeDensityBorders(rho, level,
                                                                               fillTime);
        }

        // copies the ghost nodes of rho from the patches that own them, after the border sum
        void fillChargeDensityGhosts(FieldT& rho, level_t& level, double const fillTime) override
        {
            if constexpr (sumsParticleCurrent)
                chargeDensityGhostsRefiners_.fill(rho, level.getLevelNumber(), fillTime);
            else
                HybridMessengerStrategy<HybridModel>::fillChargeDensityGhosts(rho, level, fillTime);
        }

        void fillIonBorders(IonsT& /*ions*/, level_t& level, double const fillTime) override
        {
            assert(ionFluxBorderMaxRefiners_.size() == 1);    // charge flux
            assert(ionDensityBorderMaxRefiners_.size() == 1); // charge density

            for (auto& refiner : ionFluxBorderMaxRefiners_)
                refiner.fill(level.getLevelNumber(), fillTime);
            for (auto& refiner : ionDensityBorderMaxRefiners_)
                refiner.fill(level.getLevelNumber(), fillTime);
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
                    auto& flux            = pop.flux();

                    if (level.getLevelNumber() > 0) // no levelGhost on root level
                    {
                        // grab levelGhostParticlesOld and levelGhostParticlesNew
                        // and project them with alpha and (1-alpha) coefs, respectively
                        auto& levelGhostOld = pop.levelGhostParticlesOld();
                        interpolate_(makeRange(levelGhostOld), particleDensity, flux, layout,
                                     1. - alpha);

                        auto& levelGhostNew = pop.levelGhostParticlesNew();
                        interpolate_(makeRange(levelGhostNew), particleDensity, flux, layout,
                                     alpha);
                    }
                }
            }
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
         * This method copies the electromagnetic field defined at t=n. Since prepareStep()
         * is called just before advancing the level, this operation actually saves the t=n
         * versions of B and E into the messenger. When the time comes that the next finer
         * level needs to time interpolate these quantities at its ghost nodes, this level
         * will be able to interpolate at required time because the t=n fields of previous
         * next coarser step will be in the messenger.
         */
        void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) override
        {
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::prepareStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, hybridModel.state.electromag, hybridModel.state.J,
                    hybridModel.state.ions, Bold_, Eold_);

                resourcesManager_->setTime(Bold_, *patch, currentTime);
                resourcesManager_->setTime(
                    Eold_, *patch,
                    electricTime_(hybridModel.state.electromag.E, *patch, currentTime));

                Bold_.copyData(hybridModel.state.electromag.B);
                Eold_.copyData(hybridModel.state.electromag.E);
            }
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

            // call coarsning schedules... The fully kinetic E is synchronized by reflux instead:
            // the finer model E is at another time (t - dt_f/2), its time average at the coarse E's
            // (P4 plan D4)
            if constexpr (!sumsParticleCurrent)
                electroSynchronizers_.sync(levelNumber);
            chargeDensitySynchronizers_.sync(levelNumber);
            ionChargeFluxSynchronizers_.sync(levelNumber);
        }

        // this function coarsens the fluxSum onto the corresponding coarser fluxes (E in hybrid),
        // and fills the patch ghosts, making it ready for the faraday in the solver.reflux()
        void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                    double const syncTime) override
        {
            refluxSchedules.at(fineLevelNumber)->coarsenData();
            patchGhostRefluxedSchedules.at(coarserLevelNumber)->fillData(syncTime);
        }

        // fully kinetic (B-x): coarsens the finer level's covered E, extrapolated to the coarse E
        // time by the solver, onto the coarser level's copy of its E, and refills that copy's
        // patch ghosts, as reflux does for the model E. The solver then builds the coarse B at
        // midTime, from which the finer level's B ghosts are filled at that time
        void partialReflux(int const coarserLevelNumber, int const fineLevelNumber,
                           double const midTime) override
        {
            if constexpr (sumsParticleCurrent)
            {
                if (!partialReflux_)
                    throw std::runtime_error("partialReflux: no partial reflux registered");
                partialRefluxSchedules.at(fineLevelNumber)->coarsenData();
                patchGhostPartialSchedules.at(coarserLevelNumber)->fillData(midTime);
                midMagneticTime_[fineLevelNumber] = midTime;
            }
            else
                throw std::runtime_error("partialReflux: fully kinetic only");
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

            auto& E = hybridModel.state.electromag.E;

            // Fully kinetic, above the root (three levels or more: the level synchronized here
            // has a coarser one): its E sits at its own half step, time - dt_L/2, not at the
            // synchronization time. Filled at the synchronization time, the level ghosts would be
            // interpolated from the coarser level at the wrong time, and from the wrong pair past
            // the coarser E (the ahead E is only picked by fillElectricGhosts): at the end of the
            // coarser step, time > its E^{n+1/2}, SAMRAI's time bracket fails. So they are filled
            // like every other E fill, at E's own label. On the root (no coarser level) the fill
            // is same-level only and its time is not used, as for the hybrid (E at time there).
            if constexpr (sumsParticleCurrent)
                if (levelNumber != rootLevelNumber)
                {
                    fillElectricGhosts(E, level, electricTimeOn_(E, level));
                    chargeDensityPatchGhostsRefiners_.fill(levelNumber, time);
                    return;
                }

            // this electric schedule should probably only be a patch ghost one
            // since levelghost nodes are not affected by the coarsening
            elecGhostsRefiners_.fill(E, levelNumber, time);

            chargeDensityPatchGhostsRefiners_.fill(levelNumber, time);
        }

    private:
        void buildSameLevelSchedules_(hierarchy_ptr const& hierarchy, int const levelNumber)
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magPatchGhostsRefineSchedules[levelNumber]
                = BalgoPatchGhost.createSchedule(level, &magneticRefinePatchStrategy_);

            elecPatchGhostsRefineSchedules[levelNumber] = EalgoPatchGhost.createSchedule(level);

            // technically not needed for finest as refluxing is not done onto it.
            patchGhostRefluxedSchedules[levelNumber] = PatchGhostRefluxedAlgo.createSchedule(level);
            if (partialReflux_)
                patchGhostPartialSchedules[levelNumber]
                    = PatchGhostPartialAlgo.createSchedule(level);

            // on the root level there is no coarser source, so these are same-level schedules
            if (levelNumber == rootLevelNumber)
            {
                elecGhostsRefiners_.registerLevel(hierarchy, level);
                magGhostsRefiners_.registerLevel(hierarchy, level);
            }
            if constexpr (sumsParticleCurrent)
            {
                currentGhostsRefiners_.registerLevel(hierarchy, level);
                chargeDensityGhostsRefiners_.registerLevel(hierarchy, level);
            }
            domainGhostPartRefiners_.registerLevel(hierarchy, level);

            chargeDensityPatchGhostsRefiners_.registerLevel(hierarchy, level);

            for (auto& refiner : popFluxBorderSumRefiners_)
                refiner.registerLevel(hierarchy, level);

            for (auto& refiner : popDensityBorderSumRefiners_)
                refiner.registerLevel(hierarchy, level);

            for (auto& refiner : currentBorderSumRefiners_)
                refiner.registerLevel(hierarchy, level);

            for (auto& refiner : chargeDensityBorderSumRefiners_)
                refiner.registerLevel(hierarchy, level);

            for (auto& refiner : ionFluxBorderMaxRefiners_)
                refiner.registerLevel(hierarchy, level);
            for (auto& refiner : ionDensityBorderMaxRefiners_)
                refiner.registerLevel(hierarchy, level);
        }



        void buildInterLevelSchedules_(hierarchy_ptr const& hierarchy, int const levelNumber)
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            elecGhostsRefiners_.registerLevel(hierarchy, level);
            elecGhostsAheadRefiners_.registerLevel(hierarchy, level);
            magGhostsRefiners_.registerLevel(hierarchy, level);
            magGhostsMidRefiners_.registerLevel(hierarchy, level);

            // refluxing
            auto const& coarseLevel      = hierarchy->getPatchLevel(levelNumber - 1);
            refluxSchedules[levelNumber] = RefluxAlgo.createSchedule(coarseLevel, level);
            if (partialReflux_)
                partialRefluxSchedules[levelNumber]
                    = PartialRefluxAlgo.createSchedule(coarseLevel, level);

            lvlGhostPartNewRefiners_.registerLevel(hierarchy, level);

            // and these for coarsening
            electroSynchronizers_.registerLevel(hierarchy, level);
            chargeDensitySynchronizers_.registerLevel(hierarchy, level);
            ionChargeFluxSynchronizers_.registerLevel(hierarchy, level);
        }



        // the levels lo..hi of the hierarchy, recorded against the schedules built from them
        static auto levelsOf_(hierarchy_ptr const& hierarchy, int const lo, int const hi)
        {
            std::vector<std::weak_ptr<level_t>> levels;
            for (auto ilvl = lo; ilvl <= hi; ++ilvl)
                levels.emplace_back(hierarchy->getPatchLevel(ilvl));
            return levels;
        }

        // whether the schedules of level hi were built against the current levels lo..hi
        static bool builtOnCurrentLevels_(auto const& builtOn, hierarchy_ptr const& hierarchy,
                                          int const lo, int const hi)
        {
            auto const it = builtOn.find(hi);
            if (it == builtOn.end())
                return false;
            for (auto ilvl = lo; ilvl <= hi; ++ilvl)
                if (it->second[ilvl - lo].lock() != hierarchy->getPatchLevel(ilvl))
                    return false;
            return true;
        }

        void removeInterLevelSchedules_(int const levelNumber)
        {
            elecGhostsRefiners_.removeLevel(levelNumber);
            elecGhostsAheadRefiners_.removeLevel(levelNumber);
            magGhostsRefiners_.removeLevel(levelNumber);
            magGhostsMidRefiners_.removeLevel(levelNumber);
            lvlGhostPartNewRefiners_.removeLevel(levelNumber);
            electroSynchronizers_.removeLevel(levelNumber);
            chargeDensitySynchronizers_.removeLevel(levelNumber);
            ionChargeFluxSynchronizers_.removeLevel(levelNumber);
            refluxSchedules.erase(levelNumber);
            partialRefluxSchedules.erase(levelNumber);
            midMagneticTime_.erase(levelNumber);
        }

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

            // fully kinetic, finer level E ghosts at a time past the coarse model E (E^{n+1/2}):
            // time interpolated between the model E and the E one coarse step ahead (E^{n+3/2}),
            // fillElectricGhosts picks this pool past the coarse E time. Inert until the solver
            // names the ahead E.
            if constexpr (sumsParticleCurrent)
                if (!info->aheadElectric.empty())
                    elecGhostsAheadRefiners_.addTimeRefiners(
                        info->ghostElectric, info->aheadElectric, info->modelElectric,
                        EfieldRefineOp_, vecFieldTimeOp_, nonOverwriteInteriorTFfillPattern);

            // fully kinetic, a level (re)created during the run (P5a): its E at its own half step.
            // A regrid schedule copies the source from the old level where it overlaps the new
            // one, without time interpolation (it is built without time refinement,
            // RefineSchedule.cpp:454-460), and refines from the coarser level elsewhere,
            // interpolated in time (its coarse recursion always is, :600-614) between the coarser
            // model E and the coarser E one step ahead: Eold_ is stale at regrid time. No fill
            // pattern: interior and ghosts, as electricInitRefiners_. The level ghosts are then
            // filled again from the same pair, the same-level ghosts from the new patches.
            if constexpr (sumsParticleCurrent)
                if (!info->aheadElectric.empty())
                {
                    elecRegridRefiners_.addTimeRefiner(
                        info->modelElectric, info->aheadElectric, info->modelElectric,
                        EfieldRefineOp_, vecFieldTimeOp_, info->modelElectric);
                    elecRegridGhostsRefiners_.addTimeRefiner(
                        info->modelElectric, info->aheadElectric, info->modelElectric,
                        EfieldRefineOp_, vecFieldTimeOp_, info->modelElectric,
                        nonOverwriteInteriorTFfillPattern);
                }

            // fully kinetic J (summed on patch borders, then completed by the level ghost
            // particle deposit): ghosts copied from the same-level patches owning them, never
            // refined from the coarser level, whose J is at another time and from other
            // particles; the level ghost band keeps the level ghost particle deposit (P4 plan D5).
            // the charge density of the divergence cleaning: the same, primal scalar
            if constexpr (sumsParticleCurrent)
            {
                currentGhostsRefiners_.addStaticRefiners(info->ghostCurrent, EfieldRefineOp_,
                                                         info->ghostCurrent,
                                                         nonOverwriteInteriorTFfillPattern);
                chargeDensityGhostsRefiners_.addStaticRefiners(
                    info->ghostChargeDensity, fieldRefineOp_, info->ghostChargeDensity,
                    nonOverwriteInteriorFieldFillPattern);
            }


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

            // fully kinetic (B-x): the model B ghosts at the end of a partially refluxed substep,
            // between Bold_ and the coarse mid B labelled with that time: exactly the mid B
            if constexpr (sumsParticleCurrent)
                if (!info->midMagnetic.empty())
                {
                    auto const model
                        = std::find(std::begin(info->ghostMagnetic), std::end(info->ghostMagnetic),
                                    info->modelMagnetic);
                    assert(model != std::end(info->ghostMagnetic));
                    magGhostsMidRefiners_.addTimeRefiner(
                        info->modelMagnetic, info->midMagnetic, Bold_.name(), BRefineOp_,
                        vecFieldTimeOp_, info->modelMagnetic, nonOverwriteInteriorTFfillPattern,
                        magneticPatchStratPerGhostRefiner_[static_cast<std::size_t>(
                            std::distance(std::begin(info->ghostMagnetic), model))]);
                    midMagneticGhost_ = info->modelMagnetic;
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


            for (auto const& vecfield : info->ghostFlux)
                popFluxBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        sumVec_.name(), vecfield, nullptr, sumVec_.name(),
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());

            for (auto const& field : info->sumBorderFields)
                popDensityBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        sumField_.name(), field, nullptr, sumField_.name(),
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());


            if constexpr (sumsParticleCurrent)
            {
                currentBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        sumJ_.name(), info->modelCurrent, nullptr, sumJ_.name(),
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
                chargeDensityBorderSumRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        sumField_.name(), info->modelChargeDensity, nullptr, sumField_.name(),
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());
            }

            assert(info->maxBorderFields.size() == 1); // charge density
            for (auto const& field : info->maxBorderFields)
                ionDensityBorderMaxRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        field, field, nullptr, field,
                        std::make_shared<FieldGhostInterpOverlapFillPattern<GridLayoutT>>());

            assert(info->maxBorderVecFields.size() == 1); // charge flux
            for (auto const& vecfield : info->maxBorderVecFields)
                ionFluxBorderMaxRefiners_.emplace_back(resourcesManager_)
                    .addStaticRefiner(
                        vecfield, vecfield, nullptr, vecfield,
                        std::make_shared<
                            TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
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




        // the time of Eold_: the model E's own time, fully kinetic (E at half steps, P4 plan D8),
        // currentTime for the hybrid (E at integer times)
        double electricTime_(VecFieldT const& E, SAMRAI::hier::Patch const& patch,
                             double const currentTime) const
        {
            if constexpr (sumsParticleCurrent)
            {
                auto const times = resourcesManager_->getTimes(E, patch);
                assert(std::all_of(std::begin(times), std::end(times),
                                   [&](auto const t) { return t == times[0]; }));
                return times[0];
            }
            else
                return currentTime;
        }

        // the label of E on a level, the same on every patch: read on the local ones and agreed
        // over the ranks, one of which may hold no patch of the level (the fill is collective)
        double electricTimeOn_(VecFieldT const& E, level_t const& level) const
        {
            double local = std::numeric_limits<double>::max(), localMax = -local;
            for (auto const& patch : level)
            {
                // the third argument is the hybrid's currentTime, unused here (fully kinetic)
                auto const time = electricTime_(E, *patch, local);
                local           = std::min(local, time);
                localMax        = std::max(localMax, time);
            }
            auto const time = mpi::min(local);
            // all the patches of all the ranks must agree: min alone would hide a disagreement
            PHARE_DEBUG_DO({
                if (-mpi::min(-localMax) != time)
                    throw std::runtime_error("E labels differ on level "
                                             + std::to_string(level.getLevelNumber()));
            })
            return time;
        }

        // fully kinetic: a finer level E ghost fill past the coarse E^{n+1/2}, whose time is the
        // middle of the coarse step recorded by firstStep, interpolates towards the ahead E
        // (P4 plan D1). Inert while no ahead E is registered.
        bool fillsAhead_(int const levelNumber, double const fillTime) const
        {
            if (levelNumber == rootLevelNumber or elecGhostsAheadRefiners_.empty())
                return false;
            auto const coarseElectricTime
                = 0.5
                  * (beforePushCoarseTime_.at(levelNumber) + afterPushCoarseTime_.at(levelNumber));
            return fillTime > coarseElectricTime;
        }

        // fully kinetic: the model B ghost fill at the time of the last partial reflux onto the
        // coarser level (the start of the next substep) reads the coarse mid B. Every other fill,
        // the half step one of the push included, keeps the interpolation towards the coarse
        // B^{n+1}
        bool fillsMid_(VecFieldT const& B, int const levelNumber, double const fillTime) const
        {
            if (B.name() != midMagneticGhost_)
                return false;
            auto const it = midMagneticTime_.find(levelNumber);
            return it != midMagneticTime_.end() and it->second == fillTime;
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

        TensorFieldT sumTensor_{"PHARE_sumTensor", core::HybridQuantity::Tensor::M};
        VecFieldT sumVec_{"PHARE_sumVec", core::HybridQuantity::Vector::V};
        FieldT sumField_{"PHARE_sumField", core::HybridQuantity::Scalar::rho};
        VecFieldT sumJ_{"PHARE_sumJ", core::HybridQuantity::Vector::J}; // fully kinetic only



        //! ResourceManager shared with other objects (like the HybridModel)
        std::shared_ptr<ResourcesManagerT> resourcesManager_;
        std::shared_ptr<BoundaryManagerT> boundaryManager_;


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
        using FieldGhostMaxRefinerPool    = RefinerPool<rm_t, RefinerType::PatchFieldBorderMax>;
        using VecFieldGhostMaxRefinerPool = RefinerPool<rm_t, RefinerType::PatchVecFieldBorderMax>;
        using FieldFillPattern_t          = FieldFillPattern<dimension>;
        using TensorFieldFillPattern_t    = TensorFieldFillPattern<dimension /*, rank=1*/>;

        //! += flux on ghost box overlap incomplete population moment nodes
        std::vector<VecFieldGhostSumRefinerPool> popFluxBorderSumRefiners_;
        //! += density on ghost box overlap incomplete population moment nodes
        std::vector<FieldGhostSumRefinerPool> popDensityBorderSumRefiners_;

        //! += particle current (fully kinetic only)
        std::vector<VecFieldGhostSumRefinerPool> currentBorderSumRefiners_;
        //! += charge density of the divergence cleaning (fully kinetic only)
        std::vector<FieldGhostSumRefinerPool> chargeDensityBorderSumRefiners_;

        std::vector<FieldGhostMaxRefinerPool> ionDensityBorderMaxRefiners_;
        std::vector<VecFieldGhostMaxRefinerPool> ionFluxBorderMaxRefiners_;

        InitRefinerPool electricInitRefiners_{resourcesManager_};
        //! fully kinetic only: the E of a level (re)created during the run, and its level ghosts
        //! then (empty without info->aheadElectric)
        InitRefinerPool elecRegridRefiners_{resourcesManager_};
        GhostRefinerPool elecRegridGhostsRefiners_{resourcesManager_};


        SAMRAI::xfer::RefineAlgorithm BalgoPatchGhost;
        SAMRAI::xfer::RefineAlgorithm BalgoInit;
        SAMRAI::xfer::RefineAlgorithm BregridAlgo;
        SAMRAI::xfer::RefineAlgorithm EalgoPatchGhost;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> BpredGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> elecPatchGhostsRefineSchedules;

        SAMRAI::xfer::CoarsenAlgorithm RefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::RefineAlgorithm PatchGhostRefluxedAlgo;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> refluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> patchGhostRefluxedSchedules;

        //! fully kinetic partial reflux (B-x), registered when the solver names its fields
        bool partialReflux_ = false;
        SAMRAI::xfer::CoarsenAlgorithm PartialRefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::RefineAlgorithm PatchGhostPartialAlgo;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> partialRefluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> patchGhostPartialSchedules;
        //! per finer level, the time of the coarse mid B of the last partial reflux
        std::map<int, double> midMagneticTime_;
        std::string midMagneticGhost_;

        // Level objects each level's schedules were built against: {L} for same-level schedules,
        // [root, L] for inter-level ones (coarse interpolation can read below L-1). SAMRAI
        // replaces a reinitialised level object and never changes a kept one in place.
        std::map<int, std::vector<std::weak_ptr<level_t>>> sameLevelBuiltOn_;
        std::map<int, std::vector<std::weak_ptr<level_t>>> interLevelBuiltOn_;

        //! store refiners for electric fields that need ghosts to be filled
        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};
        //! fully kinetic only: finer level E ghosts past the coarse E, between the coarse E and
        //! the E one coarse step ahead (empty without info->aheadElectric)
        GhostRefinerPool elecGhostsAheadRefiners_{resourcesManager_};

        GhostRefinerPool magGhostsRefiners_{resourcesManager_};
        //! fully kinetic only: model B ghosts from the coarse mid B (empty without
        //! info->midMagnetic)
        GhostRefinerPool magGhostsMidRefiners_{resourcesManager_};

        //! fully kinetic only: J ghosts, after the border sum (current filter). Same level only:
        //! the level ghost band holds the level ghost particle deposit, never coarse values
        PatchGhostRefinerPool currentGhostsRefiners_{resourcesManager_};
        //! fully kinetic only: rho ghosts, after the border sum (filter of the cleaning rho), as J
        PatchGhostRefinerPool chargeDensityGhostsRefiners_{resourcesManager_};

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


        //! to grab particle leaving neighboring patches and inject into domain
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

        std::shared_ptr<TensorFieldFillPattern_t> nonOverwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>();

        std::shared_ptr<TensorFieldFillPattern<dimension, 1, /*overwrite_interior=*/true>>
            overwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension, 1, /*overwrite_interior=*/true>>();

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
