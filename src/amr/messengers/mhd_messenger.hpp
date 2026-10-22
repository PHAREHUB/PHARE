#ifndef PHARE_MHD_MESSENGER_HPP
#define PHARE_MHD_MESSENGER_HPP

#include "phare_mpi.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/models/quantities/mhd_quantities.hpp"

#include "amr/data/field/coarsening/electric_field_coarsener.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/mhd_flux_coarsener.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_composite_refiner.hpp"
#include "amr/messengers/refinement_config.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/messengers/refiner.hpp"
#include "amr/messengers/composite_refine_patch_strategy.hpp"
#include "amr/messengers/refiner_pool.hpp"
#include "amr/messengers/synchronizer_pool.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/data/field/refine/field_refine_patch_strategy.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"


#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/CoarsenOperator.h"

#include <memory>
#include <string>

namespace PHARE
{
namespace amr
{
    template<typename MHDModel>
    class MHDMessenger : public IMessenger<typename MHDModel::Interface>
    {
        using amr_types     = MHDModel::amr_types;
        using level_t       = amr_types::level_t;
        using patch_t       = amr_types::patch_t;
        using hierarchy_t   = amr_types::hierarchy_t;
        using level_ptr     = amr_types::level_ptr;
        using hierarchy_ptr = amr_types::hierarchy_ptr;

        using IPhysicalModel    = MHDModel::Interface;
        using FieldT            = MHDModel::field_type;
        using VecFieldT         = MHDModel::vecfield_type;
        using MHDStateT         = MHDModel::state_type;
        using GridLayoutT       = MHDModel::gridlayout_type;
        using GridT             = MHDModel::grid_type;
        using ResourcesManagerT = MHDModel::resources_manager_type;
        using BoundaryManagerT  = MHDModel::boundary_manager_type;
        using FieldDataT        = FieldData<GridLayoutT, GridT, core::MHDQuantity::Scalar>;
        using VectorFieldDataT  = TensorFieldData<1, GridLayoutT, GridT, core::MHDQuantity>;

        static constexpr auto dimension = MHDModel::dimension;

    public:
        static constexpr std::size_t rootLevelNumber = 0;
        static inline std::string const stratName    = "MHDModel-MHDModel";

        MHDMessenger(std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager,
                     std::shared_ptr<BoundaryManagerT> boundaryManager, int const firstLevel,
                     RefinementConfig const& refinementConfig = {})
            : resourcesManager_{std::move(resourcesManager)}
            , boundaryManager_{std::move(boundaryManager)}
            , firstLevel_{firstLevel}
        {
            makeRefineOperators_(refinementConfig);

            // moment ghosts are primitive quantities
            resourcesManager_->registerResources(rhoOld_);
            resourcesManager_->registerResources(Vold_);
            resourcesManager_->registerResources(Pold_);

            resourcesManager_->registerResources(rhoVold_);
            resourcesManager_->registerResources(EtotOld_);

            resourcesManager_->registerResources(Bold_);

            // also magnetic fluxes ? or should we use static refiners instead ?
        }

        virtual ~MHDMessenger() = default;

        void allocate(patch_t& patch, double const allocateTime) const override
        {
            resourcesManager_->allocate(rhoOld_, patch, allocateTime);
            resourcesManager_->allocate(Vold_, patch, allocateTime);
            resourcesManager_->allocate(Pold_, patch, allocateTime);

            resourcesManager_->allocate(rhoVold_, patch, allocateTime);
            resourcesManager_->allocate(EtotOld_, patch, allocateTime);

            resourcesManager_->allocate(Bold_, patch, allocateTime);
        }


        void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<MHDMessengerInfo> mhdInfo{
                dynamic_cast<MHDMessengerInfo*>(fromFinerInfo.release())};

            auto b_id = resourcesManager_->getID(mhdInfo->modelMagnetic);

            if (!b_id)
            {
                throw std::runtime_error(
                    "MHDMessengerStrategy: missing magnetic field variable IDs");
            }

            magneticRefinePatchStrategy_->registerIDs(*b_id);

            auto e_id = resourcesManager_->getID(mhdInfo->modelElectric);

            if (!e_id)
            {
                throw std::runtime_error(
                    "MHDMessengerStrategy: missing electric field variable IDs");
            }

            // EalgoPatchGhost.registerRefine(*e_id, *e_id, *e_id, EfieldRefineOp_,
            //                                nonOverwriteInteriorTFfillPattern);

            // refluxing
            // we first want to coarsen the flux sum onto the coarser level
            auto rho_fx_reflux_id  = resourcesManager_->getID(mhdInfo->reflux.rho_fx);
            auto rhoV_fx_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fx);
            auto Etot_fx_reflux_id = resourcesManager_->getID(mhdInfo->reflux.Etot_fx);

            if (!rho_fx_reflux_id or !rhoV_fx_reflux_id or !Etot_fx_reflux_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing reflux variable IDs for fluxes in x direction");
            }

            auto rho_fx_fluxsum_id  = resourcesManager_->getID(mhdInfo->fluxSum.rho_fx);
            auto rhoV_fx_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fx);
            auto Etot_fx_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.Etot_fx);


            if (!rho_fx_fluxsum_id or !rhoV_fx_fluxsum_id or !Etot_fx_fluxsum_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing flux sum variable IDs for fluxes in x direction");
            }


            // fluxes, flux sums and E of all directions share the refluxing algorithms: their
            // factory types differ by centering
            registerFluxReflux_(*rho_fx_reflux_id, *rhoV_fx_reflux_id, *Etot_fx_reflux_id,
                                *rho_fx_fluxsum_id, *rhoV_fx_fluxsum_id, *Etot_fx_fluxsum_id);

            if constexpr (dimension >= 2)
            {
                auto rho_fy_reflux_id  = resourcesManager_->getID(mhdInfo->reflux.rho_fy);
                auto rhoV_fy_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fy);
                auto Etot_fy_reflux_id = resourcesManager_->getID(mhdInfo->reflux.Etot_fy);

                if (!rho_fy_reflux_id or !rhoV_fy_reflux_id or !Etot_fy_reflux_id)
                {
                    throw std::runtime_error(
                        "MHDMessenger: missing reflux variable IDs for fluxes in y direction");
                }

                auto rho_fy_fluxsum_id  = resourcesManager_->getID(mhdInfo->fluxSum.rho_fy);
                auto rhoV_fy_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fy);
                auto Etot_fy_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.Etot_fy);

                if (!rho_fy_fluxsum_id or !rhoV_fy_fluxsum_id or !Etot_fy_fluxsum_id)
                {
                    throw std::runtime_error(
                        "MHDMessenger: missing flux sum variable IDs for fluxes in y direction");
                }

                registerFluxReflux_(*rho_fy_reflux_id, *rhoV_fy_reflux_id, *Etot_fy_reflux_id,
                                    *rho_fy_fluxsum_id, *rhoV_fy_fluxsum_id, *Etot_fy_fluxsum_id);

                if constexpr (dimension == 3)
                {
                    auto rho_fz_reflux_id  = resourcesManager_->getID(mhdInfo->reflux.rho_fz);
                    auto rhoV_fz_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fz);
                    auto Etot_fz_reflux_id = resourcesManager_->getID(mhdInfo->reflux.Etot_fz);


                    if (!rho_fz_reflux_id or !rhoV_fz_reflux_id or !Etot_fz_reflux_id)
                    {
                        throw std::runtime_error(
                            "MHDMessenger: missing reflux variable IDs for fluxes in z direction");
                    }

                    auto rho_fz_fluxsum_id  = resourcesManager_->getID(mhdInfo->fluxSum.rho_fz);
                    auto rhoV_fz_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fz);
                    auto Etot_fz_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.Etot_fz);

                    if (!rho_fz_fluxsum_id or !rhoV_fz_fluxsum_id or !Etot_fz_fluxsum_id)
                    {
                        throw std::runtime_error("MHDMessenger: missing flux sum variable IDs for "
                                                 "fluxes in z direction");
                    }

                    registerFluxReflux_(*rho_fz_reflux_id, *rhoV_fz_reflux_id, *Etot_fz_reflux_id,
                                        *rho_fz_fluxsum_id, *rhoV_fz_fluxsum_id,
                                        *Etot_fz_fluxsum_id);
                }
            }

            auto e_reflux_id = resourcesManager_->getID(mhdInfo->refluxElectric);

            auto e_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSumElectric);

            if (!e_reflux_id or !e_fluxsum_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing electric refluxing field variable IDs");
            }

            refluxAlgo_.registerCoarsen(*e_reflux_id, *e_fluxsum_id, electricFieldCoarseningOp_);

            patchGhostRefluxedAlgo_.registerRefine(*e_reflux_id, *e_reflux_id, *e_reflux_id,
                                                   EfieldRefineOp_,
                                                   nonOverwriteInteriorTFfillPattern);

            registerGhostComms_(mhdInfo);
            registerInitComms_(mhdInfo);
        }



        void registerLevel(hierarchy_ptr const& hierarchy, int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            // magPatchGhostsRefineSchedules[levelNumber]
            //     = BalgoPatchGhost.createSchedule(level, &magneticRefinePatchStrategy_);

            // elecPatchGhostsRefineSchedules[levelNumber] = EalgoPatchGhost.createSchedule(level);

            auto& patchGhostRefluxed = patchGhostRefluxedSchedules_[levelNumber];
            patchGhostRefluxed       = patchGhostRefluxedAlgo_.createSchedule(level);
            checkOverwriteInteriorAgrees(*patchGhostRefluxed->getEquivalenceClasses());

            elecGhostsRefiners_.registerLevel(hierarchy, level);

            stateGhostsRefiners_.registerLevel(hierarchy, level);

            magFluxesXGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesYGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesZGhostRefiners_.registerLevel(hierarchy, level);

            magMaxRefiners_.registerLevel(hierarchy, level);

            if (levelNumber != rootLevelNumber)
            {
                // refluxing
                auto const& coarseLevel       = hierarchy->getPatchLevel(levelNumber - 1);
                refluxSchedules_[levelNumber] = refluxAlgo_.createSchedule(coarseLevel, level);

                // refinement
                stateInitRefiners_.registerLevel(hierarchy, level);
            }
        }


        void regrid(hierarchy_ptr const& hierarchy, int const levelNumber,
                    level_ptr const& oldLevel, IPhysicalModel& model,
                    double const initDataTime) override
        {
            auto& mhdModel = static_cast<MHDModel&>(model);
            auto level     = hierarchy->getPatchLevel(levelNumber);

            bool isRegriddingL0 = levelNumber == 0 and oldLevel;

            stateInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            magMaxRefiners_.fill(mhdModel.state.B, level->getLevelNumber(), initDataTime);

            // magPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);
            // elecPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);
        }


        std::string fineModelName() const override { return MHDModel::model_name; }

        std::string coarseModelName() const override { return MHDModel::model_name; }

        std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }

        void initLevel(IPhysicalModel& model, level_t& level, double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();

            auto& mhdModel = static_cast<MHDModel&>(model);

            stateInitRefiners_.fill(levelNumber, initDataTime);
        }

        void firstStep(IPhysicalModel& model, level_t& level, double const currentTime,
                       double const prevCoarserTIme, double const newCoarserTime) final
        {
        }


        void lastStep(IPhysicalModel& model, level_t& level) final {}


        void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) final
        {
            auto& mhdModel = static_cast<MHDModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, mhdModel.state.rho, mhdModel.state.V, mhdModel.state.P,
                    mhdModel.state.rhoV, mhdModel.state.Etot, mhdModel.state.B, rhoOld_, Vold_,
                    Pold_, rhoVold_, EtotOld_, Bold_);

                resourcesManager_->setTime(rhoOld_, *patch, currentTime);
                resourcesManager_->setTime(Vold_, *patch, currentTime);
                resourcesManager_->setTime(Pold_, *patch, currentTime);
                resourcesManager_->setTime(rhoVold_, *patch, currentTime);
                resourcesManager_->setTime(EtotOld_, *patch, currentTime);
                resourcesManager_->setTime(Bold_, *patch, currentTime);

                rhoOld_.copyData(mhdModel.state.rho);
                Vold_.copyData(mhdModel.state.V);
                Pold_.copyData(mhdModel.state.P);
                rhoVold_.copyData(mhdModel.state.rhoV);
                EtotOld_.copyData(mhdModel.state.Etot);
                Bold_.copyData(mhdModel.state.B);
            }
        }

        void fillRootGhosts(IPhysicalModel& model, level_t& level, double const initDataTime) final
        {
            auto& mhdModel = static_cast<MHDModel&>(model);
            fillStateGhosts(mhdModel.state, level, initDataTime);
        }

        void synchronize(level_t& level) final {}

        void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                    double const syncTime) override
        {
            refluxSchedules_[fineLevelNumber]->coarsenData();
            patchGhostRefluxedSchedules_[coarserLevelNumber]->fillData(syncTime);
        }

        void postSynchronize(IPhysicalModel& model, level_t& level, double const time) override
        {
            // The ghosts for B are obtained in the solver's reflux_euler. For B, this is because
            // refluxing is done through faraday which is computed on the ghost box for the other
            // quantities, the ghosts are filled in the end of the euler step anyways.
        }

        // Fills the ghosts of the conserved state (B, rho, rhoV, Etot), then equalizes the B
        // values shared by neighbour patches.
        void fillStateGhosts(MHDStateT& state, level_t const& level, double const fillTime)
        {
            PHARE_LOG_SCOPE(3, "MHDMessenger::fillStateGhosts");

            auto const levelNumber = level.getLevelNumber();

            setNaNsOnVecfieldGhosts(state.B, level);
            setNaNsOnFieldGhosts(state.rho, level);
            setNaNsOnVecfieldGhosts(state.rhoV, level);
            setNaNsOnFieldGhosts(state.Etot, level);
            stateGhostsRefiners_.fillByKey(state.B.name(), levelNumber, fillTime);
            magMaxRefiners_.fill(state.B, levelNumber, fillTime);
        }

        void fillMagneticFluxesXGhosts(VecFieldT& Fx_B, level_t const& level, double const fillTime)
        {
            setNaNsOnVecfieldGhosts(Fx_B, level);
            magFluxesXGhostRefiners_.fill(Fx_B, level.getLevelNumber(), fillTime);
        }

        void fillMagneticFluxesYGhosts(VecFieldT& Fy_B, level_t const& level, double const fillTime)
        {
            setNaNsOnVecfieldGhosts(Fy_B, level);
            magFluxesYGhostRefiners_.fill(Fy_B, level.getLevelNumber(), fillTime);
        }

        void fillMagneticFluxesZGhosts(VecFieldT& Fz_B, level_t const& level, double const fillTime)
        {
            setNaNsOnVecfieldGhosts(Fz_B, level);
            magFluxesZGhostRefiners_.fill(Fz_B, level.getLevelNumber(), fillTime);
        }

        void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime)
        {
            setNaNsOnVecfieldGhosts(E, level);
            elecGhostsRefiners_.fill(E, level.getLevelNumber(), fillTime);
        }

        std::string name() override { return stratName; }



    private:
        // Select the field-refinement operators once at construction. The composite runtime
        // kernels are built from the configured order.
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

            mhdFluxRefineOp_     = fieldKernel();
            mhdVecFluxRefineOp_  = vecKernel();
            mhdFieldRefineOp_    = fieldKernel();
            mhdVecFieldRefineOp_ = vecKernel();
            EfieldRefineOp_      = vecKernel();
            BfieldRegridOp_      = magKernel();
        }

        // Maybe we also need conservative ghost refiners for amr operations, actually quite
        // likely
        void registerGhostComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            // static refinement for E because in MHD it is a temporary, so keeping its
            // state updated after each regrid is not a priority. However if we do not correctly
            // refine on regrid, the post regrid state is not up to date (in our case it will be nan
            // since we nan-initialise) and thus is is better to rely on static refinement, which
            // uses the state after computation of CT.
            // The refiners for the electric field only serve for filling ghosts at physical
            // boundaries.
            registerGhostRefinePatchStrategies_(elecPatchStrats, info->ghostElectric);
            for (std::size_t i = 0; i < info->ghostElectric.size(); ++i)
                elecGhostsRefiners_.addStaticRefiner(info->ghostElectric[i], nullptr,
                                                     info->ghostElectric[i], noOverlapFillPattern,
                                                     elecPatchStrats[i]);

            // each quantity keeps its own patch strategy, so that the physical-boundary ghosts
            // are filled by its boundary condition; a state's strategies are combined below
            registerGhostRefinePatchStrategies_(rhoPatchStrats, info->ghostDensity);
            registerGhostRefinePatchStrategies_(momentumPatchStrats, info->ghostMomentum);
            registerGhostRefinePatchStrategies_(totalEnergyPatchStrats, info->ghostTotalEnergy);

            magFluxesXGhostRefiners_.addStaticRefiners(
                info->ghostMagneticFluxesX, mhdVecFluxRefineOp_, info->ghostMagneticFluxesX,
                nonOverwriteInteriorTFfillPattern);

            magFluxesYGhostRefiners_.addStaticRefiners(
                info->ghostMagneticFluxesY, mhdVecFluxRefineOp_, info->ghostMagneticFluxesY,
                nonOverwriteInteriorTFfillPattern);

            magFluxesZGhostRefiners_.addStaticRefiners(
                info->ghostMagneticFluxesZ, mhdVecFluxRefineOp_, info->ghostMagneticFluxesZ,
                nonOverwriteInteriorTFfillPattern);

            // we need a separate patch strategy for each refiner so that each one can register
            // their required ids
            registerGhostRefinePatchStrategies_(magPatchStrats, info->ghostMagnetic);

            // the i-th entry of each ghost list belongs to the same state
            auto const nStates = info->ghostMagnetic.size();
            if (info->ghostDensity.size() != nStates or info->ghostMomentum.size() != nStates
                or info->ghostTotalEnergy.size() != nStates)
                throw std::runtime_error("MHDMessenger: ghost state lists differ in size");

            for (std::size_t i = 0; i < nStates; ++i)
            {
                stateGhostsRefiners_.addRefiner(
                    info->ghostMagnetic[i],
                    {{.dst         = info->ghostMagnetic[i],
                      .src         = info->ghostMagnetic[i],
                      .told        = Bold_.name(),
                      .tnew        = info->modelMagnetic,
                      .refineOp    = BfieldRegridOp_,
                      .timeOp      = vecFieldTimeOp_,
                      .fillPattern = nonOverwriteInteriorTFfillPattern},
                     {.dst         = info->ghostDensity[i],
                      .src         = info->ghostDensity[i],
                      .told        = rhoOld_.name(),
                      .tnew        = info->modelDensity,
                      .refineOp    = mhdFieldRefineOp_,
                      .timeOp      = fieldTimeOp_,
                      .fillPattern = nonOverwriteFieldFillPattern},
                     {.dst         = info->ghostMomentum[i],
                      .src         = info->ghostMomentum[i],
                      .told        = rhoVold_.name(),
                      .tnew        = info->modelMomentum,
                      .refineOp    = mhdVecFieldRefineOp_,
                      .timeOp      = vecFieldTimeOp_,
                      .fillPattern = nonOverwriteInteriorTFfillPattern},
                     {.dst         = info->ghostTotalEnergy[i],
                      .src         = info->ghostTotalEnergy[i],
                      .told        = EtotOld_.name(),
                      .tnew        = info->modelTotalEnergy,
                      .refineOp    = mhdFieldRefineOp_,
                      .timeOp      = fieldTimeOp_,
                      .fillPattern = nonOverwriteFieldFillPattern}},
                    std::make_shared<CompositeRefinePatchStrategy>(
                        std::vector<CompositeRefinePatchStrategy::strategy_ptr>{
                            magPatchStrats[i], rhoPatchStrats[i], momentumPatchStrats[i],
                            totalEnergyPatchStrats[i]}));

                magMaxRefiners_.addStaticRefiner(
                    info->ghostMagnetic[i], info->ghostMagnetic[i], nullptr, info->ghostMagnetic[i],
                    std::make_shared<
                        TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
            }
        }




        // refluxes the fluxes of one direction, then refills their ghosts so that they agree with
        // the newly refluxed cells
        void registerFluxReflux_(int const rho, int const rhoV, int const Etot, int const rhoSum,
                                 int const rhoVSum, int const EtotSum)
        {
            refluxAlgo_.registerCoarsen(rho, rhoSum, mhdFluxCoarseningOp_);
            refluxAlgo_.registerCoarsen(rhoV, rhoVSum, mhdVecFluxCoarseningOp_);
            refluxAlgo_.registerCoarsen(Etot, EtotSum, mhdFluxCoarseningOp_);

            patchGhostRefluxedAlgo_.registerRefine(rho, rho, rho, mhdFluxRefineOp_,
                                                   nonOverwriteInteriorTFfillPattern);
            patchGhostRefluxedAlgo_.registerRefine(rhoV, rhoV, rhoV, mhdVecFluxRefineOp_,
                                                   nonOverwriteInteriorTFfillPattern);
            patchGhostRefluxedAlgo_.registerRefine(Etot, Etot, Etot, mhdFluxRefineOp_,
                                                   nonOverwriteInteriorTFfillPattern);
        }




        /**
         * @brief Register a list of refine patch strategy pointers corresponding to a list of
         * keys.
         *
         * @tparam RefinePatchStrategyT type inheriting from SAMRAI's `RefinePatchStrategy`
         * @param patchStrategies the list of refine patch strategy pointers.
         * @param keys the list of keys.
         */
        template<typename RefinePatchStrategyT>
        void registerGhostRefinePatchStrategies_(
            std::vector<std::shared_ptr<RefinePatchStrategyT>>& patchStrategies,
            std::vector<std::string> const& keys)
        {
            patchStrategies.clear(); // registerQuantities must be idempotent
            patchStrategies.reserve(keys.size());
            for (auto const& key : keys)
            {
                auto&& [id]     = resourcesManager_->getIDsList(key);
                auto patchStrat = std::make_shared<RefinePatchStrategyT>(*boundaryManager_);
                patchStrat->registerIDs(id);
                patchStrategies.push_back(patchStrat);
            }
        }


        // should this use conservative quantities ? When should we do the initial conversion ?
        // Maybe mhd_init
        void registerInitComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            // The model state's moment strategies (index 0 of each ghost-strategy list:
            // ghostX[0] == modelX == initX) fill the physical-boundary ghosts of a new or
            // regridded level, so that it carries valid moment ghosts before its first flux.
            // Moments keep the default (non-overwrite) fill pattern: overwrite is a B/face-centered
            // concern and corrupts the cell-centered moment interior fill.
            std::vector<CompositeRefinePatchStrategy::strategy_ptr> strategies{
                magneticRefinePatchStrategy_};
            if (!rhoPatchStrats.empty())
                strategies.push_back(rhoPatchStrats.front());
            if (!momentumPatchStrats.empty())
                strategies.push_back(momentumPatchStrats.front());
            if (!totalEnergyPatchStrats.empty())
                strategies.push_back(totalEnergyPatchStrats.front());

            std::vector<RefineItem> items{{.dst         = info->modelMagnetic,
                                           .src         = info->modelMagnetic,
                                           .refineOp    = BfieldRegridOp_,
                                           .fillPattern = overwriteInteriorTFfillPattern}};
            for (auto const& rho : info->initDensity)
                items.push_back({.dst = rho, .src = rho, .refineOp = mhdFieldRefineOp_});
            for (auto const& rhoV : info->initMomentum)
                items.push_back({.dst = rhoV, .src = rhoV, .refineOp = mhdVecFieldRefineOp_});
            for (auto const& Etot : info->initTotalEnergy)
                items.push_back({.dst = Etot, .src = Etot, .refineOp = mhdFieldRefineOp_});

            stateInitRefiners_.addRefiner(
                info->modelMagnetic, items,
                std::make_shared<CompositeRefinePatchStrategy>(std::move(strategies)));
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
                for (auto& component : vf)
                    setNaNsOnFieldGhosts(component);
        }


        FieldT rhoOld_{stratName + "rhoOld", core::MHDQuantity::Scalar::rho};
        VecFieldT Vold_{stratName + "Vold", core::MHDQuantity::Vector::V};
        FieldT Pold_{stratName + "Pold", core::MHDQuantity::Scalar::P};

        VecFieldT rhoVold_{stratName + "rhoVold", core::MHDQuantity::Vector::rhoV};
        FieldT EtotOld_{stratName + "EtotOld", core::MHDQuantity::Scalar::Etot};

        VecFieldT Bold_{stratName + "Bold", core::MHDQuantity::Vector::B};


        using rm_t = typename MHDModel::resources_manager_type;
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        std::shared_ptr<BoundaryManagerT> boundaryManager_;
        int const firstLevel_;

        using InitRefinerPool             = RefinerPool<rm_t, RefinerType::InitField>;
        using GhostRefinerPool            = RefinerPool<rm_t, RefinerType::GhostField>;
        using PatchGhostRefinerPool       = RefinerPool<rm_t, RefinerType::PatchGhostField>;
        using InitDomPartRefinerPool      = RefinerPool<rm_t, RefinerType::InitInteriorPart>;
        using VecFieldGhostMaxRefinerPool = RefinerPool<rm_t, RefinerType::PatchVecFieldBorderMax>;


        SAMRAI::xfer::RefineAlgorithm BalgoPatchGhost;                                         //
        SAMRAI::xfer::RefineAlgorithm EalgoPatchGhost;                                         //
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules; //
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            magPatchGhostsRefineSchedules; //
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> elecPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            magSharedNodeRefineSchedules; //

        // coarsens the flux sums and E onto the coarser level, then refills their patch ghosts
        SAMRAI::xfer::CoarsenAlgorithm refluxAlgo_{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::RefineAlgorithm patchGhostRefluxedAlgo_;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> refluxSchedules_;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> patchGhostRefluxedSchedules_;

        PatchGhostRefinerPool elecGhostsRefiners_{resourcesManager_};
        // B, rho, rhoV and Etot of each state, keyed by the state's B name
        GhostRefinerPool stateGhostsRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesXGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesYGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesZGhostRefiners_{resourcesManager_};

        // ghostMagnetic includes the model B, also equalized after regrid
        VecFieldGhostMaxRefinerPool magMaxRefiners_{resourcesManager_};

        // B, rho, rhoV and Etot of the model state on new and regridded levels
        InitRefinerPool stateInitRefiners_{resourcesManager_};

        // SynchronizerPool<rm_t> densitySynchronizers_{resourcesManager_};
        // SynchronizerPool<rm_t> momentumSynchronizers_{resourcesManager_};
        // SynchronizerPool<rm_t> magnetoSynchronizers_{resourcesManager_};
        // SynchronizerPool<rm_t> totalEnergySynchronizers_{resourcesManager_};

        using RefOp_ptr     = std::shared_ptr<SAMRAI::hier::RefineOperator>;
        using CoarsenOp_ptr = std::shared_ptr<SAMRAI::hier::CoarsenOperator>;
        using TimeOp_ptr    = std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>;

        using FieldTimeInterp = FieldLinearTimeInterpolate<GridLayoutT, GridT>;

        using VecFieldTimeInterp
            = VecFieldLinearTimeInterpolate<GridLayoutT, GridT, core::MHDQuantity>;

        template<typename Policy>
        using FieldCoarseningOp = FieldCoarsenOperator<GridLayoutT, GridT, Policy>;

        template<typename Policy>
        using VecFieldCoarsenOp
            = VecFieldCoarsenOperator<GridLayoutT, GridT, Policy, core::MHDQuantity>;

        using MHDFluxCoarsenOp       = FieldCoarseningOp<MHDFluxCoarsener<dimension>>;
        using MHDVecFluxCoarsenOp    = VecFieldCoarsenOp<MHDFluxCoarsener<dimension>>;
        using ElectricFieldCoarsenOp = VecFieldCoarsenOp<ElectricFieldCoarsener<dimension>>;

        SynchronizerPool<rm_t> electroSynchronizers_{resourcesManager_};

        // built in the ctor body (makeRefineOperators_) from the composite runtime kernels.
        RefOp_ptr mhdFluxRefineOp_;
        RefOp_ptr mhdVecFluxRefineOp_;
        RefOp_ptr mhdFieldRefineOp_;
        RefOp_ptr mhdVecFieldRefineOp_;
        RefOp_ptr EfieldRefineOp_;
        RefOp_ptr BfieldRegridOp_;

        TimeOp_ptr fieldTimeOp_{std::make_shared<FieldTimeInterp>()};
        TimeOp_ptr vecFieldTimeOp_{std::make_shared<VecFieldTimeInterp>()};

        using TensorFieldFillPattern_t = TensorFieldFillPattern<dimension /*, rank=1*/>;
        using FieldFillPattern_t       = FieldFillPattern<dimension>;

        std::shared_ptr<FieldFillPattern_t> nonOverwriteFieldFillPattern
            = std::make_shared<FieldFillPattern<dimension>>(); // stateless (mostly)

        std::shared_ptr<TensorFieldFillPattern_t> nonOverwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>();

        std::shared_ptr<NoOverlapFillPattern<dimension>> noOverlapFillPattern
            = std::make_shared<NoOverlapFillPattern<dimension>>();

        std::shared_ptr<TensorFieldFillPattern_t> overwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>(
                /*overwrite_interior=*/true);

        CoarsenOp_ptr mhdFluxCoarseningOp_{std::make_shared<MHDFluxCoarsenOp>()};
        CoarsenOp_ptr mhdVecFluxCoarseningOp_{std::make_shared<MHDVecFluxCoarsenOp>()};
        CoarsenOp_ptr electricFieldCoarseningOp_{std::make_shared<ElectricFieldCoarsenOp>()};

        using FieldRefinePatchStrategyT = FieldRefinePatchStrategy<FieldDataT, BoundaryManagerT>;
        using VectorFieldRefinePatchStrategyT
            = FieldRefinePatchStrategy<VectorFieldDataT, BoundaryManagerT>;
        using MagneticRefinePatchStrategyT
            = MagneticRefinePatchStrategy<VectorFieldDataT, BoundaryManagerT>;
        using FieldRefinePatchStrategyList
            = std::vector<std::shared_ptr<FieldRefinePatchStrategyT>>;
        using VectorFieldRefinePatchStrategyList
            = std::vector<std::shared_ptr<VectorFieldRefinePatchStrategyT>>;
        using MagneticRefinePatchStrategyList
            = std::vector<std::shared_ptr<MagneticRefinePatchStrategyT>>;

        std::shared_ptr<MagneticRefinePatchStrategyT> magneticRefinePatchStrategy_
            = std::make_shared<MagneticRefinePatchStrategyT>(*boundaryManager_);

        FieldRefinePatchStrategyList rhoPatchStrats;
        FieldRefinePatchStrategyList totalEnergyPatchStrats;
        VectorFieldRefinePatchStrategyList momentumPatchStrats;
        VectorFieldRefinePatchStrategyList elecPatchStrats;
        MagneticRefinePatchStrategyList magPatchStrats;
    };

} // namespace amr
} // namespace PHARE
#endif
