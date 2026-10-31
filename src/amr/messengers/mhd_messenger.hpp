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
#include "amr/messengers/refiner_pool.hpp"
#include "amr/messengers/synchronizer_pool.hpp"
#include "amr/messengers/field_operate_transaction.hpp"
#include "amr/physical_models/models.hpp"
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

        // schedules and temporaries are never tiled, border max works on an untiled temporary
        using Grid_rt            = MHDModel::storage_t::grid_type;
        using VectorFieldData_rt = TensorFieldData<1, GridLayoutT, Grid_rt, core::MHDQuantity>;

        using BorderMaxOp = core::FieldBorderMaxOp<typename VecFieldT::value_type>;
        using VecFieldBorderMaxTransactionFactory_t
            = FieldBorderOpTransactionFactory<VectorFieldData_rt, BorderMaxOp>;

        static constexpr auto dimension = MHDModel::dimension;

    public:
        static constexpr std::size_t rootLevelNumber = 0;
        static inline std::string const stratName    = "MHDModel-MHDModel";

        MHDMessenger(MHDModel& model, int const firstLevel,
                     RefinementConfig const& refinementConfig = {})
            : model_{model}
            , resourcesManager_{model.resourcesManager}
            , boundaryManager_{model.boundaryManager}
            , firstLevel_{firstLevel}
        {
            makeRefineOperators_(refinementConfig);

            // scratch for tiled refines, see scratchID_
            model_.requestTemporaries(
                core::QtyCentering::primal, core::MHDQuantity::Vector::VecAllPrimal,
                core::MHDQuantity::Vector::B); // B border max, see borderMaxB_

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

            magneticRefinePatchStrategy_.registerIDs(*b_id);
            magneticRefinePatchStrategy_.registerScratchID(
                scratchID_(model_.getTmp(core::MHDQuantity::Vector::VecAllPrimal)));

            registerRefineScratchID(scratchID_(model_.getTmp(core::QtyCentering::primal)),
                                    mhdFluxRefineOp_, mhdVecFluxRefineOp_, mhdFieldRefineOp_,
                                    mhdVecFieldRefineOp_, EfieldRefineOp_, BfieldRegridOp_);

            BalgoInit.registerRefine(*b_id, *b_id, *b_id, BfieldRegridOp_,
                                     overwriteInteriorTFfillPattern);

            BregridAlgo.registerRefine(*b_id, *b_id, *b_id, BfieldRegridOp_,
                                       overwriteInteriorTFfillPattern);

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


            // all of the fluxes fx are defined on the same faces no matter the component, so we
            // just need a different fill pattern per direction
            HydroXrefluxAlgo.registerCoarsen(*rho_fx_reflux_id, *rho_fx_fluxsum_id,
                                             mhdFluxCoarseningOp_);
            HydroXrefluxAlgo.registerCoarsen(*rhoV_fx_reflux_id, *rhoV_fx_fluxsum_id,
                                             mhdVecFluxCoarseningOp_);
            HydroXrefluxAlgo.registerCoarsen(*Etot_fx_reflux_id, *Etot_fx_fluxsum_id,
                                             mhdFluxCoarseningOp_);

            // we then need to refill the ghosts so that they agree with the newly refluxed
            // cells
            HydroXpatchGhostRefluxedAlgo.registerRefine(*rho_fx_reflux_id, *rho_fx_reflux_id,
                                                        *rho_fx_reflux_id, mhdFluxRefineOp_,
                                                        nonOverwriteInteriorTFfillPattern);
            HydroXpatchGhostRefluxedAlgo.registerRefine(*rhoV_fx_reflux_id, *rhoV_fx_reflux_id,
                                                        *rhoV_fx_reflux_id, mhdVecFluxRefineOp_,
                                                        nonOverwriteInteriorTFfillPattern);
            HydroXpatchGhostRefluxedAlgo.registerRefine(*Etot_fx_reflux_id, *Etot_fx_reflux_id,
                                                        *Etot_fx_reflux_id, mhdFluxRefineOp_,
                                                        nonOverwriteInteriorTFfillPattern);

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

                HydroYrefluxAlgo.registerCoarsen(*rho_fy_reflux_id, *rho_fy_fluxsum_id,
                                                 mhdFluxCoarseningOp_);
                HydroYrefluxAlgo.registerCoarsen(*rhoV_fy_reflux_id, *rhoV_fy_fluxsum_id,
                                                 mhdVecFluxCoarseningOp_);
                HydroYrefluxAlgo.registerCoarsen(*Etot_fy_reflux_id, *Etot_fy_fluxsum_id,
                                                 mhdFluxCoarseningOp_);

                HydroYpatchGhostRefluxedAlgo.registerRefine(*rho_fy_reflux_id, *rho_fy_reflux_id,
                                                            *rho_fy_reflux_id, mhdFluxRefineOp_,
                                                            nonOverwriteInteriorTFfillPattern);
                HydroYpatchGhostRefluxedAlgo.registerRefine(*rhoV_fy_reflux_id, *rhoV_fy_reflux_id,
                                                            *rhoV_fy_reflux_id, mhdVecFluxRefineOp_,
                                                            nonOverwriteInteriorTFfillPattern);
                HydroYpatchGhostRefluxedAlgo.registerRefine(*Etot_fy_reflux_id, *Etot_fy_reflux_id,
                                                            *Etot_fy_reflux_id, mhdFluxRefineOp_,
                                                            nonOverwriteInteriorTFfillPattern);

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

                    HydroZrefluxAlgo.registerCoarsen(*rho_fz_reflux_id, *rho_fz_fluxsum_id,
                                                     mhdFluxCoarseningOp_);
                    HydroZrefluxAlgo.registerCoarsen(*rhoV_fz_reflux_id, *rhoV_fz_fluxsum_id,
                                                     mhdVecFluxCoarseningOp_);
                    HydroZrefluxAlgo.registerCoarsen(*Etot_fz_reflux_id, *Etot_fz_fluxsum_id,
                                                     mhdFluxCoarseningOp_);


                    HydroZpatchGhostRefluxedAlgo.registerRefine(
                        *rho_fz_reflux_id, *rho_fz_reflux_id, *rho_fz_reflux_id, mhdFluxRefineOp_,
                        nonOverwriteInteriorTFfillPattern);
                    HydroZpatchGhostRefluxedAlgo.registerRefine(
                        *rhoV_fz_reflux_id, *rhoV_fz_reflux_id, *rhoV_fz_reflux_id,
                        mhdVecFluxRefineOp_, nonOverwriteInteriorTFfillPattern);
                    HydroZpatchGhostRefluxedAlgo.registerRefine(
                        *Etot_fz_reflux_id, *Etot_fz_reflux_id, *Etot_fz_reflux_id,
                        mhdFluxRefineOp_, nonOverwriteInteriorTFfillPattern);
                }
            }

            auto e_reflux_id = resourcesManager_->getID(mhdInfo->refluxElectric);

            auto e_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSumElectric);

            if (!e_reflux_id or !e_fluxsum_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing electric refluxing field variable IDs");
            }

            ErefluxAlgo.registerCoarsen(*e_reflux_id, *e_fluxsum_id, electricFieldCoarseningOp_);

            EpatchGhostRefluxedAlgo.registerRefine(*e_reflux_id, *e_reflux_id, *e_reflux_id,
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

            EpatchGhostRefluxedSchedules[levelNumber]
                = EpatchGhostRefluxedAlgo.createSchedule(level);
            HydroXpatchGhostRefluxedSchedules[levelNumber]
                = HydroXpatchGhostRefluxedAlgo.createSchedule(level);
            HydroYpatchGhostRefluxedSchedules[levelNumber]
                = HydroYpatchGhostRefluxedAlgo.createSchedule(level);
            HydroZpatchGhostRefluxedSchedules[levelNumber]
                = HydroZpatchGhostRefluxedAlgo.createSchedule(level);

            elecGhostsRefiners_.registerLevel(hierarchy, level);

            rhoGhostsRefiners_.registerLevel(hierarchy, level);
            // velGhostsRefiners_.registerLevel(hierarchy, level);
            // pressureGhostsRefiners_.registerLevel(hierarchy, level);

            momentumGhostsRefiners_.registerLevel(hierarchy, level);
            totalEnergyGhostsRefiners_.registerLevel(hierarchy, level);

            magFluxesXGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesYGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesZGhostRefiners_.registerLevel(hierarchy, level);

            magGhostsRefiners_.registerLevel(hierarchy, level);

            {
                auto border_max_factory = std::make_shared<VecFieldBorderMaxTransactionFactory_t>();
                magMaxRefiners_.registerLevel(hierarchy, level, border_max_factory);
            }

            if (levelNumber != rootLevelNumber)
            {
                // refluxing
                auto const& coarseLevel       = hierarchy->getPatchLevel(levelNumber - 1);
                ErefluxSchedules[levelNumber] = ErefluxAlgo.createSchedule(coarseLevel, level);
                HydroXrefluxSchedules[levelNumber]
                    = HydroXrefluxAlgo.createSchedule(coarseLevel, level);
                HydroYrefluxSchedules[levelNumber]
                    = HydroYrefluxAlgo.createSchedule(coarseLevel, level);
                HydroZrefluxSchedules[levelNumber]
                    = HydroZrefluxAlgo.createSchedule(coarseLevel, level);

                // refinement
                magInitRefineSchedules[levelNumber] = BalgoInit.createSchedule(
                    level, nullptr, levelNumber - 1, hierarchy, &magneticRefinePatchStrategy_);

                densityInitRefiners_.registerLevel(hierarchy, level);
                momentumInitRefiners_.registerLevel(hierarchy, level);
                totalEnergyInitRefiners_.registerLevel(hierarchy, level);
            }
        }


        void regrid(hierarchy_ptr const& hierarchy, int const levelNumber,
                    level_ptr const& oldLevel, IPhysicalModel& model,
                    double const initDataTime) override
        {
            auto& mhdModel = static_cast<MHDModel&>(model);
            auto level     = hierarchy->getPatchLevel(levelNumber);

            bool isRegriddingL0 = levelNumber == 0 and oldLevel;

            magneticRegriding_(hierarchy, level, oldLevel, initDataTime);
            borderMaxB_(mhdModel.state.B, *level, initDataTime);

            densityInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            momentumInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            totalEnergyInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);

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

            magInitRefineSchedules[levelNumber]->fillData(initDataTime);
            densityInitRefiners_.fill(levelNumber, initDataTime);
            momentumInitRefiners_.fill(levelNumber, initDataTime);
            totalEnergyInitRefiners_.fill(levelNumber, initDataTime);
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
            fillMagneticGhosts(mhdModel.state.B, level, initDataTime);
            fillMomentsGhosts(mhdModel.state, level, initDataTime);
        }

        void synchronize(level_t& level) final {}

        void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                    double const syncTime) override
        {
            ErefluxSchedules[fineLevelNumber]->coarsenData();
            HydroXrefluxSchedules[fineLevelNumber]->coarsenData();
            HydroYrefluxSchedules[fineLevelNumber]->coarsenData();
            HydroZrefluxSchedules[fineLevelNumber]->coarsenData();

            EpatchGhostRefluxedSchedules[coarserLevelNumber]->fillData(syncTime);
            HydroXpatchGhostRefluxedSchedules[coarserLevelNumber]->fillData(syncTime);
            HydroYpatchGhostRefluxedSchedules[coarserLevelNumber]->fillData(syncTime);
            HydroZpatchGhostRefluxedSchedules[coarserLevelNumber]->fillData(syncTime);
        }

        void postSynchronize(IPhysicalModel& model, level_t& level, double const time) override
        {
            // The ghosts for B are obtained in the solver's reflux_euler. For B, this is because
            // refluxing is done through faraday which is computed on the ghost box for the other
            // quantities, the ghosts are filled in the end of the euler step anyways.
        }

        void fillMomentsGhosts(MHDStateT& state, level_t const& level, double const fillTime)
        {
            setNaNsOnFieldGhosts(state.rho, level);
            setNaNsOnVecfieldGhosts(state.rhoV, level);
            setNaNsOnFieldGhosts(state.Etot, level);
            rhoGhostsRefiners_.fill(state.rho, level.getLevelNumber(), fillTime);
            momentumGhostsRefiners_.fill(state.rhoV, level.getLevelNumber(), fillTime);
            totalEnergyGhostsRefiners_.fill(state.Etot, level.getLevelNumber(), fillTime);
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

        void fillMagneticGhosts(VecFieldT& B, level_t const& level, double const fillTime)
        {
            PHARE_LOG_SCOPE(3, "MHDMessenger::fillMagneticGhosts");

            setNaNsOnVecfieldGhosts(B, level);
            magGhostsRefiners_.fill(B, level.getLevelNumber(), fillTime);
            borderMaxB_(B, level, fillTime);
        }

        std::string name() override { return stratName; }



    private:
        // tiled field refines borrow the untiled scratch fields of the destination patch
        int scratchID_(auto const& scratch) const
        {
            if (auto const id = resourcesManager_->getID(scratch.name()))
                return *id;
            throw std::runtime_error("MHDMessenger: unregistered scratch field " + scratch.name());
        }

        // Select the field-refinement operators once at construction. The composite runtime
        // kernels are built from the configured order.
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

            // each ghost refiner gets its own patch strategy so that physical-boundary
            // ghosts are filled by the registered boundary conditions during schedule fills
            registerGhostRefinePatchStrategies_(rhoPatchStrats, info->ghostDensity);
            for (std::size_t i = 0; i < info->ghostDensity.size(); ++i)
                rhoGhostsRefiners_.addTimeRefiner(info->ghostDensity[i], info->modelDensity,
                                                  rhoOld_.name(), mhdFieldRefineOp_, fieldTimeOp_,
                                                  info->ghostDensity[i],
                                                  nonOverwriteFieldFillPattern, rhoPatchStrats[i]);


            // velGhostsRefiners_.addTimeRefiners(info->ghostVelocity, info->modelVelocity,
            //                                    Vold_.name(), mhdVecFieldRefineOp_,
            //                                    vecFieldTimeOp_,
            //                                    nonOverwriteInteriorTFfillPattern);
            //
            // pressureGhostsRefiners_.addTimeRefiners(info->ghostPressure, info->modelPressure,
            //                                         Pold_.name(), mhdFieldRefineOp_,
            //                                         fieldTimeOp_, nonOverwriteFieldFillPattern);

            registerGhostRefinePatchStrategies_(momentumPatchStrats, info->ghostMomentum);
            for (std::size_t i = 0; i < info->ghostMomentum.size(); ++i)
                momentumGhostsRefiners_.addTimeRefiner(
                    info->ghostMomentum[i], info->modelMomentum, rhoVold_.name(),
                    mhdVecFieldRefineOp_, vecFieldTimeOp_, info->ghostMomentum[i],
                    nonOverwriteInteriorTFfillPattern, momentumPatchStrats[i]);

            registerGhostRefinePatchStrategies_(totalEnergyPatchStrats, info->ghostTotalEnergy);
            for (std::size_t i = 0; i < info->ghostTotalEnergy.size(); ++i)
                totalEnergyGhostsRefiners_.addTimeRefiner(
                    info->ghostTotalEnergy[i], info->modelTotalEnergy, EtotOld_.name(),
                    mhdFieldRefineOp_, fieldTimeOp_, info->ghostTotalEnergy[i],
                    nonOverwriteFieldFillPattern, totalEnergyPatchStrats[i]);

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
            for (auto& patchStrat : magPatchStrats)
                patchStrat->registerScratchID(
                    scratchID_(model_.getTmp(core::MHDQuantity::Vector::VecAllPrimal)));

            for (std::size_t i = 0; i < info->ghostMagnetic.size(); ++i)
            {
                magGhostsRefiners_.addTimeRefiner(
                    info->ghostMagnetic[i], info->modelMagnetic, Bold_.name(), BfieldRegridOp_,
                    vecFieldTimeOp_, info->ghostMagnetic[i], nonOverwriteInteriorTFfillPattern,
                    magPatchStrats[i]);
            }

            // any B is copied into this temporary for its border max, see borderMaxB_
            auto const tmpB = model_.getTmp(core::MHDQuantity::Vector::B).name();
            magMaxRefiners_.addStaticRefiner(
                tmpB, tmpB, nullptr, tmpB,
                std::make_shared<
                    TensorFieldGhostInterpOverlapFillPattern<GridLayoutT, /*rank_=*/1>>());
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
            // Give the init refiners the model-state moment patch strategies (index 0 of each
            // ghost-strategy list: ghostX[0] == modelX == initX). The InitField schedule already
            // fills interior + coarse-fine from the coarser level via createSchedule(level,
            // nullptr, coarser, hierarchy, patchStrat) — the same call B uses in BalgoInit — and
            // with a patch strategy it now also fills the physical-boundary ghosts. So a freshly
            // created / regridded refined level touching a physical boundary carries valid moment
            // ghosts before the first flux. Default (non-overwrite) fill pattern is kept: overwrite
            // is a B/face-centered concern and corrupts the cell-centered moment interior fill.
            std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> rhoInitStrat
                = rhoPatchStrats.empty() ? nullptr : rhoPatchStrats[0];
            std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> momentumInitStrat
                = momentumPatchStrats.empty() ? nullptr : momentumPatchStrats[0];
            std::shared_ptr<SAMRAI::xfer::RefinePatchStrategy> totalEnergyInitStrat
                = totalEnergyPatchStrats.empty() ? nullptr : totalEnergyPatchStrats[0];

            densityInitRefiners_.addStaticRefiners(info->initDensity, mhdFieldRefineOp_,
                                                   info->initDensity, nullptr, rhoInitStrat);

            momentumInitRefiners_.addStaticRefiners(info->initMomentum, mhdVecFieldRefineOp_,
                                                    info->initMomentum, nullptr, momentumInitStrat);

            totalEnergyInitRefiners_.addStaticRefiners(info->initTotalEnergy, mhdFieldRefineOp_,
                                                       info->initTotalEnergy, nullptr,
                                                       totalEnergyInitStrat);
        }


        void magneticRegriding_(hierarchy_ptr const& hierarchy, level_ptr const& level,
                                level_ptr const& oldLevel, double const initDataTime)
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
                for (auto& component : vf)
                    setNaNsOnFieldGhosts(component, *patch);
        }


        FieldT rhoOld_{stratName + "rhoOld", core::MHDQuantity::Scalar::rho};
        VecFieldT Vold_{stratName + "Vold", core::MHDQuantity::Vector::V};
        FieldT Pold_{stratName + "Pold", core::MHDQuantity::Scalar::P};

        VecFieldT rhoVold_{stratName + "rhoVold", core::MHDQuantity::Vector::rhoV};
        FieldT EtotOld_{stratName + "EtotOld", core::MHDQuantity::Scalar::Etot};

        VecFieldT Bold_{stratName + "Bold", core::MHDQuantity::Vector::B};


        using rm_t = typename MHDModel::resources_manager_type;
        MHDModel& model_; // owns the temporaries used here
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        std::shared_ptr<BoundaryManagerT> boundaryManager_;
        int const firstLevel_;

        using InitRefinerPool             = RefinerPool<rm_t, RefinerType::InitField>;
        using GhostRefinerPool            = RefinerPool<rm_t, RefinerType::GhostField>;
        using PatchGhostRefinerPool       = RefinerPool<rm_t, RefinerType::PatchGhostField>;
        using InitDomPartRefinerPool      = RefinerPool<rm_t, RefinerType::InitInteriorPart>;
        using VecFieldGhostMaxRefinerPool = RefinerPool<rm_t, RefinerType::PatchVecFieldBorderMax>;


        SAMRAI::xfer::RefineAlgorithm BalgoPatchGhost; //
        SAMRAI::xfer::RefineAlgorithm BalgoInit;
        SAMRAI::xfer::RefineAlgorithm BregridAlgo;
        SAMRAI::xfer::RefineAlgorithm EalgoPatchGhost; //
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magInitRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules; //
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            magPatchGhostsRefineSchedules; //
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> elecPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            magSharedNodeRefineSchedules; //

        SAMRAI::xfer::CoarsenAlgorithm ErefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::CoarsenAlgorithm HydroXrefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::CoarsenAlgorithm HydroYrefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::CoarsenAlgorithm HydroZrefluxAlgo{SAMRAI::tbox::Dimension{dimension}};

        SAMRAI::xfer::RefineAlgorithm EpatchGhostRefluxedAlgo;
        SAMRAI::xfer::RefineAlgorithm HydroXpatchGhostRefluxedAlgo;
        SAMRAI::xfer::RefineAlgorithm HydroYpatchGhostRefluxedAlgo;
        SAMRAI::xfer::RefineAlgorithm HydroZpatchGhostRefluxedAlgo;

        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> ErefluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> HydroXrefluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> HydroYrefluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> HydroZrefluxSchedules;

        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> EpatchGhostRefluxedSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            HydroXpatchGhostRefluxedSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            HydroYpatchGhostRefluxedSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>>
            HydroZpatchGhostRefluxedSchedules;

        PatchGhostRefinerPool elecGhostsRefiners_{resourcesManager_};
        GhostRefinerPool rhoGhostsRefiners_{resourcesManager_};
        // GhostRefinerPool velGhostsRefiners_{resourcesManager_};
        // GhostRefinerPool pressureGhostsRefiners_{resourcesManager_};
        GhostRefinerPool momentumGhostsRefiners_{resourcesManager_};
        GhostRefinerPool totalEnergyGhostsRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesXGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesYGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesZGhostRefiners_{resourcesManager_};

        GhostRefinerPool magGhostsRefiners_{resourcesManager_};
        VecFieldGhostMaxRefinerPool magMaxRefiners_{resourcesManager_};

        // border max of B through the untiled temporary: B into it per patch, the schedule on it,
        //  back into B per tile (per patch if not tiled)
        void borderMaxB_(VecFieldT& B, auto& level, double const fillTime)
        {
            auto constexpr N = core::detail::tensor_field_dim_from_rank<1>();
            auto& tmp        = model_.getTmp(core::MHDQuantity::Vector::B);

            auto accessor = make_level_accessor<GridLayoutT>(level, *resourcesManager_, B, tmp);
            patch_level_exec(accessor, [](auto& view) {
                auto& [B_p, tmp_p] = view.args;
                for (std::uint8_t c = 0; c < N; ++c)
                    core::copy_fields(tmp_p[c], B_p[c]);
            });

            magMaxRefiners_.fill(tmp, level.getLevelNumber(), fillTime);

            level_exec(
                accessor, [](auto& view) -> auto& { return std::get<0>(view.args); },
                [](auto& view, auto const&... tile_idx) {
                    auto& [B_p, tmp_p] = view.args;
                    for (std::uint8_t c = 0; c < N; ++c)
                        core::copy_fields(B_p[c], tmp_p[c], tile_idx...);
                });
        }

        InitRefinerPool densityInitRefiners_{resourcesManager_};
        InitRefinerPool momentumInitRefiners_{resourcesManager_};
        InitRefinerPool totalEnergyInitRefiners_{resourcesManager_};

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

        MagneticRefinePatchStrategyT magneticRefinePatchStrategy_{*boundaryManager_};

        FieldRefinePatchStrategyList rhoPatchStrats;
        FieldRefinePatchStrategyList totalEnergyPatchStrats;
        VectorFieldRefinePatchStrategyList momentumPatchStrats;
        VectorFieldRefinePatchStrategyList elecPatchStrats;
        MagneticRefinePatchStrategyList magPatchStrats;
    };

} // namespace amr
} // namespace PHARE
#endif
