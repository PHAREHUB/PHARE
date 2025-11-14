#ifndef PHARE_MHD_MESSENGER_HPP
#define PHARE_MHD_MESSENGER_HPP

#include "amr/data/field/coarsening/electric_field_coarsener.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/mhd_flux_coarsener.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/electric_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_regrider.hpp"
#include "amr/data/field/refine/mhd_field_refiner.hpp"
#include "amr/data/field/refine/mhd_flux_refiner.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/messengers/refiner.hpp"
#include "amr/messengers/refiner_pool.hpp"
#include "amr/messengers/synchronizer_pool.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"

#include "core/data/vecfield/vecfield.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/def/phare_mpi.hpp"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"

#include <memory>
#include <string>


#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/RefineOperator.h>
#include <SAMRAI/hier/CoarsenOperator.h>


#include <memory>
#include <string>

namespace PHARE
{
namespace amr
{
    template<typename MHDModel>
    class MHDMessenger : public IMessenger<typename MHDModel::Interface>
    {
        using amr_types   = PHARE::amr::SAMRAI_Types;
        using level_t     = amr_types::level_t;
        using patch_t     = amr_types::patch_t;
        using hierarchy_t = amr_types::hierarchy_t;

        using IPhysicalModel    = MHDModel::Interface;
        using FieldT            = MHDModel::field_type;
        using VecFieldT         = MHDModel::vecfield_type;
        using MHDStateT         = MHDModel::state_type;
        using GridLayoutT       = MHDModel::gridlayout_type;
        using GridT             = MHDModel::grid_type;
        using ResourcesManagerT = MHDModel::resources_manager_type;
        using VectorFieldDataT  = TensorFieldData<1, GridLayoutT, GridT, core::MHDQuantity>;

        static constexpr auto dimension = MHDModel::dimension;

    public:
        static constexpr std::size_t rootLevelNumber = 0;
        static inline std::string const stratName    = "MHDModel-MHDModel";

        MHDMessenger(std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager,
                     int const firstLevel)
            : resourcesManager_{std::move(resourcesManager)}
            , firstLevel_{firstLevel}
        {
            // moment ghosts are primitive quantities
            resourcesManager_->registerResources(rhoOld_);
            resourcesManager_->registerResources(Vold_);
            resourcesManager_->registerResources(Pold_);

            resourcesManager_->registerResources(rhoVold_);
            resourcesManager_->registerResources(EtotOld_);

            resourcesManager_->registerResources(Jold_); // conditionally register

            // also magnetic fluxes ? or should we use static refiners instead ?
        }

        virtual ~MHDMessenger() = default;

        void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            resourcesManager_->allocate(rhoOld_, patch, allocateTime);
            resourcesManager_->allocate(Vold_, patch, allocateTime);
            resourcesManager_->allocate(Pold_, patch, allocateTime);

            resourcesManager_->allocate(rhoVold_, patch, allocateTime);
            resourcesManager_->allocate(EtotOld_, patch, allocateTime);

            resourcesManager_->allocate(Jold_, patch, allocateTime);
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

            BalgoPatchGhost.registerRefine(*b_id, *b_id, *b_id, BfieldRefineOp_,
                                           nonOverwriteInteriorTFfillPattern);

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

            EalgoPatchGhost.registerRefine(*e_id, *e_id, *e_id, EfieldRefineOp_,
                                           nonOverwriteInteriorTFfillPattern);

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



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magPatchGhostsRefineSchedules[levelNumber]
                = BalgoPatchGhost.createSchedule(level, &magneticRefinePatchStrategy_);

            elecPatchGhostsRefineSchedules[levelNumber] = EalgoPatchGhost.createSchedule(level);

            EpatchGhostRefluxedSchedules[levelNumber]
                = EpatchGhostRefluxedAlgo.createSchedule(level);
            HydroXpatchGhostRefluxedSchedules[levelNumber]
                = HydroXpatchGhostRefluxedAlgo.createSchedule(level);
            HydroYpatchGhostRefluxedSchedules[levelNumber]
                = HydroYpatchGhostRefluxedAlgo.createSchedule(level);
            HydroZpatchGhostRefluxedSchedules[levelNumber]
                = HydroZpatchGhostRefluxedAlgo.createSchedule(level);

            elecGhostsRefiners_.registerLevel(hierarchy, level);
            currentGhostsRefiners_.registerLevel(hierarchy, level);

            rhoGhostsRefiners_.registerLevel(hierarchy, level);
            // velGhostsRefiners_.registerLevel(hierarchy, level);
            // pressureGhostsRefiners_.registerLevel(hierarchy, level);

            momentumGhostsRefiners_.registerLevel(hierarchy, level);
            totalEnergyGhostsRefiners_.registerLevel(hierarchy, level);

            magFluxesXGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesYGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesZGhostRefiners_.registerLevel(hierarchy, level);

            magGhostsRefiners_.registerLevel(hierarchy, level);

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


        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                    int const levelNumber,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                    IPhysicalModel& model, double const initDataTime) override
        {
            auto& mhdModel = static_cast<MHDModel&>(model);
            auto level     = hierarchy->getPatchLevel(levelNumber);

            bool isRegriddingL0 = levelNumber == 0 and oldLevel;

            magneticRegriding_(hierarchy, level, oldLevel, initDataTime);
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

        void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                       double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();

            auto& mhdModel = static_cast<MHDModel&>(model);

            magInitRefineSchedules[levelNumber]->fillData(initDataTime);
            densityInitRefiners_.fill(levelNumber, initDataTime);
            momentumInitRefiners_.fill(levelNumber, initDataTime);
            totalEnergyInitRefiners_.fill(levelNumber, initDataTime);
        }

        void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                       std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                       double const currentTime, double const prevCoarserTIme,
                       double const newCoarserTime) final
        {
        }


        void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) final {}


        void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                         double currentTime) final
        {
            auto& mhdModel = static_cast<MHDModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, mhdModel.state.rho, mhdModel.state.V, mhdModel.state.P,
                    mhdModel.state.rhoV, mhdModel.state.Etot, mhdModel.state.J, rhoOld_, Vold_,
                    Pold_, rhoVold_, EtotOld_, Jold_);

                resourcesManager_->setTime(rhoOld_, *patch, currentTime);
                resourcesManager_->setTime(Vold_, *patch, currentTime);
                resourcesManager_->setTime(Pold_, *patch, currentTime);
                resourcesManager_->setTime(rhoVold_, *patch, currentTime);
                resourcesManager_->setTime(EtotOld_, *patch, currentTime);
                resourcesManager_->setTime(Jold_, *patch, currentTime);

                rhoOld_.copyData(mhdModel.state.rho);
                Vold_.copyData(mhdModel.state.V);
                Pold_.copyData(mhdModel.state.P);
                rhoVold_.copyData(mhdModel.state.rhoV);
                EtotOld_.copyData(mhdModel.state.Etot);
                Jold_.copyData(mhdModel.state.J);
            }
        }

        void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                            double const initDataTime) final
        {
        }

        void synchronize(SAMRAI::hier::PatchLevel& level) final {}

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

        void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                             double const time) override
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
            PHARE_LOG_SCOPE(3, "HybridHybridMessengerStrategy::fillMagneticGhosts");

            setNaNsOnVecfieldGhosts(B, level);
            magGhostsRefiners_.fill(B, level.getLevelNumber(), fillTime);
        }

        void fillCurrentGhosts(VecFieldT& J, level_t const& level, double const fillTime)
        {
            setNaNsOnVecfieldGhosts(J, level);
            currentGhostsRefiners_.fill(J, level.getLevelNumber(), fillTime);
        }

        std::string name() override { return stratName; }



    private:
        // Maybe we also need conservative ghost refiners for amr operations, actually quite
        // likely
        void registerGhostComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            // static refinement for J and E because in MHD they are temporaries, so keeping there
            // state updated after each regrid is not a priority. However if we do not correctly
            // refine on regrid, the post regrid state is not up to date (in our case it will be nan
            // since we nan-initialise) and thus is is better to rely on static refinement, which
            // uses the state after computation of ampere or CT.
            elecGhostsRefiners_.addStaticRefiners(info->ghostElectric, EfieldRefineOp_,
                                                  info->ghostElectric,
                                                  nonOverwriteInteriorTFfillPattern);

            currentGhostsRefiners_.addStaticRefiners(info->ghostCurrent, EfieldRefineOp_,
                                                     info->ghostCurrent,
                                                     nonOverwriteInteriorTFfillPattern);


            rhoGhostsRefiners_.addTimeRefiners(info->ghostDensity, info->modelDensity,
                                               rhoOld_.name(), mhdFieldRefineOp_, fieldTimeOp_,
                                               nonOverwriteFieldFillPattern);


            // velGhostsRefiners_.addTimeRefiners(info->ghostVelocity, info->modelVelocity,
            //                                    Vold_.name(), mhdVecFieldRefineOp_,
            //                                    vecFieldTimeOp_,
            //                                    nonOverwriteInteriorTFfillPattern);
            //
            // pressureGhostsRefiners_.addTimeRefiners(info->ghostPressure, info->modelPressure,
            //                                         Pold_.name(), mhdFieldRefineOp_,
            //                                         fieldTimeOp_, nonOverwriteFieldFillPattern);

            momentumGhostsRefiners_.addTimeRefiners(
                info->ghostMomentum, info->modelMomentum, rhoVold_.name(), mhdVecFieldRefineOp_,
                vecFieldTimeOp_, nonOverwriteInteriorTFfillPattern);

            totalEnergyGhostsRefiners_.addTimeRefiners(
                info->ghostTotalEnergy, info->modelTotalEnergy, EtotOld_.name(), mhdFieldRefineOp_,
                fieldTimeOp_, nonOverwriteFieldFillPattern);

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
            magneticPatchStratPerGhostRefiner_ = [&]() {
                std::vector<std::shared_ptr<
                    MagneticRefinePatchStrategy<ResourcesManagerT, VectorFieldDataT>>>
                    result;

                result.reserve(info->ghostMagnetic.size());

                for (auto const& key : info->ghostMagnetic)
                {
                    auto&& [id] = resourcesManager_->getIDsList(key);

                    auto patch_strat = std::make_shared<
                        MagneticRefinePatchStrategy<ResourcesManagerT, VectorFieldDataT>>(
                        *resourcesManager_);

                    patch_strat->registerIDs(id);

                    result.push_back(patch_strat);
                }
                return result;
            }();

            for (size_t i = 0; i < info->ghostMagnetic.size(); ++i)
            {
                magGhostsRefiners_.addStaticRefiner(
                    info->ghostMagnetic[i], BfieldRegridOp_, info->ghostMagnetic[i],
                    nonOverwriteInteriorTFfillPattern, magneticPatchStratPerGhostRefiner_[i]);
            }
        }




        // should this use conservative quantities ? When should we do the initial conversion ?
        // Maybe mhd_init
        void registerInitComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            densityInitRefiners_.addStaticRefiners(info->initDensity, mhdFieldRefineOp_,
                                                   info->initDensity);

            momentumInitRefiners_.addStaticRefiners(info->initMomentum, mhdVecFieldRefineOp_,
                                                    info->initMomentum);

            totalEnergyInitRefiners_.addStaticRefiners(info->initTotalEnergy, mhdFieldRefineOp_,
                                                       info->initTotalEnergy);
        }


        void magneticRegriding_(std::shared_ptr<hierarchy_t> const& hierarchy,
                                std::shared_ptr<level_t> const& level,
                                std::shared_ptr<level_t> const& oldLevel, double const initDataTime)
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
            // not gbox is a fieldBox (thanks to the layout)

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


        FieldT rhoOld_{stratName + "rhoOld", core::MHDQuantity::Scalar::rho};
        VecFieldT Vold_{stratName + "Vold", core::MHDQuantity::Vector::V};
        FieldT Pold_{stratName + "Pold", core::MHDQuantity::Scalar::P};

        VecFieldT rhoVold_{stratName + "rhoVold", core::MHDQuantity::Vector::rhoV};
        FieldT EtotOld_{stratName + "EtotOld", core::MHDQuantity::Scalar::Etot};

        VecFieldT Jold_{stratName + "Jold", core::MHDQuantity::Vector::J};


        using rm_t = typename MHDModel::resources_manager_type;
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        int const firstLevel_;

        using InitRefinerPool        = RefinerPool<rm_t, RefinerType::InitField>;
        using GhostRefinerPool       = RefinerPool<rm_t, RefinerType::GhostField>;
        using InitDomPartRefinerPool = RefinerPool<rm_t, RefinerType::InitInteriorPart>;


        SAMRAI::xfer::RefineAlgorithm BalgoPatchGhost;
        SAMRAI::xfer::RefineAlgorithm BalgoInit;
        SAMRAI::xfer::RefineAlgorithm BregridAlgo;
        SAMRAI::xfer::RefineAlgorithm EalgoPatchGhost;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magInitRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> elecPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magSharedNodeRefineSchedules;

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

        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};
        GhostRefinerPool currentGhostsRefiners_{resourcesManager_};
        GhostRefinerPool rhoGhostsRefiners_{resourcesManager_};
        // GhostRefinerPool velGhostsRefiners_{resourcesManager_};
        // GhostRefinerPool pressureGhostsRefiners_{resourcesManager_};
        GhostRefinerPool momentumGhostsRefiners_{resourcesManager_};
        GhostRefinerPool totalEnergyGhostsRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesXGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesYGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesZGhostRefiners_{resourcesManager_};

        GhostRefinerPool magGhostsRefiners_{resourcesManager_};

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

        template<typename Policy>
        using FieldRefineOp = FieldRefineOperator<GridLayoutT, GridT, Policy>;

        template<typename Policy>
        using VecFieldRefineOp = VecFieldRefineOperator<GridLayoutT, GridT, Policy>;

        using DefaultVecFieldRefineOp = VecFieldRefineOp<DefaultFieldRefiner<dimension>>;
        using MagneticFieldRefineOp   = VecFieldRefineOp<MagneticFieldRefiner<dimension>>;
        using MagneticFieldRegridOp   = VecFieldRefineOp<MagneticFieldRegrider<dimension>>;
        using ElectricFieldRefineOp   = VecFieldRefineOp<ElectricFieldRefiner<dimension>>;

        using MHDFluxRefineOp     = FieldRefineOp<MHDFluxRefiner<dimension>>;
        using MHDVecFluxRefineOp  = VecFieldRefineOp<MHDFluxRefiner<dimension>>;
        using MHDFieldRefineOp    = FieldRefineOp<MHDFieldRefiner<dimension>>;
        using MHDVecFieldRefineOp = VecFieldRefineOp<MHDFieldRefiner<dimension>>;

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

        RefOp_ptr mhdFluxRefineOp_{std::make_shared<MHDFluxRefineOp>()};
        RefOp_ptr mhdVecFluxRefineOp_{std::make_shared<MHDVecFluxRefineOp>()};
        RefOp_ptr mhdFieldRefineOp_{std::make_shared<MHDFieldRefineOp>()};
        RefOp_ptr mhdVecFieldRefineOp_{std::make_shared<MHDVecFieldRefineOp>()};
        RefOp_ptr EfieldRefineOp_{std::make_shared<ElectricFieldRefineOp>()};
        RefOp_ptr BfieldRefineOp_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr BfieldRegridOp_{std::make_shared<MagneticFieldRegridOp>()};

        TimeOp_ptr fieldTimeOp_{std::make_shared<FieldTimeInterp>()};
        TimeOp_ptr vecFieldTimeOp_{std::make_shared<VecFieldTimeInterp>()};

        using TensorFieldFillPattern_t = TensorFieldFillPattern<dimension /*, rank=1*/>;
        using FieldFillPattern_t       = FieldFillPattern<dimension>;

        std::shared_ptr<FieldFillPattern_t> nonOverwriteFieldFillPattern
            = std::make_shared<FieldFillPattern<dimension>>(); // stateless (mostly)

        std::shared_ptr<TensorFieldFillPattern_t> nonOverwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>();

        std::shared_ptr<TensorFieldFillPattern_t> overwriteInteriorTFfillPattern
            = std::make_shared<TensorFieldFillPattern<dimension /*, rank=1*/>>(
                /*overwrite_interior=*/true);

        CoarsenOp_ptr mhdFluxCoarseningOp_{std::make_shared<MHDFluxCoarsenOp>()};
        CoarsenOp_ptr mhdVecFluxCoarseningOp_{std::make_shared<MHDVecFluxCoarsenOp>()};
        CoarsenOp_ptr electricFieldCoarseningOp_{std::make_shared<ElectricFieldCoarsenOp>()};

        MagneticRefinePatchStrategy<ResourcesManagerT, VectorFieldDataT>
            magneticRefinePatchStrategy_{*resourcesManager_};

        std::vector<
            std::shared_ptr<MagneticRefinePatchStrategy<ResourcesManagerT, VectorFieldDataT>>>
            magneticPatchStratPerGhostRefiner_;
    };

} // namespace amr
} // namespace PHARE
#endif
