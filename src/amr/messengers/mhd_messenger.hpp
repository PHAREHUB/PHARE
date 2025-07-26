#ifndef PHARE_MHD_MESSENGER_HPP
#define PHARE_MHD_MESSENGER_HPP

#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/field/coarsening/electric_field_coarsener.hpp"
#include "amr/data/field/coarsening/magnetic_field_coarsener.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/mhd_flux_coarsener.hpp"
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
#include "amr/messengers/component_variable_field_paterns.hpp"
#include "amr/data/field/refine/magnetic_refine_patch_strategy.hpp"

#include "core/data/vecfield/vecfield.hpp"
#include "core/mhd/mhd_quantities.hpp"
#include "core/def/phare_mpi.hpp"

#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"

#include <memory>
#include <string>

namespace PHARE
{
namespace amr
{
    template<typename MHDModel>
    class MHDMessenger : public IMessenger<typename MHDModel::Interface>
    {
        using IPhysicalModel    = typename MHDModel::Interface;
        using FieldT            = typename MHDModel::field_type;
        using VecFieldT         = typename MHDModel::vecfield_type;
        using MHDStateT         = typename MHDModel::state_type;
        using GridLayoutT       = typename MHDModel::gridlayout_type;
        using GridT             = typename MHDModel::grid_type;
        using ResourcesManagerT = typename MHDModel::resources_manager_type;
        using FieldDataT        = FieldData<GridLayoutT, GridT>;

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

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> xVariableFillPattern
                = std::make_shared<XVariableFillPattern>();

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> yVariableFillPattern
                = std::make_shared<YVariableFillPattern>();

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> zVariableFillPattern
                = std::make_shared<ZVariableFillPattern>();


            auto bx_id = resourcesManager_->getID(mhdInfo->modelMagnetic.xName);
            auto by_id = resourcesManager_->getID(mhdInfo->modelMagnetic.yName);
            auto bz_id = resourcesManager_->getID(mhdInfo->modelMagnetic.zName);

            if (!bx_id or !by_id or !bz_id)
            {
                throw std::runtime_error("MHDMessenger: missing magnetic field variable IDs");
            }

            Balgo.registerRefine(*bx_id, *bx_id, *bx_id, BfieldRefineOp_x_, xVariableFillPattern);
            Balgo.registerRefine(*by_id, *by_id, *by_id, BfieldRefineOp_y_, yVariableFillPattern);
            Balgo.registerRefine(*bz_id, *bz_id, *bz_id, BfieldRefineOp_z_, zVariableFillPattern);

            BalgoNode.registerRefine(*bx_id, *bx_id, *bx_id, BfieldNodeRefineOp_x_,
                                     xVariableFillPattern);
            BalgoNode.registerRefine(*by_id, *by_id, *by_id, BfieldNodeRefineOp_y_,
                                     yVariableFillPattern);
            BalgoNode.registerRefine(*bz_id, *bz_id, *bz_id, BfieldNodeRefineOp_z_,
                                     zVariableFillPattern);

            BcopyAlgo.registerRefine(*bx_id, *bx_id, *bx_id, BfieldRefineOp_x_,
                                     xVariableFillPattern);
            BcopyAlgo.registerRefine(*by_id, *by_id, *by_id, BfieldRefineOp_y_,
                                     yVariableFillPattern);
            BcopyAlgo.registerRefine(*bz_id, *bz_id, *bz_id, BfieldRefineOp_z_,
                                     zVariableFillPattern);

            BregridAlgo.registerRefine(*bx_id, *bx_id, *bx_id, BfieldRegridOp_x_,
                                       xVariableFillPattern);
            BregridAlgo.registerRefine(*by_id, *by_id, *by_id, BfieldRegridOp_y_,
                                       yVariableFillPattern);
            BregridAlgo.registerRefine(*bz_id, *bz_id, *bz_id, BfieldRegridOp_z_,
                                       zVariableFillPattern);

            magneticRefinePatchStrategy_.registerIDs(*bx_id, *by_id, *bz_id);

            auto ex_id = resourcesManager_->getID(mhdInfo->modelElectric.xName);
            auto ey_id = resourcesManager_->getID(mhdInfo->modelElectric.yName);
            auto ez_id = resourcesManager_->getID(mhdInfo->modelElectric.zName);

            if (!ex_id or !ey_id or !ez_id)
            {
                throw std::runtime_error("MHDMessenger: missing electric field variable IDs");
            }

            Ealgo.registerRefine(*ex_id, *ex_id, *ex_id, EfieldRefineOp_x_, xVariableFillPattern);
            Ealgo.registerRefine(*ey_id, *ey_id, *ey_id, EfieldRefineOp_y_, yVariableFillPattern);
            Ealgo.registerRefine(*ez_id, *ez_id, *ez_id, EfieldRefineOp_z_, zVariableFillPattern);

            // refluxing
            // we first want to coarsen the flux sum onto the coarser level
            auto rho_fx_reflux_id   = resourcesManager_->getID(mhdInfo->reflux.rho_fx);
            auto rhoVx_fx_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fx.xName);
            auto rhoVy_fx_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fx.yName);
            auto rhoVz_fx_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fx.zName);
            auto Etot_fx_reflux_id  = resourcesManager_->getID(mhdInfo->reflux.Etot_fx);

            if (!rho_fx_reflux_id or !rhoVx_fx_reflux_id or !rhoVy_fx_reflux_id
                or !rhoVz_fx_reflux_id or !Etot_fx_reflux_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing reflux variable IDs for fluxes in x direction");
            }

            auto rho_fx_fluxsum_id   = resourcesManager_->getID(mhdInfo->fluxSum.rho_fx);
            auto rhoVx_fx_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fx.xName);
            auto rhoVy_fx_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fx.yName);
            auto rhoVz_fx_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fx.zName);
            auto Etot_fx_fluxsum_id  = resourcesManager_->getID(mhdInfo->fluxSum.Etot_fx);


            if (!rho_fx_fluxsum_id or !rhoVx_fx_fluxsum_id or !rhoVy_fx_fluxsum_id
                or !rhoVz_fx_fluxsum_id or !Etot_fx_fluxsum_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing flux sum variable IDs for fluxes in x direction");
            }


            // all of the fluxes fx are defined on the same faces no matter the component, so we
            // just need a different fill pattern per direction
            RefluxAlgo.registerCoarsen(*rho_fx_reflux_id, *rho_fx_fluxsum_id,
                                       mhdFluxCoarseningOp_x_, xVariableFillPattern);
            RefluxAlgo.registerCoarsen(*rhoVx_fx_reflux_id, *rhoVx_fx_fluxsum_id,
                                       mhdFluxCoarseningOp_x_, xVariableFillPattern);
            RefluxAlgo.registerCoarsen(*rhoVy_fx_reflux_id, *rhoVy_fx_fluxsum_id,
                                       mhdFluxCoarseningOp_x_, xVariableFillPattern);
            RefluxAlgo.registerCoarsen(*rhoVz_fx_reflux_id, *rhoVz_fx_fluxsum_id,
                                       mhdFluxCoarseningOp_x_, xVariableFillPattern);
            RefluxAlgo.registerCoarsen(*Etot_fx_reflux_id, *Etot_fx_fluxsum_id,
                                       mhdFluxCoarseningOp_x_, xVariableFillPattern);

            // we then need to refill the ghosts so that they agree with the newly refluxed
            // cells
            PatchGhostRefluxedAlgo.registerRefine(*rho_fx_reflux_id, *rho_fx_reflux_id,
                                                  *rho_fx_reflux_id, mhdFluxRefineOp_x_,
                                                  xVariableFillPattern);
            PatchGhostRefluxedAlgo.registerRefine(*rhoVx_fx_reflux_id, *rhoVx_fx_reflux_id,
                                                  *rhoVx_fx_reflux_id, mhdFluxRefineOp_x_,
                                                  xVariableFillPattern);
            PatchGhostRefluxedAlgo.registerRefine(*rhoVy_fx_reflux_id, *rhoVy_fx_reflux_id,
                                                  *rhoVy_fx_reflux_id, mhdFluxRefineOp_x_,
                                                  xVariableFillPattern);
            PatchGhostRefluxedAlgo.registerRefine(*rhoVz_fx_reflux_id, *rhoVz_fx_reflux_id,
                                                  *rhoVz_fx_reflux_id, mhdFluxRefineOp_x_,
                                                  xVariableFillPattern);
            PatchGhostRefluxedAlgo.registerRefine(*Etot_fx_reflux_id, *Etot_fx_reflux_id,
                                                  *Etot_fx_reflux_id, mhdFluxRefineOp_x_,
                                                  xVariableFillPattern);

            if constexpr (dimension >= 2)
            {
                auto rho_fy_reflux_id   = resourcesManager_->getID(mhdInfo->reflux.rho_fy);
                auto rhoVx_fy_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fy.xName);
                auto rhoVy_fy_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fy.yName);
                auto rhoVz_fy_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rhoV_fy.zName);
                auto Etot_fy_reflux_id  = resourcesManager_->getID(mhdInfo->reflux.Etot_fy);

                if (!rho_fy_reflux_id or !rhoVx_fy_reflux_id or !rhoVy_fy_reflux_id
                    or !rhoVz_fy_reflux_id or !Etot_fy_reflux_id)
                {
                    throw std::runtime_error(
                        "MHDMessenger: missing reflux variable IDs for fluxes in y direction");
                }

                auto rho_fy_fluxsum_id   = resourcesManager_->getID(mhdInfo->fluxSum.rho_fy);
                auto rhoVx_fy_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fy.xName);
                auto rhoVy_fy_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fy.yName);
                auto rhoVz_fy_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fy.zName);
                auto Etot_fy_fluxsum_id  = resourcesManager_->getID(mhdInfo->fluxSum.Etot_fy);

                if (!rho_fy_fluxsum_id or !rhoVx_fy_fluxsum_id or !rhoVy_fy_fluxsum_id
                    or !rhoVz_fy_fluxsum_id or !Etot_fy_fluxsum_id)
                {
                    throw std::runtime_error(
                        "MHDMessenger: missing flux sum variable IDs for fluxes in y direction");
                }

                RefluxAlgo.registerCoarsen(*rho_fy_reflux_id, *rho_fy_fluxsum_id,
                                           mhdFluxCoarseningOp_y_, yVariableFillPattern);
                RefluxAlgo.registerCoarsen(*rhoVx_fy_reflux_id, *rhoVx_fy_fluxsum_id,
                                           mhdFluxCoarseningOp_y_, yVariableFillPattern);
                RefluxAlgo.registerCoarsen(*rhoVy_fy_reflux_id, *rhoVy_fy_fluxsum_id,
                                           mhdFluxCoarseningOp_y_, yVariableFillPattern);
                RefluxAlgo.registerCoarsen(*rhoVz_fy_reflux_id, *rhoVz_fy_fluxsum_id,
                                           mhdFluxCoarseningOp_y_, yVariableFillPattern);
                RefluxAlgo.registerCoarsen(*Etot_fy_reflux_id, *Etot_fy_fluxsum_id,
                                           mhdFluxCoarseningOp_y_, yVariableFillPattern);

                PatchGhostRefluxedAlgo.registerRefine(*rho_fy_reflux_id, *rho_fy_reflux_id,
                                                      *rho_fy_reflux_id, mhdFluxRefineOp_y_,
                                                      yVariableFillPattern);
                PatchGhostRefluxedAlgo.registerRefine(*rhoVx_fy_reflux_id, *rhoVx_fy_reflux_id,
                                                      *rhoVx_fy_reflux_id, mhdFluxRefineOp_y_,
                                                      yVariableFillPattern);
                PatchGhostRefluxedAlgo.registerRefine(*rhoVy_fy_reflux_id, *rhoVy_fy_reflux_id,
                                                      *rhoVy_fy_reflux_id, mhdFluxRefineOp_y_,
                                                      yVariableFillPattern);
                PatchGhostRefluxedAlgo.registerRefine(*rhoVz_fy_reflux_id, *rhoVz_fy_reflux_id,
                                                      *rhoVz_fy_reflux_id, mhdFluxRefineOp_y_,
                                                      yVariableFillPattern);
                PatchGhostRefluxedAlgo.registerRefine(*Etot_fy_reflux_id, *Etot_fy_reflux_id,
                                                      *Etot_fy_reflux_id, mhdFluxRefineOp_y_,
                                                      yVariableFillPattern);

                if constexpr (dimension == 3)
                {
                    auto rho_fz_reflux_id = resourcesManager_->getID(mhdInfo->reflux.rho_fz);
                    auto rhoVx_fz_reflux_id
                        = resourcesManager_->getID(mhdInfo->reflux.rhoV_fz.xName);
                    auto rhoVy_fz_reflux_id
                        = resourcesManager_->getID(mhdInfo->reflux.rhoV_fz.yName);
                    auto rhoVz_fz_reflux_id
                        = resourcesManager_->getID(mhdInfo->reflux.rhoV_fz.zName);
                    auto Etot_fz_reflux_id = resourcesManager_->getID(mhdInfo->reflux.Etot_fz);


                    if (!rho_fz_reflux_id or !rhoVx_fz_reflux_id or !rhoVy_fz_reflux_id
                        or !rhoVz_fz_reflux_id or !Etot_fz_reflux_id)
                    {
                        throw std::runtime_error(
                            "MHDMessenger: missing reflux variable IDs for fluxes in z direction");
                    }

                    auto rho_fz_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.rho_fz);
                    auto rhoVx_fz_fluxsum_id
                        = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fz.xName);
                    auto rhoVy_fz_fluxsum_id
                        = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fz.yName);
                    auto rhoVz_fz_fluxsum_id
                        = resourcesManager_->getID(mhdInfo->fluxSum.rhoV_fz.zName);
                    auto Etot_fz_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSum.Etot_fz);

                    if (!rho_fz_fluxsum_id or !rhoVx_fz_fluxsum_id or !rhoVy_fz_fluxsum_id
                        or !rhoVz_fz_fluxsum_id or !Etot_fz_fluxsum_id)
                    {
                        throw std::runtime_error("MHDMessenger: missing flux sum variable IDs for "
                                                 "fluxes in z direction");
                    }

                    RefluxAlgo.registerCoarsen(*rho_fz_reflux_id, *rho_fz_fluxsum_id,
                                               mhdFluxCoarseningOp_z_, zVariableFillPattern);
                    RefluxAlgo.registerCoarsen(*rhoVx_fz_reflux_id, *rhoVx_fz_fluxsum_id,
                                               mhdFluxCoarseningOp_z_, zVariableFillPattern);
                    RefluxAlgo.registerCoarsen(*rhoVy_fz_reflux_id, *rhoVy_fz_fluxsum_id,
                                               mhdFluxCoarseningOp_z_, zVariableFillPattern);
                    RefluxAlgo.registerCoarsen(*rhoVz_fz_reflux_id, *rhoVz_fz_fluxsum_id,
                                               mhdFluxCoarseningOp_z_, zVariableFillPattern);
                    RefluxAlgo.registerCoarsen(*Etot_fz_reflux_id, *Etot_fz_fluxsum_id,
                                               mhdFluxCoarseningOp_z_, zVariableFillPattern);



                    PatchGhostRefluxedAlgo.registerRefine(*rho_fz_reflux_id, *rho_fz_reflux_id,
                                                          *rho_fz_reflux_id, mhdFluxRefineOp_z_,
                                                          zVariableFillPattern);
                    PatchGhostRefluxedAlgo.registerRefine(*rhoVx_fz_reflux_id, *rhoVx_fz_reflux_id,
                                                          *rhoVx_fz_reflux_id, mhdFluxRefineOp_z_,
                                                          zVariableFillPattern);
                    PatchGhostRefluxedAlgo.registerRefine(*rhoVy_fz_reflux_id, *rhoVy_fz_reflux_id,
                                                          *rhoVy_fz_reflux_id, mhdFluxRefineOp_z_,
                                                          zVariableFillPattern);
                    PatchGhostRefluxedAlgo.registerRefine(*rhoVz_fz_reflux_id, *rhoVz_fz_reflux_id,
                                                          *rhoVz_fz_reflux_id, mhdFluxRefineOp_z_,
                                                          zVariableFillPattern);
                    PatchGhostRefluxedAlgo.registerRefine(*Etot_fz_reflux_id, *Etot_fz_reflux_id,
                                                          *Etot_fz_reflux_id, mhdFluxRefineOp_z_,
                                                          zVariableFillPattern);
                }
            }

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> yzVariableFillPattern
                = std::make_shared<YZVariableFillPattern>();

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> xzVariableFillPattern
                = std::make_shared<XZVariableFillPattern>();

            std::shared_ptr<SAMRAI::xfer::VariableFillPattern> xyVariableFillPattern
                = std::make_shared<XYVariableFillPattern>();

            auto ex_reflux_id = resourcesManager_->getID(mhdInfo->refluxElectric.xName);
            auto ey_reflux_id = resourcesManager_->getID(mhdInfo->refluxElectric.yName);
            auto ez_reflux_id = resourcesManager_->getID(mhdInfo->refluxElectric.zName);

            auto ex_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSumElectric.xName);
            auto ey_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSumElectric.yName);
            auto ez_fluxsum_id = resourcesManager_->getID(mhdInfo->fluxSumElectric.zName);

            if (!ex_reflux_id or !ey_reflux_id or !ez_reflux_id or !ex_fluxsum_id or !ey_fluxsum_id
                or !ez_fluxsum_id)
            {
                throw std::runtime_error(
                    "MHDMessenger: missing electric refluxing field variable IDs");
            }

            RefluxAlgo.registerCoarsen(*ex_reflux_id, *ex_fluxsum_id, electricFieldCoarseningOp_x_,
                                       yzVariableFillPattern);
            RefluxAlgo.registerCoarsen(*ey_reflux_id, *ey_fluxsum_id, electricFieldCoarseningOp_y_,
                                       xzVariableFillPattern);
            RefluxAlgo.registerCoarsen(*ez_reflux_id, *ez_fluxsum_id, electricFieldCoarseningOp_z_,
                                       xyVariableFillPattern);

            PatchGhostRefluxedAlgo.registerRefine(*ex_reflux_id, *ex_reflux_id, *ex_reflux_id,
                                                  EfieldRefineOp_, yzVariableFillPattern);
            PatchGhostRefluxedAlgo.registerRefine(*ey_reflux_id, *ey_reflux_id, *ey_reflux_id,
                                                  EfieldRefineOp_, xzVariableFillPattern);
            PatchGhostRefluxedAlgo.registerRefine(*ez_reflux_id, *ez_reflux_id, *ez_reflux_id,
                                                  EfieldRefineOp_, xyVariableFillPattern);

            registerGhostComms_(mhdInfo);
            registerInitComms_(mhdInfo);
            // registerSyncComms_(mhdInfo);
        }



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magSharedNodeRefineSchedules[levelNumber]
                = BalgoNode.createSchedule(level, &magneticRefinePatchStrategy_);

            magPatchGhostsRefineSchedules[levelNumber]
                = Balgo.createSchedule(level, &magneticRefinePatchStrategy_);

            elecPatchGhostsRefineSchedules[levelNumber] = Ealgo.createSchedule(level);

            magGhostsRefineSchedules[levelNumber] = Balgo.createSchedule(
                level, levelNumber - 1, hierarchy, &magneticRefinePatchStrategy_);

            patchGhostRefluxedSchedules[levelNumber] = PatchGhostRefluxedAlgo.createSchedule(level);

            elecSharedNodesRefiners_.registerLevel(hierarchy, level);
            elecGhostsRefiners_.registerLevel(hierarchy, level);
            currentGhostsRefiners_.registerLevel(hierarchy, level);

            rhoGhostsRefiners_.registerLevel(hierarchy, level);
            velGhostsRefiners_.registerLevel(hierarchy, level);
            pressureGhostsRefiners_.registerLevel(hierarchy, level);

            momentumGhostsRefiners_.registerLevel(hierarchy, level);
            totalEnergyGhostsRefiners_.registerLevel(hierarchy, level);

            magFluxesXSharedNodesRefiners_.registerLevel(hierarchy, level);
            magFluxesYSharedNodesRefiners_.registerLevel(hierarchy, level);
            magFluxesZSharedNodesRefiners_.registerLevel(hierarchy, level);

            magFluxesXGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesYGhostRefiners_.registerLevel(hierarchy, level);
            magFluxesZGhostRefiners_.registerLevel(hierarchy, level);

            if (levelNumber != rootLevelNumber)
            {
                // refluxing
                auto const& coarseLevel      = hierarchy->getPatchLevel(levelNumber - 1);
                refluxSchedules[levelNumber] = RefluxAlgo.createSchedule(coarseLevel, level);

                // refinement
                magInitRefineSchedules[levelNumber] = Balgo.createSchedule(
                    level, nullptr, levelNumber - 1, hierarchy, &magneticRefinePatchStrategy_);

                densityInitRefiners_.registerLevel(hierarchy, level);
                momentumInitRefiners_.registerLevel(hierarchy, level);
                totalEnergyInitRefiners_.registerLevel(hierarchy, level);

                // densitySynchronizers_.registerLevel(hierarchy, level);
                // momentumSynchronizers_.registerLevel(hierarchy, level);
                // magnetoSynchronizers_.registerLevel(hierarchy, level);
                // totalEnergySynchronizers_.registerLevel(hierarchy, level);
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

            magneticRegriding_(hierarchy, level, oldLevel, mhdModel, initDataTime);
            densityInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            momentumInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            totalEnergyInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);

            magPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);
            elecPatchGhostsRefineSchedules[levelNumber]->fillData(initDataTime);
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

        void synchronize(SAMRAI::hier::PatchLevel& level) final
        {
            // auto levelNumber = level.getLevelNumber();
            //
            // densitySynchronizers_.sync(levelNumber);
            // momentumSynchronizers_.sync(levelNumber);
            // magnetoSynchronizers_.sync(levelNumber);
            // totalEnergySynchronizers_.sync(levelNumber);
        }

        void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                    double const syncTime) override
        {
            refluxSchedules[fineLevelNumber]->coarsenData();
            patchGhostRefluxedSchedules[coarserLevelNumber]->fillData(syncTime);
        }

        void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                             double const time) override
        {
            auto levelNumber = level.getLevelNumber();
            auto& mhdModel   = static_cast<MHDModel&>(model);

            std::cout << "MHDMessenger::postSynchronize: levelNumber = " << levelNumber
                      << ", time = " << time << std::endl;
            magSharedNodeRefineSchedules[levelNumber]->fillData(time);
            // magPatchGhostsRefineSchedules[levelNumber]->fillData(time);
        }

        void fillMomentsGhosts(MHDStateT& state, int const levelNumber, double const fillTime)
        {
            rhoGhostsRefiners_.fill(state.rho, levelNumber, fillTime);
            velGhostsRefiners_.fill(state.V, levelNumber, fillTime);
            pressureGhostsRefiners_.fill(state.P, levelNumber, fillTime);
        }

        void fillMagneticFluxesXGhosts(VecFieldT& Fx_B, int const levelNumber,
                                       double const fillTime)
        {
            magFluxesXSharedNodesRefiners_.fill(Fx_B, levelNumber, fillTime);
            magFluxesXGhostRefiners_.fill(Fx_B, levelNumber, fillTime);
        }

        void fillMagneticFluxesYGhosts(VecFieldT& Fy_B, int const levelNumber,
                                       double const fillTime)
        {
            magFluxesYSharedNodesRefiners_.fill(Fy_B, levelNumber, fillTime);
            magFluxesYGhostRefiners_.fill(Fy_B, levelNumber, fillTime);
        }

        void fillMagneticFluxesZGhosts(VecFieldT& Fz_B, int const levelNumber,
                                       double const fillTime)
        {
            magFluxesZSharedNodesRefiners_.fill(Fz_B, levelNumber, fillTime);
            magFluxesZGhostRefiners_.fill(Fz_B, levelNumber, fillTime);
        }

        void fillElectricGhosts(VecFieldT& E, int const levelNumber, double const fillTime)
        {
            elecSharedNodesRefiners_.fill(E, levelNumber, fillTime);
            elecGhostsRefiners_.fill(E, levelNumber, fillTime);
        }

        std::string name() override { return stratName; }



    private:
        // Maybe we also need conservative ghost refiners for amr operations, actually quite
        // likely
        void registerGhostComms_(std::unique_ptr<MHDMessengerInfo> const& info)
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

            currentGhostsRefiners_.addTimeRefiners(info->ghostCurrent, info->modelCurrent,
                                                   core::VecFieldNames{Jold_}, EfieldRefineOp_,
                                                   fieldTimeOp_);

            rhoGhostsRefiners_.addTimeRefiners(info->ghostDensity, info->modelDensity,
                                               rhoOld_.name(), mhdFieldRefineOp_, fieldTimeOp_);


            velGhostsRefiners_.addTimeRefiners(info->ghostVelocity, info->modelVelocity,
                                               core::VecFieldNames{Vold_}, mhdFieldRefineOp_,
                                               fieldTimeOp_);

            pressureGhostsRefiners_.addTimeRefiners(info->ghostPressure, info->modelPressure,
                                                    Pold_.name(), mhdFieldRefineOp_, fieldTimeOp_);

            momentumGhostsRefiners_.addTimeRefiners(info->ghostMomentum, info->modelMomentum,
                                                    core::VecFieldNames{rhoVold_},
                                                    mhdFieldRefineOp_, fieldTimeOp_);

            totalEnergyGhostsRefiners_.addTimeRefiners(info->ghostTotalEnergy,
                                                       info->modelTotalEnergy, EtotOld_.name(),
                                                       mhdFieldRefineOp_, fieldTimeOp_);


            magFluxesXSharedNodesRefiners_.addStaticRefiners(info->ghostMagneticFluxesX,
                                                             mhdFluxNodeRefineOp_,
                                                             makeKeys(info->ghostMagneticFluxesX));

            magFluxesYSharedNodesRefiners_.addStaticRefiners(info->ghostMagneticFluxesY,
                                                             mhdFluxNodeRefineOp_,
                                                             makeKeys(info->ghostMagneticFluxesY));

            magFluxesZSharedNodesRefiners_.addStaticRefiners(info->ghostMagneticFluxesZ,
                                                             mhdFluxNodeRefineOp_,
                                                             makeKeys(info->ghostMagneticFluxesZ));

            magFluxesXGhostRefiners_.addStaticRefiners(info->ghostMagneticFluxesX, mhdFluxRefineOp_,
                                                       makeKeys(info->ghostMagneticFluxesX));

            magFluxesYGhostRefiners_.addStaticRefiners(info->ghostMagneticFluxesY, mhdFluxRefineOp_,
                                                       makeKeys(info->ghostMagneticFluxesY));

            magFluxesZGhostRefiners_.addStaticRefiners(info->ghostMagneticFluxesZ, mhdFluxRefineOp_,
                                                       makeKeys(info->ghostMagneticFluxesZ));
        }




        // should this use conservative quantities ? When should we do the initial conversion ?
        // Maybe mhd_init
        void registerInitComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& descriptor) {
                std::vector<std::string> keys;
                std::transform(std::begin(descriptor), std::end(descriptor),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            densityInitRefiners_.addStaticRefiners(info->initDensity, mhdFieldRefineOp_,
                                                   info->initDensity);

            momentumInitRefiners_.addStaticRefiners(info->initMomentum, mhdFieldRefineOp_,
                                                    makeKeys(info->initMomentum));

            totalEnergyInitRefiners_.addStaticRefiners(info->initTotalEnergy, mhdFieldRefineOp_,
                                                       info->initTotalEnergy);
        }

        // void registerSyncComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        // {
        //     densitySynchronizers_.add(info->modelDensity, fieldCoarseningOp_,
        //     info->modelDensity);
        //
        //     momentumSynchronizers_.add(info->modelMomentum, fieldCoarseningOp_,
        //                                info->modelMomentum.vecName);
        //
        //     magnetoSynchronizers_.add(info->modelMagnetic, magneticCoarseningOp_,
        //                               info->modelMagnetic.vecName);
        //
        //     totalEnergySynchronizers_.add(info->modelTotalEnergy, fieldCoarseningOp_,
        //                                   info->modelTotalEnergy);
        // }


        void magneticRegriding_(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                std::shared_ptr<SAMRAI::hier::PatchLevel> const& level,
                                std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                                MHDModel& mhdModel, double const initDataTime)
        {
            for (auto& patch : *level)
            {
                auto& B          = mhdModel.state.B;
                auto layout      = amr::layoutFromPatch<GridLayoutT>(*patch);
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, B);

                layout.evalOnGhostBox(B(core::Component::X), [&](auto const&... args) mutable {
                    B(core::Component::X)(args...) = std::numeric_limits<double>::quiet_NaN();
                });

                layout.evalOnGhostBox(B(core::Component::Y), [&](auto const&... args) mutable {
                    B(core::Component::Y)(args...) = std::numeric_limits<double>::quiet_NaN();
                });

                layout.evalOnGhostBox(B(core::Component::Z), [&](auto const&... args) mutable {
                    B(core::Component::Z)(args...) = std::numeric_limits<double>::quiet_NaN();
                });
            }

            std::cout << "old level ? "
                      << (oldLevel ? std::to_string(oldLevel->getLevelNumber()) : "null") << "\n";

            auto regridCopySchedule = BcopyAlgo.createSchedule(level, oldLevel);
            regridCopySchedule->fillData(initDataTime);

            auto magSchedule = BregridAlgo.createSchedule(
                level, oldLevel, level->getNextCoarserHierarchyLevelNumber(), hierarchy,
                &magneticRefinePatchStrategy_);
            magSchedule->fillData(initDataTime);

            for (auto& patch : *level)
            {
                auto& B          = mhdModel.state.B;
                auto layout      = amr::layoutFromPatch<GridLayoutT>(*patch);
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, B);

                auto notNan = [&](auto& b, core::MeshIndex<dimension> idx) {
                    auto check = [&](auto&&... indices) {
                        if (std::isnan(b(indices...)))
                        {
                            std::string index_str;
                            ((index_str
                              += (index_str.empty() ? "" : ", ") + std::to_string(indices)),
                             ...);
                            std::cout << "NaN found in magnetic field " + b.name() + " at index ("
                                             + index_str + ")"
                                      << "\n";
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

        FieldT rhoOld_{stratName + "rhoOld", core::MHDQuantity::Scalar::rho};
        VecFieldT Vold_{stratName + "Vold", core::MHDQuantity::Vector::V};
        FieldT Pold_{stratName + "Pold", core::MHDQuantity::Scalar::P};

        VecFieldT rhoVold_{stratName + "rhoVold", core::MHDQuantity::Vector::rhoV};
        FieldT EtotOld_{stratName + "EtotOld", core::MHDQuantity::Scalar::Etot};

        VecFieldT Jold_{stratName + "Jold", core::MHDQuantity::Vector::J};

        using rm_t = typename MHDModel::resources_manager_type;
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        int const firstLevel_;

        using InitRefinerPool           = RefinerPool<rm_t, RefinerType::InitField>;
        using SharedNodeRefinerPool     = RefinerPool<rm_t, RefinerType::SharedBorder>;
        using GhostRefinerPool          = RefinerPool<rm_t, RefinerType::GhostField>;
        using PatchGhostRefinerPool     = RefinerPool<rm_t, RefinerType::PatchGhostField>;
        using InitDomPartRefinerPool    = RefinerPool<rm_t, RefinerType::InitInteriorPart>;
        using PatchGhostPartRefinerPool = RefinerPool<rm_t, RefinerType::InteriorGhostParticles>;

        SAMRAI::xfer::RefineAlgorithm Balgo;
        SAMRAI::xfer::RefineAlgorithm BalgoNode;
        SAMRAI::xfer::RefineAlgorithm BcopyAlgo;
        SAMRAI::xfer::RefineAlgorithm BregridAlgo;
        SAMRAI::xfer::RefineAlgorithm Ealgo;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magInitRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> elecPatchGhostsRefineSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> magSharedNodeRefineSchedules;

        SAMRAI::xfer::CoarsenAlgorithm RefluxAlgo{SAMRAI::tbox::Dimension{dimension}};
        SAMRAI::xfer::RefineAlgorithm PatchGhostRefluxedAlgo;
        std::map<int, std::shared_ptr<SAMRAI::xfer::CoarsenSchedule>> refluxSchedules;
        std::map<int, std::shared_ptr<SAMRAI::xfer::RefineSchedule>> patchGhostRefluxedSchedules;

        SharedNodeRefinerPool elecSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool elecGhostsRefiners_{resourcesManager_};
        GhostRefinerPool currentGhostsRefiners_{resourcesManager_};
        GhostRefinerPool rhoGhostsRefiners_{resourcesManager_};
        GhostRefinerPool velGhostsRefiners_{resourcesManager_};
        GhostRefinerPool pressureGhostsRefiners_{resourcesManager_};
        GhostRefinerPool momentumGhostsRefiners_{resourcesManager_};
        GhostRefinerPool totalEnergyGhostsRefiners_{resourcesManager_};
        SharedNodeRefinerPool magFluxesXSharedNodesRefiners_{resourcesManager_};
        SharedNodeRefinerPool magFluxesYSharedNodesRefiners_{resourcesManager_};
        SharedNodeRefinerPool magFluxesZSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesXGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesYGhostRefiners_{resourcesManager_};
        GhostRefinerPool magFluxesZGhostRefiners_{resourcesManager_};

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
        using BaseRefineOp          = FieldRefineOperator<GridLayoutT, GridT, Policy>;
        using DefaultFieldRefineOp  = BaseRefineOp<DefaultFieldRefiner<dimension>>;
        using MagneticFieldRefineOp = BaseRefineOp<MagneticFieldRefiner<dimension>>;
        using MagneticFieldRegridOp = BaseRefineOp<MagneticFieldRegrider<dimension>>;
        using ElectricFieldRefineOp = BaseRefineOp<ElectricFieldRefiner<dimension>>;
        using MHDFluxRefineOp       = BaseRefineOp<MHDFluxRefiner<dimension>>;
        using MHDFieldRefineOp      = BaseRefineOp<MHDFieldRefiner<dimension>>;
        using FieldTimeInterp       = FieldLinearTimeInterpolate<GridLayoutT, GridT>;

        template<typename Policy>
        using BaseCoarsenOp          = FieldCoarsenOperator<GridLayoutT, GridT, Policy>;
        using DefaultCoarsenOp       = BaseCoarsenOp<DefaultFieldCoarsener<dimension>>;
        using ElectricFieldCoarsenOp = BaseCoarsenOp<ElectricFieldCoarsener<dimension>>;
        using MHDFluxCoarsenOp       = BaseCoarsenOp<MHDFluxCoarsener<dimension>>;
        // using MagneticCoarsenOp      = BaseCoarsenOp<MagneticFieldCoarsener<dimension>>; //

        RefOp_ptr mhdFluxRefineOp_{std::make_shared<MHDFluxRefineOp>()};
        RefOp_ptr mhdFluxNodeRefineOp_{std::make_shared<MHDFieldRefineOp>(/*node_only=*/true)};

        RefOp_ptr mhdFieldRefineOp_{std::make_shared<MHDFieldRefineOp>()};

        RefOp_ptr EfieldRefineOp_{std::make_shared<ElectricFieldRefineOp>()};
        RefOp_ptr EfieldNodeRefineOp_{std::make_shared<ElectricFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr fieldNodeRefineOp_{std::make_shared<DefaultFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr fieldRefineOp_{std::make_shared<DefaultFieldRefineOp>()};

        TimeOp_ptr fieldTimeOp_{std::make_shared<FieldTimeInterp>()};

        CoarsenOp_ptr fieldCoarseningOp_{std::make_shared<DefaultCoarsenOp>()};

        // one refine operator per component
        RefOp_ptr mhdFluxRefineOp_x_{std::make_shared<MHDFluxRefineOp>()};
        RefOp_ptr mhdFluxRefineOp_y_{std::make_shared<MHDFluxRefineOp>()};
        RefOp_ptr mhdFluxRefineOp_z_{std::make_shared<MHDFluxRefineOp>()};

        RefOp_ptr BfieldRefineOp_x_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr BfieldRefineOp_y_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr BfieldRefineOp_z_{std::make_shared<MagneticFieldRefineOp>()};

        RefOp_ptr BfieldNodeRefineOp_x_{
            std::make_shared<MagneticFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr BfieldNodeRefineOp_y_{
            std::make_shared<MagneticFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr BfieldNodeRefineOp_z_{
            std::make_shared<MagneticFieldRefineOp>(/*node_only=*/true)};

        RefOp_ptr BfieldRegridOp_x_{std::make_shared<MagneticFieldRegridOp>()};
        RefOp_ptr BfieldRegridOp_y_{std::make_shared<MagneticFieldRegridOp>()};
        RefOp_ptr BfieldRegridOp_z_{std::make_shared<MagneticFieldRegridOp>()};

        RefOp_ptr EfieldRefineOp_x_{std::make_shared<ElectricFieldRefineOp>()};
        RefOp_ptr EfieldRefineOp_y_{std::make_shared<ElectricFieldRefineOp>()};
        RefOp_ptr EfieldRefineOp_z_{std::make_shared<ElectricFieldRefineOp>()};

        CoarsenOp_ptr mhdFluxCoarseningOp_x_{std::make_shared<MHDFluxCoarsenOp>()};
        CoarsenOp_ptr mhdFluxCoarseningOp_y_{std::make_shared<MHDFluxCoarsenOp>()};
        CoarsenOp_ptr mhdFluxCoarseningOp_z_{std::make_shared<MHDFluxCoarsenOp>()};

        CoarsenOp_ptr electricFieldCoarseningOp_x_{std::make_shared<ElectricFieldCoarsenOp>()};
        CoarsenOp_ptr electricFieldCoarseningOp_y_{std::make_shared<ElectricFieldCoarsenOp>()};
        CoarsenOp_ptr electricFieldCoarseningOp_z_{std::make_shared<ElectricFieldCoarsenOp>()};

        // CoarsenOp_ptr magneticCoarseningOp_{std::make_shared<MagneticCoarsenOp>()}; //

        MagneticRefinePatchStrategy<ResourcesManagerT, FieldDataT> magneticRefinePatchStrategy_{
            *resourcesManager_};
    };

} // namespace amr
} // namespace PHARE
#endif
