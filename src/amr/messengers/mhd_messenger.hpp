
#ifndef PHARE_MHD_MESSENGER_HPP
#define PHARE_MHD_MESSENGER_HPP

#include <memory>
#include <string>
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/field/coarsening/field_coarsen_operator.hpp"
#include "amr/data/field/coarsening/magnetic_field_coarsener.hpp"
#include "amr/data/field/refine/electric_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_refiner.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/messengers/refiner.hpp"
#include "amr/messengers/refiner_pool.hpp"
#include "amr/messengers/synchronizer_pool.hpp"
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/RefineOperator.h>

#include "core/mhd/mhd_quantities.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"

namespace PHARE
{
namespace amr
{
    template<typename MHDModel>
    class MHDMessenger : public IMessenger<typename MHDModel::Interface>
    {
        using IPhysicalModel            = typename MHDModel::Interface;
        using FieldT                    = typename MHDModel::field_type;
        using VecFieldT                 = typename MHDModel::vecfield_type;
        using MHDStateT                 = typename MHDModel::state_type;
        using GridLayoutT               = typename MHDModel::gridlayout_type;
        using GridT                     = typename MHDModel::grid_type;
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

            registerGhostComms_(mhdInfo);
            registerInitComms_(mhdInfo);
            registerSyncComms_(mhdInfo);
        }



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magSharedNodesRefiners_.registerLevel(hierarchy, level);
            magPatchGhostsRefiners_.registerLevel(hierarchy, level);
            magGhostsRefiners_.registerLevel(hierarchy, level);
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
                densityInitRefiners_.registerLevel(hierarchy, level);
                momentumInitRefiners_.registerLevel(hierarchy, level);
                magneticInitRefiners_.registerLevel(hierarchy, level);
                totalEnergyInitRefiners_.registerLevel(hierarchy, level);

                densitySynchronizers_.registerLevel(hierarchy, level);
                momentumSynchronizers_.registerLevel(hierarchy, level);
                magnetoSynchronizers_.registerLevel(hierarchy, level);
                totalEnergySynchronizers_.registerLevel(hierarchy, level);
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

            densityInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            momentumInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            magneticInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            totalEnergyInitRefiners_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);

            if (!isRegriddingL0)
            {
                auto& B = mhdModel.state.B;
                magGhostsRefiners_.fill(B, levelNumber, initDataTime);

                fix_magnetic_divergence_(*hierarchy, levelNumber, B);
            }
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

            densityInitRefiners_.fill(levelNumber, initDataTime);
            momentumInitRefiners_.fill(levelNumber, initDataTime);
            magneticInitRefiners_.fill(levelNumber, initDataTime);
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
            auto levelNumber = level.getLevelNumber();

            densitySynchronizers_.sync(levelNumber);
            momentumSynchronizers_.sync(levelNumber);
            magnetoSynchronizers_.sync(levelNumber);
            totalEnergySynchronizers_.sync(levelNumber);
        }

        void postSynchronize(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                             double const time) override
        {
            auto levelNumber = level.getLevelNumber();
            auto& mhdModel   = static_cast<MHDModel&>(model);

            magSharedNodesRefiners_.fill(mhdModel.state.B, levelNumber, time);
            magPatchGhostsRefiners_.fill(mhdModel.state.B, levelNumber, time);
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
        // Maybe we also need conservative ghost refiners for amr operations, actually quite likely
        void registerGhostComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& vecFieldNames) {
                std::vector<std::string> keys;
                std::transform(std::begin(vecFieldNames), std::end(vecFieldNames),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            magSharedNodesRefiners_.addStaticRefiners(info->ghostMagnetic, BfieldNodeRefineOp_,
                                                      makeKeys(info->ghostMagnetic));

            magGhostsRefiners_.addStaticRefiners(info->ghostMagnetic, BfieldRefineOp_,
                                                 makeKeys(info->ghostMagnetic));

            magPatchGhostsRefiners_.addStaticRefiner(info->modelMagnetic, BfieldRefineOp_,
                                                     info->modelMagnetic.vecName);

            elecSharedNodesRefiners_.addStaticRefiners(info->ghostElectric, EfieldNodeRefineOp_,
                                                       makeKeys(info->ghostElectric));

            elecGhostsRefiners_.addStaticRefiners(info->ghostElectric, EfieldRefineOp_,
                                                  makeKeys(info->ghostElectric));

            currentGhostsRefiners_.addTimeRefiners(info->ghostCurrent, info->modelCurrent,
                                                   core::VecFieldNames{Jold_}, EfieldRefineOp_,
                                                   fieldTimeOp_);

            rhoGhostsRefiners_.addTimeRefiners(info->ghostDensity, info->modelDensity,
                                               rhoOld_.name(), fieldRefineOp_, fieldTimeOp_);


            velGhostsRefiners_.addTimeRefiners(info->ghostVelocity, info->modelVelocity,
                                               core::VecFieldNames{Vold_}, fieldRefineOp_,
                                               fieldTimeOp_);

            pressureGhostsRefiners_.addTimeRefiners(info->ghostPressure, info->modelPressure,
                                                    Pold_.name(), fieldRefineOp_, fieldTimeOp_);

            momentumGhostsRefiners_.addTimeRefiners(info->ghostMomentum, info->modelMomentum,
                                                    core::VecFieldNames{rhoVold_}, fieldRefineOp_,
                                                    fieldTimeOp_);

            totalEnergyGhostsRefiners_.addTimeRefiners(info->ghostTotalEnergy,
                                                       info->modelTotalEnergy, EtotOld_.name(),
                                                       fieldRefineOp_, fieldTimeOp_);


            magFluxesXSharedNodesRefiners_.addStaticRefiners(info->ghostMagneticFluxesX,
                                                             fieldNodeRefineOp_,
                                                             makeKeys(info->ghostMagneticFluxesX));

            magFluxesYSharedNodesRefiners_.addStaticRefiners(info->ghostMagneticFluxesY,
                                                             fieldNodeRefineOp_,
                                                             makeKeys(info->ghostMagneticFluxesY));

            magFluxesZSharedNodesRefiners_.addStaticRefiners(info->ghostMagneticFluxesZ,
                                                             fieldNodeRefineOp_,
                                                             makeKeys(info->ghostMagneticFluxesZ));

            magFluxesXGhostRefiners_.addStaticRefiners(info->ghostMagneticFluxesX, fieldRefineOp_,
                                                       makeKeys(info->ghostMagneticFluxesX));

            magFluxesYGhostRefiners_.addStaticRefiners(info->ghostMagneticFluxesY, fieldRefineOp_,
                                                       makeKeys(info->ghostMagneticFluxesY));

            magFluxesZGhostRefiners_.addStaticRefiners(info->ghostMagneticFluxesZ, fieldRefineOp_,
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

            densityInitRefiners_.addStaticRefiners(info->initDensity, fieldRefineOp_,
                                                   info->initDensity);

            momentumInitRefiners_.addStaticRefiners(info->initMomentum, fieldRefineOp_,
                                                    makeKeys(info->initMomentum));

            magneticInitRefiners_.addStaticRefiners(info->initMagnetic, BfieldRefineOp_,
                                                    makeKeys(info->initMagnetic));

            totalEnergyInitRefiners_.addStaticRefiners(info->initTotalEnergy, fieldRefineOp_,
                                                       info->initTotalEnergy);
        }




        // should this use conservative quantities ? Quite likely as these are what we store on the
        // grid
        void registerSyncComms_(std::unique_ptr<MHDMessengerInfo> const& info)
        {
            densitySynchronizers_.add(info->modelDensity, fieldCoarseningOp_, info->modelDensity);

            momentumSynchronizers_.add(info->modelMomentum, fieldCoarseningOp_,
                                       info->modelMomentum.vecName);

            magnetoSynchronizers_.add(info->modelMagnetic, magneticCoarseningOp_,
                                      info->modelMagnetic.vecName);

            totalEnergySynchronizers_.add(info->modelTotalEnergy, fieldCoarseningOp_,
                                          info->modelTotalEnergy);
        }

        void debug_print(VecFieldT const& B, GridLayoutT const& layout, int loc, int ix, int iy,
                         std::string const& aftbef)
        {
            auto& Bx       = B(core::Component::X);
            auto& By       = B(core::Component::Y);
            auto& Bz       = B(core::Component::Z);
            auto const& dx = layout.meshSize()[0];
            auto const& dy = layout.meshSize()[1];

            if (loc == 3) // w hi, y hi
            {
                std::cout << aftbef << "\n";
                std::cout << "cell 4 : "
                          << (Bx(ix, iy) - Bx(ix - 1, iy)) / dx
                                 + (By(ix - 1, iy + 1) - By(ix - 1, iy)) / dy
                          << "\n";

                std::cout << "cell 2 : "
                          << (Bx(ix + 1, iy - 1) - Bx(ix, iy - 1)) / dx
                                 + (By(ix, iy) - By(ix, iy - 1)) / dy
                          << "\n";

                std::cout << "cell 1 : "
                          << (Bx(ix + 1, iy) - Bx(ix, iy)) / dx + (By(ix, iy + 1) - By(ix, iy)) / dy
                          << "\n";

                std::cout << "cell 3 : "
                          << (Bx(ix, iy - 1) - Bx(ix - 1, iy - 1)) / dx
                                 + (By(ix - 1, iy) - By(ix - 1, iy - 1)) / dy
                          << "\n";
            }
        }

        /*
         *
         * */
        void fix_magnetic_divergence_(SAMRAI::hier::PatchHierarchy const& hierarchy,
                                      int levelNumber, VecFieldT& B)
        {
            auto lvlBoundary = SAMRAI::hier::CoarseFineBoundary(
                hierarchy, levelNumber,
                SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension}, 0});

            for (auto& patch : *hierarchy.getPatchLevel(levelNumber))
            {
                if constexpr (dimension == 2)
                {
                    auto _         = resourcesManager_->setOnPatch(*patch, B);
                    auto layout    = layoutFromPatch<GridLayoutT>(*patch);
                    auto& Bx       = B(core::Component::X);
                    auto& By       = B(core::Component::Y);
                    auto& Bz       = B(core::Component::Z);
                    auto const& dx = layout.meshSize()[0];
                    auto const& dy = layout.meshSize()[1];

                    auto boundaries = lvlBoundary.getEdgeBoundaries(patch->getGlobalId());
                    for (auto& boundary : boundaries)
                    {
                        int loc         = boundary.getLocationIndex();
                        auto const& box = boundary.getBox();

                        if (loc == 0) // x_lo
                        {
                            // we're on the left side border outside the domain
                            // we need to get to the first domain cell
                            auto ixAMR = box.lower()[0] + 1;

                            // this is a X edge we need to fix By
                            for (int iyAMR = box.lower()[1]; iyAMR <= box.upper()[1]; ++iyAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(iyAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    By(ix, iy)
                                        = By(ix, iy + 1) + dy / dx * (Bx(ix + 1, iy) - Bx(ix, iy));
                                }
                            }
                        }
                        else if (loc == 1) // x_hi
                        {
                            // we're on the right side border outside the domain
                            // we need to get to the last domain cell
                            auto ixAMR = box.upper()[0] - 1;

                            // this is a X edge we need to fix By
                            for (int iyAMR = box.lower()[1]; iyAMR <= box.upper()[1]; ++iyAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(iyAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    By(ix, iy)
                                        = By(ix, iy + 1) + dy / dx * (Bx(ix + 1, iy) - Bx(ix, iy));
                                }
                            }
                        }
                        else if (loc == 2) // y_lo
                        {
                            // we're on the bottom edge, we need the first domain cell
                            auto iyAMR = box.lower()[1] + 1;

                            // this is a Y edge we need to fix Bx
                            for (int ixAMR = box.lower()[0]; ixAMR <= box.upper()[0]; ++ixAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(ixAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    Bx(ix, iy)
                                        = Bx(ix + 1, iy) + dx / dy * (By(ix, iy + 1) - By(ix, iy));
                                }
                            }
                        }
                        else if (loc == 3) // y_hi
                        {
                            // we're on the top edge, we need the last domain cell
                            auto iyAMR = box.upper()[1] - 1;

                            // this is a Y edge we need to fix Bx
                            for (int ixAMR = box.lower()[0]; ixAMR <= box.upper()[0]; ++ixAMR)
                            {
                                // we want to change By only at fine faces not shared with
                                // coarse faces.
                                if (!(ixAMR % 2 == 0))
                                {
                                    auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                                    auto ix       = localIdx[0];
                                    auto iy       = localIdx[1];
                                    Bx(ix, iy)
                                        = Bx(ix + 1, iy) + dx / dy * (By(ix, iy + 1) - By(ix, iy));
                                }
                            }
                        }
                    } // end boundary loop

                    // above we have treated boundaries as 1D lines
                    // this is not treating corners well and we need to deal with them separatly
                    //
                    //__________________
                    //      |    |     |
                    //      | 4  BX 1  |
                    //      |    |     |
                    //      |_By_|_BY___
                    //      |    |     |
                    //      | 3 Bx  2  |
                    //      |    |     |
                    //      |____|_____
                    //                 |
                    //                 |
                    //                 |
                    //
                    // above are the 4 fine top right corner cells and
                    // above code have changed BX and BY but not Bx and By
                    // faces. But these cells are not independant.
                    // To fix divB in these 4 cells we will re-assigne the four
                    // fine faces.
                    // The idea is to give them the same values they would have had
                    // if refined. By and BY would have the same value, equal to the
                    // average of the coarse By on top and bottom faces of the coarse cell
                    // we do not have access to the coarse cell here but we know that because
                    // of the previous coarsening, these coarse faces have the same flux
                    // as the average of the 2 shared fine faces
                    // So By and BY will be 1/4 * sum of top and bottom fine By
                    // Similarly, BX and Bx will take 1/4 the four fine Bx faces share with the
                    // coarse ones
                    //
                    auto corners = lvlBoundary.getNodeBoundaries(patch->getGlobalId());
                    for (auto& corner : corners)
                    {
                        int loc = corner.getLocationIndex();

                        // the box we get should consist of just 1 cell
                        // i.e. with lower==upper so it should not matter which
                        // one we take in the following.
                        auto const& box = corner.getBox();

                        //* x_lo, y_lo: 0
                        //* x_hi, y_lo: 1
                        //* x_lo, y_hi: 2
                        // * x_hi, y_hi: 3

                        if (loc == 3) // x_hi, y_hi
                        {
                            // we're on the top right corner
                            // and we want the domain cell
                            // that is labeled cell 1 in above drawing
                            // we fix BX and BY
                            auto ixAMR    = box.lower()[0] - 1;
                            auto iyAMR    = box.lower()[1] - 1;
                            auto localIdx = layout.AMRToLocal(core::Point{ixAMR, iyAMR});
                            auto ix       = localIdx[0];
                            auto iy       = localIdx[1];

                            // maybe we should keep these for some time
                            // as comments in case they are useful again
                            PHARE_DEBUG_DO(std::string const before = "BEFORE";
                                           debug_print(B, layout, loc, ix, iy, before);)
                            Bx(ix, iy - 1)
                                = Bx(ix + 1, iy - 1) + dx / dy * (By(ix, iy) - By(ix, iy - 1));

                            By(ix - 1, iy)
                                = By(ix - 1, iy + 1) + dy / dx * (Bx(ix, iy) - Bx(ix - 1, iy));

                            PHARE_DEBUG_DO(std::string const after = "AFTER";
                                           debug_print(B, layout, loc, ix, iy, after);)
                        }
                    } // end corner loops
                } // end if 2D
            } // end patch loop
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

        SharedNodeRefinerPool magSharedNodesRefiners_{resourcesManager_};
        GhostRefinerPool magGhostsRefiners_{resourcesManager_};
        GhostRefinerPool magPatchGhostsRefiners_{resourcesManager_};
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
        InitRefinerPool magneticInitRefiners_{resourcesManager_};
        InitRefinerPool totalEnergyInitRefiners_{resourcesManager_};

        SynchronizerPool<rm_t> densitySynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> momentumSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> magnetoSynchronizers_{resourcesManager_};
        SynchronizerPool<rm_t> totalEnergySynchronizers_{resourcesManager_};

        using RefOp_ptr     = std::shared_ptr<SAMRAI::hier::RefineOperator>;
        using CoarsenOp_ptr = std::shared_ptr<SAMRAI::hier::CoarsenOperator>;
        using TimeOp_ptr    = std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator>;

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

        RefOp_ptr BfieldNodeRefineOp_{std::make_shared<MagneticFieldRefineOp>(/*node_only=*/
                                                                              true)};
        RefOp_ptr BfieldRefineOp_{std::make_shared<MagneticFieldRefineOp>()};
        RefOp_ptr EfieldNodeRefineOp_{std::make_shared<ElectricFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr EfieldRefineOp_{std::make_shared<ElectricFieldRefineOp>()};
        RefOp_ptr fieldNodeRefineOp_{std::make_shared<DefaultFieldRefineOp>(/*node_only=*/true)};
        RefOp_ptr fieldRefineOp_{std::make_shared<DefaultFieldRefineOp>()};

        TimeOp_ptr fieldTimeOp_{std::make_shared<FieldTimeInterp>()};

        CoarsenOp_ptr fieldCoarseningOp_{std::make_shared<DefaultCoarsenOp>()};
        CoarsenOp_ptr magneticCoarseningOp_{std::make_shared<MagneticCoarsenOp>()};
    };

} // namespace amr
} // namespace PHARE
#endif
