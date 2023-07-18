#ifndef PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_HYBRID_HYBRID_MESSENGER_STRATEGY_HPP

#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/IntVector.h"
#include "refiner_pool.hpp"
#include "synchronizer_pool.hpp"
#include "amr/data/field/coarsening/default_field_coarsener.hpp"
#include "amr/data/field/coarsening/magnetic_field_coarsener.hpp"
#include "amr/data/field/refine/field_refiner.hpp"
#include "amr/data/field/refine/magnetic_field_refiner.hpp"
#include "amr/data/field/refine/electric_field_refiner.hpp"
#include "amr/data/field/time_interpolate/field_linear_time_interpolate.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "core/numerics/interpolator/interpolator.hpp"
#include "core/hybrid/hybrid_quantities.hpp"



#include <SAMRAI/xfer/RefineAlgorithm.h>
#include <SAMRAI/xfer/RefineSchedule.h>


#include <iterator>
#include <optional>
#include <utility>
#include <iomanip>


namespace PHARE
{
namespace amr
{

    // this structure is a wrapper of a field*
    // so to serve as a ResourcesUser for the ResourcesManager
    template<typename FieldT>
    struct FieldUser
    {
        struct Property
        {
            std::string name;
            typename HybridQuantity::Scalar qty;
        };
        FieldUser(std::string fieldName, FieldT* ptr, typename HybridQuantity::Scalar qty)
            : name{fieldName}
            , f{ptr}
            , quantity{qty}
        {
        }
        std::string name;
        FieldT* f;
        typename HybridQuantity::Scalar quantity;
        using field_type = FieldT;

        std::vector<Property> getFieldNamesAndQuantities() const { return {{name, quantity}}; }

        void setBuffer(std::string const& /*bufferName*/, FieldT* field) { f = field; }
        void copyData(FieldT const& source) { f->copyData(source); }
    };



    /** \brief An HybridMessenger is the specialization of a HybridMessengerStrategy for hybrid to
     * hybrid data communications.
     */
    template<typename HybridModel, typename RefinementParams>
    class HybridHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
    {
        using IonsT                              = typename HybridModel::ions_type;
        using ElectromagT                        = typename HybridModel::electromag_type;
        using VecFieldT                          = typename HybridModel::vecfield_type;
        using GridLayoutT                        = typename HybridModel::gridlayout_type;
        using FieldT                             = typename VecFieldT::field_type;
        using ResourcesManagerT                  = typename HybridModel::resources_manager_type;
        static constexpr std::size_t dimension   = GridLayoutT::dimension;
        static constexpr std::size_t interpOrder = GridLayoutT::interp_order;
        using IPhysicalModel                     = typename HybridModel::Interface;

        using InteriorParticleRefineOp = typename RefinementParams::InteriorParticleRefineOp;
        using CoarseToFineRefineOpOld  = typename RefinementParams::CoarseToFineRefineOpOld;
        using CoarseToFineRefineOpNew  = typename RefinementParams::CoarseToFineRefineOpNew;


    public:
        static const std::string stratName;
        static constexpr std::size_t rootLevelNumber = 0;


        HybridHybridMessengerStrategy(std::shared_ptr<ResourcesManagerT> manager,
                                      int const firstLevel)
            : HybridMessengerStrategy<HybridModel>{stratName}
            , resourcesManager_{std::move(manager)}
            , firstLevel_{firstLevel}
        {
            // resourcesManager_->registerResources(EM_old_);
            resourcesManager_->registerResources(Jold_);
            resourcesManager_->registerResources(NiOldUser_);
            resourcesManager_->registerResources(ViOld_);
        }

        virtual ~HybridHybridMessengerStrategy() = default;



        /* ------------------------------------------------------------------------
                        methods used for the IMessenger interface
           ------------------------------------------------------------------------ */


        /**
         * @brief allocate the messenger strategy internal variables to the model resourceManager
         */
        void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            // resourcesManager_->allocate(Eavg_, patch, allocateTime);
            resourcesManager_->allocate(Jold_, patch, allocateTime);
            resourcesManager_->allocate(NiOldUser_, patch, allocateTime);
            resourcesManager_->allocate(ViOld_, patch, allocateTime);
        }



        /**
         * @brief setup creates all SAMRAI algorithms to communicate data involved in a messenger
         * between the coarse and fine levels.
         *
         * This method creates the SAMRAI algorithms for communications associated between pairs of
         * variables. The function does not create the SAMRAI schedules since they depend on the
         * levels
         */
        void registerQuantities([[maybe_unused]] std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<HybridMessengerInfo> hybridInfo{
                dynamic_cast<HybridMessengerInfo*>(fromFinerInfo.release())};

            registerGhostComms_(hybridInfo);
            registerInitComms(hybridInfo);
            registerSyncComms(hybridInfo);
        }



        /**
         * @brief registerLevel registers the level for all Communicators
         *
         * The level must always be registered to ghost Communicators but only to init Communicators
         * if it is not the root level since root is not initialized by coarser data.
         *
         * Ghost communicators that need to know this level :
         *
         *  - magnetic fields
         *  - electric fields
         *  - patch ghost particles
         *
         *  ion moments :
         *  do not need to be filled on shared border nodes by SAMRAI schedules
         *  since they will be filled with levelGhostParticles[old,new] on level ghost nodes
         *  and computed by ghost particles on interior patch ghost nodes
         *  however they need to be filled on pure ghost nodes.
         *  Note : they are filled for all pure ghost nodes for now although
         *  only the first one is needed (for the coarsening stencil).
         *
         * Init communicators that need to know this level :
         *
         *  - magnetic fields
         *  - electric fields
         *  - ion interior particle arrays
         *  - ion levelGhostParticlesOld particle arrays
         */
        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                           int const levelNumber) override
        {
            auto const level = hierarchy->getPatchLevel(levelNumber);

            magneticSharedNodes_.registerLevel(hierarchy, level);
            electricSharedNodes_.registerLevel(hierarchy, level);
            currentSharedNodes_.registerLevel(hierarchy, level);

            magneticPatchGhosts_.registerLevel(hierarchy, level);
            // magneticLevelGhosts_.registerLevel(hierarchy, level);
            magneticGhosts_.registerLevel(hierarchy, level);
            electricGhosts_.registerLevel(hierarchy, level);
            currentGhosts_.registerLevel(hierarchy, level);

            densityGhosts_.registerLevel(hierarchy, level);
            bulkVelGhosts_.registerLevel(hierarchy, level);

            patchGhostParticles_.registerLevel(hierarchy, level);

            // root level is not initialized with a schedule using coarser level data
            // so we don't create these schedules if root level
            if (levelNumber != rootLevelNumber)
            {
                // those are for refinement
                magneticInit_.registerLevel(hierarchy, level);
                electricInit_.registerLevel(hierarchy, level);
                interiorParticles_.registerLevel(hierarchy, level);
                levelGhostParticlesOld_.registerLevel(hierarchy, level);
                levelGhostParticlesNew_.registerLevel(hierarchy, level);

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
                    const int levelNumber,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                    IPhysicalModel& model, double const initDataTime) override
        {
            auto& hybridModel = dynamic_cast<HybridModel&>(model);
            auto level        = hierarchy->getPatchLevel(levelNumber);
            magneticInit_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            electricInit_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            interiorParticles_.regrid(hierarchy, levelNumber, oldLevel, initDataTime);
            patchGhostParticles_.fill(levelNumber, initDataTime);

            // regriding will fill the new level wherever it has points that overlap
            // old level. This will include its level border points.
            // These new level border points will thus take values that where previous
            // domain values. Magnetic flux is thus not necessarily consistent with
            // the Loring et al. method to sync the induction between coarse and fine faces.
            // Specifically, we need all fine faces to have equal magnetic field and also
            // equal to that of the shared coarse face.
            // This means that we now need to fill ghosts and border included
            auto& B = hybridModel.state.electromag.B;
            auto& E = hybridModel.state.electromag.E;
            // magneticSharedNodes_.fill(B, levelNumber, initDataTime);
            magneticGhosts_.fill(B, levelNumber, initDataTime);
            // electricSharedNodes_.fill(E, levelNumber, initDataTime);
            electricGhosts_.fill(E, levelNumber, initDataTime);

            auto lvlBoundary = SAMRAI::hier::CoarseFineBoundary(
                *hierarchy, levelNumber,
                SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{dimension}, 0});

            for (auto& patch : *hierarchy->getPatchLevel(levelNumber))
            {
                if constexpr (dimension == 2)
                {
                    auto _         = resourcesManager_->setOnPatch(*patch, B);
                    auto layout    = layoutFromPatch<GridLayoutT>(*patch);
                    auto& Bx       = B(Component::X);
                    auto& By       = B(Component::Y);
                    auto& Bz       = B(Component::Z);
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
                                    auto localIdx = layout.AMRToLocal(Point{ixAMR, iyAMR});
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
                                    auto localIdx = layout.AMRToLocal(Point{ixAMR, iyAMR});
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
                                    auto localIdx = layout.AMRToLocal(Point{ixAMR, iyAMR});
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
                                    auto localIdx = layout.AMRToLocal(Point{ixAMR, iyAMR});
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
                            auto localIdx = layout.AMRToLocal(Point{ixAMR, iyAMR});
                            auto ix       = localIdx[0];
                            auto iy       = localIdx[1];

                            std::cout << "BEFORE\n";
                            std::cout << "cell 4 : "
                                      << (Bx(ix, iy) - Bx(ix - 1, iy)) / dx
                                             + (By(ix - 1, iy + 1) - By(ix - 1, iy)) / dy
                                      << "\n";

                            std::cout << "cell 2 : "
                                      << (Bx(ix + 1, iy - 1) - Bx(ix, iy - 1)) / dx
                                             + (By(ix, iy) - By(ix, iy - 1)) / dy
                                      << "\n";

                            std::cout << "cell 1 : "
                                      << (Bx(ix + 1, iy) - Bx(ix, iy)) / dx
                                             + (By(ix, iy + 1) - By(ix, iy)) / dy
                                      << "\n";

                            std::cout << "cell 3 : "
                                      << (Bx(ix, iy - 1) - Bx(ix - 1, iy - 1)) / dx
                                             + (By(ix - 1, iy) - By(ix - 1, iy - 1)) / dy
                                      << "\n";

                            Bx(ix, iy - 1)
                                = Bx(ix + 1, iy - 1) + dx / dy * (By(ix, iy) - By(ix, iy - 1));

                            By(ix - 1, iy)
                                = By(ix - 1, iy + 1) + dy / dx * (Bx(ix, iy) - Bx(ix - 1, iy));


                            std::cout << "AFTER";
                            std::cout << "cell 4 : "
                                      << (Bx(ix, iy) - Bx(ix - 1, iy)) / dx
                                             + (By(ix - 1, iy + 1) - By(ix - 1, iy)) / dy
                                      << "\n";

                            std::cout << "cell 2 : "
                                      << (Bx(ix + 1, iy - 1) - Bx(ix, iy - 1)) / dx
                                             + (By(ix, iy) - By(ix, iy - 1)) / dy
                                      << "\n";

                            std::cout << "cell 1 : "
                                      << (Bx(ix + 1, iy) - Bx(ix, iy)) / dx
                                             + (By(ix, iy + 1) - By(ix, iy)) / dy
                                      << "\n";

                            std::cout << "cell 3 : "
                                      << (Bx(ix, iy - 1) - Bx(ix - 1, iy - 1)) / dx
                                             + (By(ix - 1, iy) - By(ix - 1, iy - 1)) / dy
                                      << "\n";

                            /*
                                                        Bx(ix, iy) = 0.25
                                                                     * (Bx(ix + 1, iy) + Bx(ix + 1,
                               iy - 1) + Bx(ix - 1, iy)
                                                                        + Bx(ix - 1, iy - 1));

                                                        By(ix, iy) = 0.25
                                                                     * (By(ix, iy + 1) + By(ix - 1,
                               iy + 1) + By(ix - 1, iy - 1)
                                                                        + By(ix, iy - 1));
                            */



                            // now we need to fix Bx between cells 2 and 3
                            // and By between 3 and 4.
                            // Bx should be strictly equal to BX
                            // and By should be strictly equal to BY
                            //                           Bx(ix, iy - 1) = Bx(ix, iy);
                            //                            By(ix - 1, iy) = By(ix, iy);
                        }
                    } // end corner loops
                }     // end if 2D
            }         // end patch loop




            // we now call only levelGhostParticlesOld.fill() and not .regrid()
            // regrid() would refine from next coarser in regions of level not overlaping
            // oldLevel, but copy from domain particles of oldLevel where there is an
            // overlap while we do not a priori see why this could be wrong,but this led to
            // occasional failures of the SAMRAI MPI module. See
            // https://github.com/PHAREHUB/PHARE/issues/604 calling .fill() ensures that
            // levelGhostParticlesOld particles are filled exclusively from spliting next
            // coarser domain ones like when a new finest level is created.
            levelGhostParticlesOld_.fill(levelNumber, initDataTime);
            copyLevelGhostOldToPushable_(*level, model);

            // computeIonMoments_(*level, model);
            // levelGhostNew will be refined in next firstStep
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

            magneticInit_.fill(levelNumber, initDataTime);
            electricInit_.fill(levelNumber, initDataTime);

            // no need to call these :
            // magneticGhosts_.fill(levelNumber, initDataTime);
            // electricGhosts_.fill(levelNumber, initDataTime);
            // because the SAMRAI schedules in the 'init' communicators
            // already fill the patch ghost box from the neighbor interior box.
            // so ghost nodes are already filled .

            PHARE_LOG_START("hybhybmessengerStrat::initLevel : interior part fill schedule");
            interiorParticles_.fill(levelNumber, initDataTime);
            PHARE_LOG_STOP("hybhybmessengerStrat::initLevel : interior part fill schedule");
            // however we need to call the ghost communicator for patch ghost particles
            // since the interior schedules have a restriction to the interior of the patch.
            PHARE_LOG_START("hybhybmessengerStrat::initLevel : patch ghost part fill schedule");
            patchGhostParticles_.fill(levelNumber, initDataTime);
            PHARE_LOG_STOP("hybhybmessengerStrat::initLevel : patch ghost part fill schedule");


            levelGhostParticlesOld_.fill(levelNumber, initDataTime);


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
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillElectricGhosts");
            std::cout << "filling electric ghosts\n";
            electricSharedNodes_.fill(E, levelNumber, fillTime);
            electricGhosts_.fill(E, levelNumber, fillTime);
        }




        void fillCurrentGhosts(VecFieldT& J, int const levelNumber, double const fillTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillCurrentGhosts");
            currentSharedNodes_.fill(J, levelNumber, fillTime);
            currentGhosts_.fill(J, levelNumber, fillTime);
        }




        /**
         * @brief fillIonGhostParticles will fill the interior ghost particle array from
         * neighbor patches of the same level. Before doing that, it empties the array for
         * all populations
         */
        void fillIonGhostParticles(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                   double const fillTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillIonGhostParticles");

            for (auto patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(*patch, ions);
                for (auto& pop : ions)
                {
                    empty(pop.patchGhostParticles());
                }
            }

            patchGhostParticles_.fill(level.getLevelNumber(), fillTime);
        }




        /**
         * @brief fillIonMomentGhosts works on moment ghost nodes
         *
         * patch border node moments are completed by the deposition of patch ghost
         * particles for all populations level border nodes are completed by the deposition
         * of level ghost [old,new] particles for all populations, linear time interpolation
         * is used to get the contribution of old/new particles
         */
        void fillIonPopMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                    double const afterPushTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::fillIonMomentGhosts");

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
            densityGhosts_.fill(level.getLevelNumber(), afterPushTime);
            bulkVelGhosts_.fill(level.getLevelNumber(), afterPushTime);
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
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::firstStep");

            auto levelNumber = level.getLevelNumber();
            if (newCoarserTime < prevCoarserTime)
                throw std::runtime_error(
                    "Error : prevCoarserTime (" + std::to_string(prevCoarserTime)
                    + ") should be < newCoarserTime (" + std::to_string(prevCoarserTime) + ")");

            // root level has no levelghost particles
            if (levelNumber != 0)
            {
                PHARE_LOG_START("HybridHybridMessengerStrategy::firstStep.fill");
                levelGhostParticlesNew_.fill(levelNumber, currentTime);
                PHARE_LOG_STOP("HybridHybridMessengerStrategy::firstStep.fill");

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
            if (level.getLevelNumber() > 0)
            {
                PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::lastStep");

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

                        core::swap(levelGhostParticlesNew, levelGhostParticlesOld);
                        core::empty(levelGhostParticlesNew);
                        core::empty(levelGhostParticles);
                        std::copy(std::begin(levelGhostParticlesOld),
                                  std::end(levelGhostParticlesOld),
                                  std::back_inserter(levelGhostParticles));

                        if (level.getLevelNumber() == 0)
                        {
                            if (levelGhostParticlesNew.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesNew detected in root level : "
                                    + std::to_string(levelGhostParticlesNew.size()));
                            if (levelGhostParticles.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticles detected in root level : "
                                    + std::to_string(levelGhostParticles.size()));
                            if (levelGhostParticlesOld.size() != 0)
                                throw std::runtime_error(
                                    "levelGhostParticlesOld detected in root level : "
                                    + std::to_string(levelGhostParticlesOld.size()));
                        }
                    }
                }
            }
        }



        /**
         * @brief prepareStep is the concrete implementation of the
         * HybridMessengerStrategy::prepareStep method For hybrid-Hybrid communications.
         * This method copies the current model electromagnetic field E,B, the current
         * density J, and the density and bulk velocity, defined at t=n. Since prepareStep()
         * is called just before advancing the level, this operation actually saves the t=n
         * versions of E,B,J,Ni,Vi into the messenger. When the time comes that the next
         * finer level needs to time interpolate the electromagnetic field and current at
         * its ghost nodes, this level will be able to interpolate at required time because
         * the t=n Vi,Ni,E,B,J fields of previous next coarser step will be in the
         * messenger.
         */
        void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                         double currentTime) override
        {
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::prepareStep");

            auto& hybridModel = static_cast<HybridModel&>(model);
            for (auto& patch : level)
            {
                auto dataOnPatch = resourcesManager_->setOnPatch(
                    *patch, hybridModel.state.electromag, hybridModel.state.J,
                    hybridModel.state.ions, Jold_, NiOldUser_, ViOld_);

                resourcesManager_->setTime(Jold_, *patch, currentTime);
                resourcesManager_->setTime(NiOldUser_, *patch, currentTime);
                resourcesManager_->setTime(ViOld_, *patch, currentTime);

                auto& EM = hybridModel.state.electromag;
                auto& J  = hybridModel.state.J;
                auto& Vi = hybridModel.state.ions.velocity();
                auto& Ni = hybridModel.state.ions.density();

                Jold_.copyData(J);
                ViOld_.copyData(Vi);
                NiOldUser_.copyData(Ni);
            }
        }




        void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                            double const initDataTime) override
        {
            auto levelNumber = level.getLevelNumber();
            assert(levelNumber == 0);

            auto& hybridModel = static_cast<HybridModel&>(model);

            electricSharedNodes_.fill(hybridModel.state.electromag.E, levelNumber, initDataTime);

            electricGhosts_.fill(hybridModel.state.electromag.E, levelNumber, initDataTime);
            patchGhostParticles_.fill(levelNumber, initDataTime);

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
            PHARE_LOG_SCOPE("HybridHybridMessengerStrategy::synchronize");

            auto levelNumber = level.getLevelNumber();
            std::cout << "synchronizing level " << levelNumber << "\n";

            // call coarsning schedules...
            magnetoSynchronizers_.sync(levelNumber);
            // electroSynchronizers_.sync(levelNumber);
            densitySynchronizers_.sync(levelNumber);
            ionBulkVelSynchronizers_.sync(levelNumber);
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

            std::cout << "postSynchronize level " << levelNumber << "\n";

            magneticSharedNodes_.fill(hybridModel.state.electromag.B, levelNumber, time);
            electricSharedNodes_.fill(hybridModel.state.electromag.E, levelNumber, time);

            // we fill magnetic field ghosts only on patch ghost nodes and not on level
            // ghosts the reason is that 1/ filling ghosts is necessary to prevent mismatch
            // between ghost and overlaped neighboring patch domain nodes resulting from
            // former coarsening which does not occur for level ghosts and 2/ overwriting
            // level border with next coarser model B would invalidate divB on the first
            // fine domain cell since its border face only received a fraction of the
            // induction that has occured on the shared coarse face.
            magneticPatchGhosts_.fill(hybridModel.state.electromag.B, levelNumber, time);
            electricGhosts_.fill(hybridModel.state.electromag.E, levelNumber, time);
            densityGhosts_.fill(levelNumber, time);
            bulkVelGhosts_.fill(hybridModel.state.ions.velocity(), levelNumber, time);
        }

    private:
        void registerGhostComms_(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& descriptor) {
                std::vector<std::string> keys;
                std::transform(std::begin(descriptor), std::end(descriptor),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };
            fillStaticRefiners_(info->ghostMagnetic, BfieldNodeRefineOp_, magneticSharedNodes_,
                                makeKeys(info->ghostMagnetic), /*ghosts=*/true);

            fillStaticRefiners_(info->ghostMagnetic, BfieldRefineOp_, magneticGhosts_,
                                makeKeys(info->ghostMagnetic), /*ghosts=*/true);

            fillStaticRefiners_(info->modelMagnetic, BfieldRefineOp_, magneticPatchGhosts_,
                                info->modelMagnetic.vecName, /*ghosts=*/true);

            fillStaticRefiners_(info->ghostElectric, EfieldNodeRefineOp_, electricSharedNodes_,
                                makeKeys(info->ghostElectric), /*ghosts=*/true);

            fillStaticRefiners_(info->ghostElectric, EfieldRefineOp_, electricGhosts_,
                                makeKeys(info->ghostElectric), /*ghosts=*/true);

            fillTimeRefiners_(info->ghostCurrent, info->modelCurrent, VecFieldDescriptor{Jold_},
                              currentSharedNodes_, EfieldNodeRefineOp_);

            fillTimeRefiners_(info->ghostCurrent, info->modelCurrent, VecFieldDescriptor{Jold_},
                              currentGhosts_, EfieldRefineOp_);

            densityGhosts_.addTimeRefiner(info->modelIonDensity, info->modelIonDensity,
                                          NiOldUser_.name, resourcesManager_, fieldRefineOp_,
                                          fieldTimeOp_, info->modelIonDensity);

            fillTimeRefiners_(info->ghostBulkVelocity, info->modelIonBulkVelocity,
                              VecFieldDescriptor{ViOld_}, bulkVelGhosts_);
        }




        void registerInitComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            auto makeKeys = [](auto const& descriptor) {
                std::vector<std::string> keys;
                std::transform(std::begin(descriptor), std::end(descriptor),
                               std::back_inserter(keys), [](auto const& d) { return d.vecName; });
                return keys;
            };

            fillStaticRefiners_(info->initMagnetic, BfieldRefineOp_, magneticInit_,
                                makeKeys(info->initMagnetic));

            fillStaticRefiners_(info->initElectric, EfieldRefineOp_, electricInit_,
                                makeKeys(info->initElectric));


            fillStaticRefiners_(info->interiorParticles, interiorParticleRefineOp_,
                                interiorParticles_, info->interiorParticles);


            fillStaticRefiners_(info->levelGhostParticlesOld, levelGhostParticlesOldOp_,
                                levelGhostParticlesOld_, info->levelGhostParticlesOld);


            fillStaticRefiners_(info->levelGhostParticlesNew, levelGhostParticlesNewOp_,
                                levelGhostParticlesNew_, info->levelGhostParticlesNew);


            fillStaticRefiners_(info->patchGhostParticles, nullptr, patchGhostParticles_,
                                info->patchGhostParticles);
        }




        void registerSyncComms(std::unique_ptr<HybridMessengerInfo> const& info)
        {
            magnetoSynchronizers_.add(info->modelMagnetic, resourcesManager_, magneticCoarseningOp_,
                                      info->modelMagnetic.vecName);

            electroSynchronizers_.add(info->modelElectric, resourcesManager_, fieldCoarseningOp_,
                                      info->modelElectric.vecName);

            ionBulkVelSynchronizers_.add(info->modelIonBulkVelocity, resourcesManager_,
                                         fieldCoarseningOp_, info->modelIonBulkVelocity.vecName);


            densitySynchronizers_.add(info->modelIonDensity, fieldCoarseningOp_,
                                      info->modelIonDensity, resourcesManager_);
        }


        /**
         * @brief fill the given pool of refiners with a new refiner per VecFieldDescriptor
         * in ghostVecs. Data will be spatially refined using the specified refinement
         * operator, and time interpolated between time n and n+1 of next coarser data,
         * represented by modelVec and oldModelVec descriptors.
         */
        template<typename RefinerT, RefinerT RefineType>
        void fillTimeRefiners_(std::vector<VecFieldDescriptor> const& ghostVecs,
                               VecFieldDescriptor const& modelVec,
                               VecFieldDescriptor const& oldModelVec,
                               RefinerPool<RefineType>& refiners,
                               std::shared_ptr<SAMRAI::hier::RefineOperator>& refineOp)
        {
            for (auto const& ghostVec : ghostVecs)
            {
                refiners.addTimeRefiner(ghostVec, modelVec, oldModelVec, resourcesManager_,
                                        refineOp, fieldTimeOp_, ghostVec.vecName);
            }
        }


        /** @brief convenience overload of above with defauled field refinement operator*/
        template<typename RefinerT, RefinerT RefineType>
        void fillTimeRefiners_(std::vector<VecFieldDescriptor> const& ghostVecs,
                               VecFieldDescriptor const& modelVec,
                               VecFieldDescriptor const& oldModelVec,
                               RefinerPool<RefineType>& refiners)
        {
            fillTimeRefiners_(ghostVecs, modelVec, oldModelVec, refiners, fieldRefineOp_);
        }



        template<typename Descriptors, typename RefinerPool>
        void fillStaticRefiners_(Descriptors const& dest, Descriptors const& source,
                                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                                 RefinerPool& refiners, std::vector<std::string> keys,
                                 bool ghosts = false)
        {
            assert(dest.size() == source.size());
            auto key = std::begin(keys);
            for (std::size_t i = 0; i < dest.size(); ++i)
            {
                refiners.addStaticRefiner(dest[i], source[i], refineOp, *key++, resourcesManager_,
                                          ghosts);
            }
        }

        template<typename Descriptors, typename RefinerPool>
        void fillStaticRefiners_(Descriptors const& descriptors,
                                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                                 RefinerPool& refiners, std::vector<std::string> keys,
                                 bool ghosts = false)
        {
            fillStaticRefiners_(descriptors, descriptors, refineOp, refiners, keys, ghosts);
        }


        template<typename RefinerPool>
        void fillStaticRefiners_(VecFieldDescriptor const& descriptor,
                                 std::shared_ptr<SAMRAI::hier::RefineOperator> refineOp,
                                 RefinerPool& refiners, std::string key, bool ghosts = false)
        {
            refiners.addStaticRefiner(descriptor, descriptor, refineOp, key, resourcesManager_,
                                      ghosts);
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

                    core::empty(levelGhostParticles);
                    std::copy(std::begin(levelGhostParticlesOld), std::end(levelGhostParticlesOld),
                              std::back_inserter(levelGhostParticles));
                }
            }
        }




        double timeInterpCoef_(double const afterPushTime, std::size_t levelNumber)
        {
            return (afterPushTime - beforePushCoarseTime_[levelNumber])
                   / (afterPushCoarseTime_[levelNumber] - beforePushCoarseTime_[levelNumber]);
        }



        //! keeps a copy of the model electromagnetic field at t=n
        //        ElectromagT EM_old_{stratName + "_EM_old"}; // TODO needs to be allocated
        //        somewhere and
        // updated to t=n before advanceLevel()

        // VecFieldT Eavg_{stratName + "_avg", core::HybridQuantity::Vector::E};
        VecFieldT Jold_{stratName + "_Jold", core::HybridQuantity::Vector::J};
        VecFieldT ViOld_{stratName + "_VBulkOld", core::HybridQuantity::Vector::V};
        FieldT* NiOld_{nullptr};
        FieldUser<FieldT> NiOldUser_{stratName + "_NiOld", NiOld_,
                                     core::HybridQuantity::Scalar::rho};


        //! ResourceManager shared with other objects (like the HybridModel)
        std::shared_ptr<ResourcesManagerT> resourcesManager_;


        int const firstLevel_;
        std::unordered_map<std::size_t, double> beforePushCoarseTime_;
        std::unordered_map<std::size_t, double> afterPushCoarseTime_;

        core::Interpolator<dimension, interpOrder> interpolate_;

        //! store communicators for magnetic fields that need to be initialized
        RefinerPool<RefinerType::InitField> magneticInit_;
        //! store communicators for electric fields that need to be initializes
        RefinerPool<RefinerType::InitField> electricInit_;



        //! store communicators for magnetic fields that need ghosts to be filled
        RefinerPool<RefinerType::SharedBorder> magneticSharedNodes_;
        RefinerPool<RefinerType::GhostField> magneticGhosts_;
        RefinerPool<RefinerType::PatchGhostField> magneticPatchGhosts_;

        //! store refiners for electric fields that need ghosts to be filled
        RefinerPool<RefinerType::SharedBorder> electricSharedNodes_;
        RefinerPool<RefinerType::GhostField> electricGhosts_;

        RefinerPool<RefinerType::GhostField> currentSharedNodes_;
        RefinerPool<RefinerType::GhostField> currentGhosts_;

        // moment ghosts
        // these do not need sharedNode refiners. The reason is that
        // the border node is already complete by the deposit of ghost particles
        // these refiners are used to fill ghost nodes, and therefore, owing to
        // the GhostField tag, will only assign pur ghost nodes. Border nodes will
        // be overwritten only on level borders, which does not seem to be an issue.
        RefinerPool<RefinerType::GhostField> densityGhosts_;
        RefinerPool<RefinerType::GhostField> bulkVelGhosts_;

        // algo and schedule used to initialize domain particles
        // from coarser level using particleRefineOp<domain>
        RefinerPool<RefinerType::InitInteriorPart> interiorParticles_;

        //! store communicators for coarse to fine particles old
        RefinerPool<RefinerType::LevelBorderParticles> levelGhostParticlesOld_;

        //! store communicators for coarse to fine particles new
        RefinerPool<RefinerType::LevelBorderParticles> levelGhostParticlesNew_;

        // keys : model particles (initialization and 2nd push), temporaryParticles
        // (firstPush)
        RefinerPool<RefinerType::InteriorGhostParticles> patchGhostParticles_;

        SynchronizerPool<dimension> densitySynchronizers_;
        SynchronizerPool<dimension> ionBulkVelSynchronizers_;
        SynchronizerPool<dimension> electroSynchronizers_;
        SynchronizerPool<dimension> magnetoSynchronizers_;



        std::shared_ptr<SAMRAI::hier::RefineOperator> fieldRefineOp_{std::make_shared<
            FieldRefineOperator<GridLayoutT, FieldT, DefaultFieldRefiner<dimension>>>()};

        // see field_variable_fill_pattern.hpp for explanation about this "node_only" flag
        // Note that refinement operator, via the boolean argument, serve as a relay for the
        // the RefineAlgorithm to get the correct VariableFillPattern
        std::shared_ptr<SAMRAI::hier::RefineOperator> BfieldNodeRefineOp_{std::make_shared<
            FieldRefineOperator<GridLayoutT, FieldT, MagneticFieldRefiner<dimension>>>(
            /*node_only*/ true)};

        std::shared_ptr<SAMRAI::hier::RefineOperator> BfieldRefineOp_{std::make_shared<
            FieldRefineOperator<GridLayoutT, FieldT, MagneticFieldRefiner<dimension>>>()};

        std::shared_ptr<SAMRAI::hier::RefineOperator> EfieldNodeRefineOp_{std::make_shared<
            FieldRefineOperator<GridLayoutT, FieldT, ElectricFieldRefiner<dimension>>>(
            /*node_only*/ true)};

        std::shared_ptr<SAMRAI::hier::RefineOperator> EfieldRefineOp_{std::make_shared<
            FieldRefineOperator<GridLayoutT, FieldT, ElectricFieldRefiner<dimension>>>()};

        // field data time op
        std::shared_ptr<SAMRAI::hier::TimeInterpolateOperator> fieldTimeOp_{
            std::make_shared<FieldLinearTimeInterpolate<GridLayoutT, FieldT>>()};


        std::shared_ptr<SAMRAI::hier::RefineOperator> interiorParticleRefineOp_{
            std::make_shared<InteriorParticleRefineOp>()};

        std::shared_ptr<SAMRAI::hier::RefineOperator> levelGhostParticlesOldOp_{
            std::make_shared<CoarseToFineRefineOpOld>()};

        std::shared_ptr<SAMRAI::hier::RefineOperator> levelGhostParticlesNewOp_{
            std::make_shared<CoarseToFineRefineOpNew>()};


        std::shared_ptr<SAMRAI::hier::CoarsenOperator> fieldCoarseningOp_{std::make_shared<
            FieldCoarsenOperator<GridLayoutT, FieldT, DefaultFieldCoarsener<dimension>>>()};

        std::shared_ptr<SAMRAI::hier::CoarsenOperator> magneticCoarseningOp_{std::make_shared<
            FieldCoarsenOperator<GridLayoutT, FieldT, MagneticFieldCoarsener<dimension>>>()};
    };

    template<typename HybridModel, typename RefinementParams>
    const std::string HybridHybridMessengerStrategy<HybridModel, RefinementParams>::stratName
        = "HybridModel-HybridModel";

} // namespace amr

} // namespace PHARE

#endif
