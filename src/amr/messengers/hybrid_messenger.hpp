#ifndef PHARE_HYBRID_MESSENGER_HPP
#define PHARE_HYBRID_MESSENGER_HPP


#include "core/def.hpp"
#include <core/logger.hpp>

#include "amr/messengers/messenger.hpp"
#include "amr/messengers/messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"



namespace PHARE
{
namespace amr
{
    enum class HybridMessengerStrategyType { HybridHybrid, MHDHybrid };



    /** @brief HybridMessenger is a concrete implementation of the IMessenger interface which
     * represents all data messengers towards a Hybrid level.
     *
     * Most of the operations needed to communicate data cannot however be defined here since it is
     * not known whether the origin level is Hybrid or not. The class therefore relies on a pointer
     * to a strategy that performs the operations.
     *
     * Possible HybridMessengerStrategy are:
     *
     * - HybridHybridMessengerStrategy
     * - MHDHybridMessengerStrategy
     *
     *
     * In addition to implement the interface of a IMessenger to be used by the
     * MultiPhysicsIntegrator, HybridMessenger also provides Hybrid ISolver objects with an
     * interface specific to Hybrid systems. Those methods are:
     *
     *
     * Ghost filling methods tuned for Hybrid quantities:
     *
     * - fillElectricGhosts()
     * - fillIonGhostParticles()
     * - fillIonPopMomentGhosts()
     *
     * Synchronization methods of Hybrid quantities:
     *
     * - syncMagnetic()
     * - syncElectric()
     * - syncIonMOments()
     *
     */
    template<typename HybridModel>
    class HybridMessenger : public IMessenger<typename HybridModel::Interface>
    {
    private:
        using IonsT          = decltype(std::declval<HybridModel>().state.ions);
        using VecFieldT      = decltype(std::declval<HybridModel>().state.electromag.E);
        using FieldT         = VecFieldT::field_type;
        using IPhysicalModel = HybridModel::Interface;
        using amr_types      = HybridModel::amr_types;
        using patch_t        = amr_types::patch_t;
        using level_t        = amr_types::level_t;
        using hierarchy_t    = amr_types::hierarchy_t;
        using level_ptr      = amr_types::level_ptr;
        using hierarchy_ptr  = amr_types::hierarchy_ptr;


        using stratT = HybridMessengerStrategy<HybridModel>;

    public:
        explicit HybridMessenger(std::unique_ptr<stratT> strat)
            : strat_{std::move(strat)}
        {
        }


        /* -------------------------------------------------------------------------
                                Start IMessenger interface
           -------------------------------------------------------------------------*/



        /**
         * @brief see IMessenger::allocate. Allocate calls the abstract HybridMessengerStrategy to
         * perform the allocation of its internal resources
         */
        void allocate(patch_t& patch, double const allocateTime) const override
        {
            strat_->allocate(patch, allocateTime);
        }



        /**
         * @brief see IMessenger::registerQuantities. Prepares a HybridMessenger to communicates
         * variables given in the two information structures. Since the HybridMessenger does not
         * know which concrete strategy it is using, it cannot directly use these informations. It
         * thus passes them to its strategy.
         *
         * @param fromCoarserInfo see IMessenger
         * @param fromFinerInfo see IMessenger
         */
        void registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                                std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            strat_->registerQuantities(std::move(fromCoarserInfo), std::move(fromFinerInfo));
        }




        /**
         * @brief see IMessenger::registerLevel
         */
        void registerLevel(hierarchy_ptr const& hierarchy, int const levelNumber) override
        {
            strat_->registerLevel(hierarchy, levelNumber);
        }

        /**
         * @brief see IMessenger::onHierarchyChanged
         */
        void onHierarchyChanged(hierarchy_ptr const& hierarchy, int const levelNumber) override
        {
            strat_->onHierarchyChanged(hierarchy, levelNumber);
        }



        /**
         * @brief see IMessenger::registerLevel
         */
        void regrid(hierarchy_ptr const& hierarchy, int const levelNumber,
                    level_ptr const& oldLevel, IPhysicalModel& model,
                    double const initDataTime) final
        {
            strat_->regrid(hierarchy, levelNumber, oldLevel, model, initDataTime);
        }




        /**
         * @brief see IMessenger::initLevel
         * @param levelNumber
         * @param initDataTime
         */
        void initLevel(hierarchy_ptr const& hierarchy, IPhysicalModel& model, level_t& level,
                       double const initDataTime) override
        {
            strat_->initLevel(hierarchy, model, level, initDataTime);
        }


        /**
         * @brief see IMessenger::firstStep
         * @param model
         */
        void firstStep(IPhysicalModel& model, level_t& level, double const currentTime,
                       double const prevCoarserTime, double const newCoarserTime) final
        {
            strat_->firstStep(model, level, currentTime, prevCoarserTime, newCoarserTime);
        }



        /**
         * @brief see IMessenger::lastStep
         * @param model
         */
        void lastStep(IPhysicalModel& model, level_t& level) override
        {
            strat_->lastStep(model, level);
        }



        /**
         * @brief prepareStep see IMessenger::prepareStep
         */
        void prepareStep(IPhysicalModel& model, level_t& level, double currentTime) override
        {
            strat_->prepareStep(model, level, currentTime);
        }



        void fillRootGhosts(IPhysicalModel& model, level_t& level,
                            double const initDataTime) override
        {
            strat_->fillRootGhosts(model, level, initDataTime);
        }



        void synchronize(level_t& level) override { strat_->synchronize(level); }


        void reflux(int const coarserLevelNumber, int const fineLevelNumber,
                    double const syncTime) override
        {
            strat_->reflux(coarserLevelNumber, fineLevelNumber, syncTime);
        }

        // fully kinetic: coarsens the finer level's extrapolated covered E onto the coarser
        // level's copy of its E, for the coarse B at midTime (SolverFK::partialReflux)
        void partialReflux(int const coarserLevelNumber, int const fineLevelNumber,
                           double const midTime)
        {
            strat_->partialReflux(coarserLevelNumber, fineLevelNumber, midTime);
        }

        // fully kinetic: the E of a level (re)created during the run, at its own half step
        // (FKLevelInitializer, P5a)
        void regridElectric(VecFieldT& E, hierarchy_ptr const& hierarchy, int const levelNumber,
                            level_ptr const& oldLevel, double const fillTime)
        {
            strat_->regridElectric(E, hierarchy, levelNumber, oldLevel, fillTime);
        }

        // fully kinetic: the E ghosts other patches of the level own, before a cleaning pass
        void fillElectricPatchGhosts(VecFieldT& E, level_t const& level, double const fillTime)
        {
            strat_->fillElectricPatchGhosts(E, level, fillTime);
        }


        void postSynchronize(IPhysicalModel& model, level_t& level, double const time) override
        {
            strat_->postSynchronize(model, level, time);
        }



        /**
         * @brief IMessenger::fineModelName
         * @return
         */
        NO_DISCARD std::string fineModelName() const override { return strat_->fineModelName(); }




        /**
         * @brief see IMessenger::coarseModelName
         * @return
         */
        NO_DISCARD std::string coarseModelName() const override
        {
            return strat_->coarseModelName();
        }




        /**
         * @brief see IMessenger::emptyInfoFromCoarser
         */
        NO_DISCARD std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return strat_->emptyInfoFromCoarser();
        }



        /**
         * @brief see IMessenger::emptyInfoFromFiner
         */
        NO_DISCARD std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return strat_->emptyInfoFromFiner();
        }



        /**
         * @brief returns the name of the concrete IMessenger, which in the case of a
         * HybridMessenger is just the name of its strategy.
         */
        NO_DISCARD std::string name() override
        {
            if (strat_ != nullptr)
            {
                return strat_->name();
            }
            else
            {
                throw std::runtime_error("invalid strat");
            }
        }

        /* -------------------------------------------------------------------------
                                End IMessenger interface

                            Start HybridMessenger Interface
           -------------------------------------------------------------------------*/

        void fillMagneticGhosts(VecFieldT& B, level_t const& level, double const fillTime)
        {
            strat_->fillMagneticGhosts(B, level, fillTime);
        }


        /**
         * @brief fillElectricGhosts is called by a ISolver solving a hybrid equatons to fill
         * the ghost nodes of the electric field
         * @param E is the electric field for which ghost nodes will be filled
         * @param levelNumber
         * @param fillTime
         */
        void fillElectricGhosts(VecFieldT& E, level_t const& level, double const fillTime)
        {
            strat_->fillElectricGhosts(E, level, fillTime);
        }



        /**
         * @brief fillIonGhostParticles is called by a ISolver solving hybrid equations to fill the
         * ghosts particles
         * @param ions for which ghost particles will be filled
         * @param levelNumber
         * @param fillTime
         */
        void fillIonGhostParticles(IonsT& ions, level_t& level, double const fillTime)
        {
            strat_->fillIonGhostParticles(ions, level, fillTime);
        }

        void fillBandLeavers(IonsT& ions, level_t& level, double const fillTime)
        {
            strat_->fillBandLeavers(ions, level, fillTime);
        }



        /**
         * @brief fillIonPopMomentGhosts is called by a ISolver solving hybrid equations to fill the
         * ion moments
         * @param ions
         * @param levelNumber
         * @param fillTime
         */
        void fillIonPopMomentGhosts(IonsT& ions, level_t& level, double const fillTime)
        {
            strat_->fillIonPopMomentGhosts(ions, level, fillTime);
        }


        // synchronization/coarsening methods

        void syncMagnetic(VecFieldT& B) { strat_->syncMagnetic(B); }
        void syncElectric(VecFieldT& E) { strat_->syncElectric(E); }
        void syncIonMoments(IonsT& ions) { strat_->syncIonMoments(ions); }


        void fillFluxBorders(IonsT& ions, level_t& level, double const fillTime)
        {
            strat_->fillFluxBorders(ions, level, fillTime);
        }

        void fillDensityBorders(IonsT& ions, level_t& level, double const fillTime)
        {
            strat_->fillDensityBorders(ions, level, fillTime);
        }

        void fillIonBorders(IonsT& ions, level_t& level, double const fillTime)
        {
            strat_->fillIonBorders(ions, level, fillTime);
        }

        void fillCurrentBorders(VecFieldT& J, level_t& level, double const fillTime)
        {
            strat_->fillCurrentBorders(J, level, fillTime);
        }

        void equalizeCurrentBorders(VecFieldT& J, level_t& level, double const fillTime)
        {
            strat_->equalizeCurrentBorders(J, level, fillTime);
        }

        void fillCurrentGhosts(VecFieldT& J, level_t& level, double const fillTime)
        {
            strat_->fillCurrentGhosts(J, level, fillTime);
        }

        void ownSharedElectromag(VecFieldT& E, VecFieldT& B, level_t& level)
        {
            strat_->ownSharedElectromag(E, B, level);
        }

        void checkSharedElectromagCopies(VecFieldT& E, VecFieldT& B, level_t& level)
        {
            strat_->checkSharedElectromagCopies(E, B, level);
        }

        void fillChargeDensityBorders(FieldT& rho, level_t& level, double const fillTime)
        {
            strat_->fillChargeDensityBorders(rho, level, fillTime);
        }

        void equalizeChargeDensityBorders(FieldT& rho, level_t& level, double const fillTime)
        {
            strat_->equalizeChargeDensityBorders(rho, level, fillTime);
        }

        void fillChargeDensityGhosts(FieldT& rho, level_t& level, double const fillTime)
        {
            strat_->fillChargeDensityGhosts(rho, level, fillTime);
        }


        /* -------------------------------------------------------------------------
                            End HybridMessenger Interface
           -------------------------------------------------------------------------*/



        virtual ~HybridMessenger() = default;


    private:
        std::unique_ptr<stratT> const strat_;
    };


} // namespace amr

} // namespace PHARE
#endif
