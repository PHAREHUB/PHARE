#ifndef PHARE_MHD_HYBRID_MESSENGER_STRATEGY_HPP
#define PHARE_MHD_HYBRID_MESSENGER_STRATEGY_HPP

#include "amr/messengers/hybrid_messenger_info.hpp"
#include "amr/messengers/hybrid_messenger_strategy.hpp"
#include "amr/messengers/mhd_messenger_info.hpp"

#include <string>

namespace PHARE
{
namespace amr
{
    template<typename MHDModel, typename HybridModel>
    class MHDHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel>
    {
        using IonsT          = decltype(std::declval<HybridModel>().state.ions);
        using VecFieldT      = decltype(std::declval<HybridModel>().state.electromag.E);
        using IPhysicalModel = typename HybridModel::Interface;

        // we will probably need both resources managers if we still have 2 in the future
        using resources_manager_type = HybridModel::resources_manager_type;

    public:
        static inline std::string const stratName = "MHDModel-HybridModel";

        MHDHybridMessengerStrategy(std::shared_ptr<resources_manager_type> const& resourcesManager,
                                   int const firstLevel)
            : HybridMessengerStrategy<HybridModel>{stratName}
            , resourcesManager_{resourcesManager}
            , firstLevel_{firstLevel}
        {
            resourcesManager_->registerResources(EM_old_);
        }

        /**
         * @brief allocate the internal resources to the hybrid and MHD
         * resourcesManagers
         */
        void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            // hybModel.resourcesManager->allocate(EM_old_.E, patch, allocateTime);
            // hybModel.resourcesManager->allocate(EM_old_.B, patch, allocateTime);
            resourcesManager_->allocate(EM_old_, patch, allocateTime);
        }

        void registerQuantities(
            std::unique_ptr<IMessengerInfo> /*fromCoarserInfo*/,
            [[maybe_unused]] std::unique_ptr<IMessengerInfo> /*fromFinerInfo*/) override
        {
        }

        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                           int const /*levelNumber*/) override
        {
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<HybridMessengerInfo>();
        }

        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                    int const /*levelNumber*/,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*oldLevel*/,
                    IPhysicalModel& /*model*/, double const /*initDataTime*/) override
        {
            //
        }

        std::string fineModelName() const override { return HybridModel::model_name; }

        std::string coarseModelName() const override { return MHDModel::model_name; }

        void initLevel(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       double const /*initDataTime*/) override
        {
        }

        virtual ~MHDHybridMessengerStrategy() = default;

        void fillMagneticGhosts(VecFieldT& /*B*/, SAMRAI::hier::PatchLevel const& /*level*/,
                                double const /*fillTime*/) override
        {
        }

        void fillElectricGhosts(VecFieldT& /*E*/, SAMRAI::hier::PatchLevel const& /*level*/,
                                double const /*fillTime*/) override
        {
        }

        void fillCurrentGhosts(VecFieldT& /*J*/, SAMRAI::hier::PatchLevel const& /*level*/,
                               double const /*fillTime*/) override
        {
        }

        void fillIonGhostParticles(IonsT& /*ions*/, SAMRAI::hier::PatchLevel& /*level*/,
                                   double const /*fillTime*/) override
        {
        }
        void fillIonPopMomentGhosts(IonsT& /*ions*/, SAMRAI::hier::PatchLevel& /*level*/,
                                    double const /*fillTime*/) override
        {
        }

        // void fillIonMomentGhosts(IonsT& /*ions*/, SAMRAI::hier::PatchLevel& /*level*/,
        //                          double const /*fillTime*/) override
        // {
        // }


        void fillFluxBorders(IonsT& /*ions*/, SAMRAI::hier::PatchLevel& /*level*/,
                             double const /*fillTime*/) override
        {
        }
        void fillDensityBorders(IonsT& /*ions*/, SAMRAI::hier::PatchLevel& /*level*/,
                                double const /*fillTime*/) override
        {
        }
        void fillIonBorders(IonsT& /*ions*/, SAMRAI::hier::PatchLevel& /*level*/,
                            double const /*fillTime*/) override
        {
        }



        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                       double const /*currentTime*/, double const /*prevCoarserTime*/,
                       double const /*newCoarserTime*/) override
        {
        }

        void lastStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/) override {}

        void prepareStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                         double /*currentTime*/) final
        {
        }

        void fillRootGhosts(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                            double const /*initDataTime*/) final
        {
        }

        void synchronize(SAMRAI::hier::PatchLevel& /*level*/) final
        {
            // call coarsning schedules...
        }

        void reflux(int const /*coarserLevelNumber*/, int const /*fineLevelNumber*/,
                    double const /*syncTime*/) override
        {
        }

        void postSynchronize(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                             double const /*time*/) override
        {
        }

    private:
        using Electromag = decltype(std::declval<HybridModel>().state.electromag);

        std::shared_ptr<resources_manager_type> resourcesManager_;
        int const firstLevel_;
        Electromag EM_old_{stratName + "_EM_old"};
    };


} // namespace amr
} // namespace PHARE

#endif
