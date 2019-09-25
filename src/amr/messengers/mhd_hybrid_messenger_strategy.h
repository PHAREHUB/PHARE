
#ifndef PHARE_MHD_HYBRID_MESSENGER_STRATEGY_H
#define PHARE_MHD_HYBRID_MESSENGER_STRATEGY_H

#include "messengers/hybrid_messenger_info.h"
#include "messengers/hybrid_messenger_strategy.h"
#include "messengers/mhd_messenger_info.h"

#include <string>

namespace PHARE
{
namespace amr
{
    template<typename MHDModel, typename HybridModel, typename IPhysicalModel>
    class MHDHybridMessengerStrategy : public HybridMessengerStrategy<HybridModel, IPhysicalModel>
    {
        using IonsT     = decltype(std::declval<HybridModel>().state.ions);
        using VecFieldT = decltype(std::declval<HybridModel>().state.electromag.E);

    public:
        static const std::string stratName;

        MHDHybridMessengerStrategy(
            std::shared_ptr<typename MHDModel::resources_manager_type> mhdResourcesManager,
            std::shared_ptr<typename HybridModel::resources_manager_type> hybridResourcesManager,
            int const firstLevel)
            : HybridMessengerStrategy<HybridModel, IPhysicalModel>{stratName}
            , mhdResourcesManager_{std::move(mhdResourcesManager)}
            , hybridResourcesManager_{std::move(hybridResourcesManager)}
            , firstLevel_{firstLevel}
        {
            hybridResourcesManager_->registerResources(EM_old_);
        }

        /**
         * @brief allocate allocate the internal resources to the hybrid and MHD resourcesManagers
         */
        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
            // hybModel.resourcesManager->allocate(EM_old_.E, patch, allocateTime);
            // hybModel.resourcesManager->allocate(EM_old_.B, patch, allocateTime);
            hybridResourcesManager_->allocate(EM_old_, patch, allocateTime);
        }



        virtual void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
        }


        virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   int const levelNumber) override
        {
        }

        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }

        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<HybridMessengerInfo>();
        }




        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            const int levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            IPhysicalModel& model, double const initDataTime) override
        {
            //
        }




        virtual std::string fineModelName() const override { return HybridModel::model_name; }

        virtual std::string coarseModelName() const override { return MHDModel::model_name; }



        virtual void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               double const initDataTime) override
        {
        }

        virtual ~MHDHybridMessengerStrategy() = default;


        virtual void fillMagneticGhosts(VecFieldT& B, int const levelNumber,
                                        double const fillTime) override
        {
        }
        virtual void fillElectricGhosts(VecFieldT& E, int const levelNumber,
                                        double const fillTime) override
        {
        }
        virtual void fillIonGhostParticles(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                           double const fillTime) override
        {
        }
        virtual void fillIonMomentGhosts(IonsT& ions, SAMRAI::hier::PatchLevel& level,
                                         double const currentTime, double const fillTime) override
        {
        }



        virtual void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                               double time) override
        {
        }

        virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) override {}


        virtual void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) final {}

        virtual void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                    double const initDataTime) final
        {
        }



        virtual void synchronize(SAMRAI::hier::PatchLevel& level) final
        {
            // call coarsning schedules...
        }




    private:
        using Electromag = decltype(std::declval<HybridModel>().state.electromag);

        std::shared_ptr<typename MHDModel::resources_manager_type> mhdResourcesManager_;
        std::shared_ptr<typename HybridModel::resources_manager_type> hybridResourcesManager_;
        int const firstLevel_;
        Electromag EM_old_{stratName + "_EM_old"};
    };

    template<typename MHDModel, typename HybridModel, typename IPhysicalModel>
    const std::string MHDHybridMessengerStrategy<MHDModel, HybridModel, IPhysicalModel>::stratName
        = "MHDModel-HybridModel";


} // namespace amr
} // namespace PHARE


#endif
