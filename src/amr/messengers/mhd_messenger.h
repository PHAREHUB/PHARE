
#ifndef PHARE_MHD_MESSENGER_H
#define PHARE_MHD_MESSENGER_H

#include <memory>
#include <string>

#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/RefineOperator.h>

#include "core/hybrid/hybrid_quantities.h"
#include "amr/messengers/messenger.h"
#include "amr/messengers/messenger_info.h"
#include "amr/messengers/mhd_messenger_info.h"

namespace PHARE
{
namespace amr
{
    template<typename MHDModel>
    class MHDMessenger : public IMessenger<typename MHDModel::Interface>
    {
    public:
        using IPhysicalModel = typename MHDModel::Interface;
        MHDMessenger(std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager,
                     int const firstLevel)
            : resourcesManager_{std::move(resourcesManager)}
            , firstLevel_{firstLevel}
        {
        }

        void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<MHDMessengerInfo> mhdInfo{
                dynamic_cast<MHDMessengerInfo*>(fromCoarserInfo.release())};
        }



        void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                           int const /*levelNumber*/) override
        {
        }


        static const std::string stratName;

        std::string fineModelName() const override { return MHDModel::model_name; }

        std::string coarseModelName() const override { return MHDModel::model_name; }

        void allocate(SAMRAI::hier::Patch& /*patch*/, double const /*allocateTime*/) const override
        {
        }

        void initLevel(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       double const /*initDataTime*/) override
        {
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }

        std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }




        void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                    const int /*levelNumber*/,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*oldLevel*/,
                    IPhysicalModel& /*model*/, double const /*initDataTime*/) override
        {
        }


        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& /*hierarchy*/,
                       double /*time*/) final
        {
        }


        void lastStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/) final {}


        void prepareStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/) final {}

        void fillRootGhosts(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                            double const /*initDataTime*/) final
        {
        }



        void synchronize(SAMRAI::hier::PatchLevel& /*level*/) final
        {
            // call coarsning schedules...
        }



        std::string name() override { return stratName; }

        virtual ~MHDMessenger() = default;


    private:
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        int const firstLevel_;
    };


    template<typename MHDModel>
    const std::string MHDMessenger<MHDModel>::stratName = "MHDModel-MHDModel";
} // namespace amr
} // namespace PHARE
#endif
