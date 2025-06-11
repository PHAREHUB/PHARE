
#ifndef PHARE_MHD_MESSENGER_HPP
#define PHARE_MHD_MESSENGER_HPP

#include <memory>
#include <string>
#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/RefineOperator.h>

#include "core/hybrid/hybrid_quantities.hpp"
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


        static std::string const stratName;

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
                    int const /*levelNumber*/,
                    std::shared_ptr<SAMRAI::hier::PatchLevel> const& /*oldLevel*/,
                    IPhysicalModel& /*model*/, double const /*initDataTime*/) override
        {
        }


        void firstStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/,
                       std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& /*hierarchy*/,
                       double const /*currentTime*/, double const /*prevCoarserTIme*/,
                       double const /*newCoarserTime*/) final
        {
        }


        void lastStep(IPhysicalModel& /*model*/, SAMRAI::hier::PatchLevel& /*level*/) final {}


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


        std::string name() override { return stratName; }

        virtual ~MHDMessenger() = default;


    private:
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        int const firstLevel_;
    };


    template<typename MHDModel>
    std::string const MHDMessenger<MHDModel>::stratName = "MHDModel-MHDModel";
} // namespace amr
} // namespace PHARE
#endif
