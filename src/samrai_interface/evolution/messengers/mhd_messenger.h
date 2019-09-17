
#ifndef PHARE_MHD_MESSENGER_H
#define PHARE_MHD_MESSENGER_H

#include <memory>
#include <string>

#include <SAMRAI/hier/CoarsenOperator.h>
#include <SAMRAI/hier/PatchLevel.h>
#include <SAMRAI/hier/RefineOperator.h>

#include "evolution/messengers/messenger.h"
#include "evolution/messengers/messenger_info.h"
#include "hybrid/hybrid_quantities.h"
#include "physical_models/mhd_model.h"
#include "physical_models/physical_model.h"

namespace PHARE
{
namespace amr_interface
{
    template<typename MHDModel>
    class MHDMessenger : public IMessenger
    {
    public:
        MHDMessenger(std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager,
                     int const firstLevel)
            : resourcesManager_{std::move(resourcesManager)}
            , firstLevel_{firstLevel}
        {
        }

        virtual void
        registerQuantities(std::unique_ptr<IMessengerInfo> fromCoarserInfo,
                           [[maybe_unused]] std::unique_ptr<IMessengerInfo> fromFinerInfo) override
        {
            std::unique_ptr<MHDMessengerInfo> mhdInfo{
                dynamic_cast<MHDMessengerInfo*>(fromCoarserInfo.release())};
        }



        virtual void registerLevel(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                                   int const levelNumber) override
        {
        }


        static const std::string stratName;

        virtual std::string fineModelName() const override { return MHDModel::model_name; }

        virtual std::string coarseModelName() const override { return MHDModel::model_name; }

        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) const override
        {
        }

        virtual void initLevel(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               double const initDataTime) override
        {
        }

        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromCoarser() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }

        virtual std::unique_ptr<IMessengerInfo> emptyInfoFromFiner() override
        {
            return std::make_unique<MHDMessengerInfo>();
        }




        virtual void regrid(std::shared_ptr<SAMRAI::hier::PatchHierarchy> const& hierarchy,
                            const int levelNumber,
                            std::shared_ptr<SAMRAI::hier::PatchLevel> const& oldLevel,
                            IPhysicalModel& model, double const initDataTime) override
        {
        }


        virtual void firstStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                               const std::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                               double time) final
        {
        }


        virtual void lastStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) final {}


        virtual void prepareStep(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level) final {}

        virtual void fillRootGhosts(IPhysicalModel& model, SAMRAI::hier::PatchLevel& level,
                                    double const initDataTime) final
        {
        }



        virtual void synchronize(SAMRAI::hier::PatchLevel& level) final
        {
            (void)level;
            // call coarsning schedules...
        }



        virtual std::string name() override { return stratName; }

        virtual ~MHDMessenger() = default;


    private:
        std::shared_ptr<typename MHDModel::resources_manager_type> resourcesManager_;
        int const firstLevel_;
    };


    template<typename MHDModel>
    const std::string MHDMessenger<MHDModel>::stratName = "MHDModel-MHDModel";
} // namespace amr_interface
} // namespace PHARE
#endif
