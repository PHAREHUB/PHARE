#ifndef PHARE_MESSENGER_MANAGER_H
#define PHARE_MESSENGER_MANAGER_H


#include <algorithm>
#include <memory>
#include <optional>
#include <string>

#include "evolution/messengers/hybrid_hybrid_messenger_strategy.h"
#include "evolution/messengers/hybrid_messenger.h"
#include "evolution/messengers/messenger.h"
#include "evolution/messengers/mhd_hybrid_messenger_strategy.h"
#include "evolution/messengers/mhd_messenger.h"
#include "physical_models/physical_model.h"


namespace PHARE
{
namespace amr_interface
{
    struct MessengerDescriptor
    {
        std::string coarseModel;
        std::string fineModel;
    };


    template<typename MHDModel, typename HybridModel>
    class MessengerFactory
    {
    public:
        MessengerFactory(std::vector<MessengerDescriptor> messengerDescriptors)
            : descriptors_{messengerDescriptors}
        {
        }




        std::optional<std::string> name(IPhysicalModel const& coarseModel,
                                        IPhysicalModel const& fineModel) const
        {
            auto finder = [&coarseModel, &fineModel](MessengerDescriptor const& desc) {
                return desc.coarseModel == coarseModel.name() && desc.fineModel == fineModel.name();
            };

            auto messenger = std::find_if(std::begin(descriptors_), std::end(descriptors_), finder);

            if (messenger != std::end(descriptors_))
            {
                return coarseModel.name() + "-" + fineModel.name();
            }
            else
            {
                return {};
            }
        }




        std::unique_ptr<IMessenger> create(std::string messengerName,
                                           IPhysicalModel const& coarseModel,
                                           IPhysicalModel const& fineModel,
                                           int const firstLevel) const
        {
            if (messengerName == HybridHybridMessengerStrategy<HybridModel>::stratName)
            {
                auto resourcesManager
                    = dynamic_cast<HybridModel const&>(coarseModel).resourcesManager;

                auto messengerStrategy
                    = std::make_unique<HybridHybridMessengerStrategy<HybridModel>>(
                        std::move(resourcesManager), firstLevel);

                return std::make_unique<HybridMessenger<HybridModel>>(std::move(messengerStrategy));
            }



            else if (messengerName == MHDHybridMessengerStrategy<MHDModel, HybridModel>::stratName)
            {
                // caution we move them so don't put a ref
                auto mhdResourcesManager
                    = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;
                auto hybridResourcesManager
                    = dynamic_cast<HybridModel const&>(fineModel).resourcesManager;

                auto messengerStrategy
                    = std::make_unique<MHDHybridMessengerStrategy<MHDModel, HybridModel>>(
                        std::move(mhdResourcesManager), std::move(hybridResourcesManager),
                        firstLevel);

                return std::make_unique<HybridMessenger<HybridModel>>(std::move(messengerStrategy));
            }




            else if (messengerName == MHDMessenger<MHDModel>::stratName)
            {
                auto mhdResourcesManager
                    = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;

                return std::make_unique<MHDMessenger<MHDModel>>(std::move(mhdResourcesManager),
                                                                firstLevel);
            }
            else
                return {};
        }


    private:
        std::vector<MessengerDescriptor> descriptors_;
    };

} // namespace amr_interface
} // namespace PHARE


#endif
