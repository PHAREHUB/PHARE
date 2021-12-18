#ifndef PHARE_MESSENGER_MANAGER_H
#define PHARE_MESSENGER_MANAGER_H



#include "amr/messengers/hybrid_hybrid_messenger_strategy.hpp"
#include "amr/messengers/hybrid_messenger.hpp"
#include "amr/messengers/messenger.hpp"
#include "amr/messengers/mhd_hybrid_messenger_strategy.hpp"
#include "amr/messengers/mhd_messenger.hpp"

#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace PHARE::amr
{
struct MessengerDescriptor
{
    std::string coarseModel;
    std::string fineModel;
};




std::vector<MessengerDescriptor> makeDescriptors(std::vector<std::string> modelNames);




template<typename MHDModel, typename HybridModel, typename RefinementParams>
class MessengerFactory
{
    using HybridHybridMessengerStrategy_t
        = HybridHybridMessengerStrategy<HybridModel, RefinementParams>;
    using IPhysicalModel = typename HybridModel::Interface;
    static_assert(std::is_same_v<typename HybridModel::Interface, typename MHDModel::Interface>,
                  "MHD and Hybrid model need to have the same interface");

public:
    static constexpr auto dimension = HybridModel::dimension;
    static_assert(dimension == MHDModel::dimension,
                  "MHDModel::dimension != HybridModel::dimension");


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




    std::unique_ptr<IMessenger<IPhysicalModel>> create(std::string messengerName,
                                                       IPhysicalModel const& coarseModel,
                                                       IPhysicalModel const& fineModel,
                                                       int const firstLevel) const
    {
        if (messengerName == HybridHybridMessengerStrategy_t::stratName)
        {
            auto resourcesManager = dynamic_cast<HybridModel const&>(coarseModel).resourcesManager;

            auto messengerStrategy = std::make_unique<HybridHybridMessengerStrategy_t>(
                std::move(resourcesManager), firstLevel);

            return std::make_unique<HybridMessenger<HybridModel>>(std::move(messengerStrategy));
        }



        else if (messengerName == MHDHybridMessengerStrategy<MHDModel, HybridModel>::stratName)
        {
            // caution we move them so don't put a ref
            auto mhdResourcesManager = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;
            auto hybridResourcesManager
                = dynamic_cast<HybridModel const&>(fineModel).resourcesManager;

            auto messengerStrategy
                = std::make_unique<MHDHybridMessengerStrategy<MHDModel, HybridModel>>(
                    std::move(mhdResourcesManager), std::move(hybridResourcesManager), firstLevel);

            return std::make_unique<HybridMessenger<HybridModel>>(std::move(messengerStrategy));
        }




        else if (messengerName == MHDMessenger<MHDModel>::stratName)
        {
            auto mhdResourcesManager = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;

            return std::make_unique<MHDMessenger<MHDModel>>(std::move(mhdResourcesManager),
                                                            firstLevel);
        }
        else
            return {};
    }


private:
    std::vector<MessengerDescriptor> descriptors_;
};

} // namespace PHARE::amr



#endif
