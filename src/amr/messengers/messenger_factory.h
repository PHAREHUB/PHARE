#ifndef PHARE_MESSENGER_MANAGER_H
#define PHARE_MESSENGER_MANAGER_H


#include <algorithm>
#include <memory>
#include <optional>
#include <string>

#include "amr/messengers/hybrid_hybrid_messenger_strategy.h"
#include "amr/messengers/hybrid_messenger.h"
#include "amr/messengers/messenger.h"
#include "amr/messengers/mhd_hybrid_messenger_strategy.h"
#include "amr/messengers/mhd_messenger.h"


namespace PHARE
{
namespace amr
{
    struct MessengerDescriptor
    {
        std::string coarseModel;
        std::string fineModel;
    };




    std::vector<MessengerDescriptor> makeDescriptors(std::vector<std::string> modelNames)
    {
        if (modelNames.size() == 1)
            return {{modelNames[0], modelNames[0]}};

        else if (modelNames.size() == 2)
        {
            return {{modelNames[0], modelNames[0]},
                    {modelNames[0], modelNames[1]},
                    {modelNames[1], modelNames[1]}};
        }
        else
            throw std::runtime_error("Error max number of models is 2");
    }


    template<typename MHDModel, typename HybridModel, typename IPhysicalModel>
    class MessengerFactory
    {
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
            if (messengerName
                == HybridHybridMessengerStrategy<HybridModel, IPhysicalModel>::stratName)
            {
                auto resourcesManager
                    = dynamic_cast<HybridModel const&>(coarseModel).resourcesManager;

                auto messengerStrategy
                    = std::make_unique<HybridHybridMessengerStrategy<HybridModel, IPhysicalModel>>(
                        std::move(resourcesManager), firstLevel);

                return std::make_unique<HybridMessenger<HybridModel, IPhysicalModel>>(
                    std::move(messengerStrategy));
            }



            else if (messengerName
                     == MHDHybridMessengerStrategy<MHDModel, HybridModel,
                                                   IPhysicalModel>::stratName)
            {
                // caution we move them so don't put a ref
                auto mhdResourcesManager
                    = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;
                auto hybridResourcesManager
                    = dynamic_cast<HybridModel const&>(fineModel).resourcesManager;

                auto messengerStrategy = std::make_unique<
                    MHDHybridMessengerStrategy<MHDModel, HybridModel, IPhysicalModel>>(
                    std::move(mhdResourcesManager), std::move(hybridResourcesManager), firstLevel);

                return std::make_unique<HybridMessenger<HybridModel, IPhysicalModel>>(
                    std::move(messengerStrategy));
            }




            else if (messengerName == MHDMessenger<MHDModel, IPhysicalModel>::stratName)
            {
                auto mhdResourcesManager
                    = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;

                return std::make_unique<MHDMessenger<MHDModel, IPhysicalModel>>(
                    std::move(mhdResourcesManager), firstLevel);
            }
            else
                return {};
        }


    private:
        std::vector<MessengerDescriptor> descriptors_;
    };

} // namespace amr
} // namespace PHARE


#endif
