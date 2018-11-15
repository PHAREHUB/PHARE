#ifndef PHARE_TRANSACTION_MANAGER_H
#define PHARE_TRANSACTION_MANAGER_H


#include <algorithm>
#include <memory>
#include <optional>
#include <string>

#include "evolution/transactions/hybrid_hybrid_transaction_strategy.h"
#include "evolution/transactions/hybrid_transaction.h"
#include "evolution/transactions/mhd_hybrid_transaction_strategy.h"
#include "evolution/transactions/mhd_transaction.h"
#include "evolution/transactions/transaction.h"
#include "physical_models/physical_model.h"


namespace PHARE
{
struct TransactionDescriptor
{
    std::string coarseModel;
    std::string fineModel;
};


template<typename MHDModel, typename HybridModel>
class TransactionFactory
{
public:
    TransactionFactory(std::vector<TransactionDescriptor> transactionDescriptors)
        : descriptors_{transactionDescriptors}
    {
    }

    std::optional<std::string> name(IPhysicalModel const& coarseModel,
                                    IPhysicalModel const& fineModel) const
    {
        auto finder = [&coarseModel, &fineModel](TransactionDescriptor const& desc) {
            return desc.coarseModel == coarseModel.name() && desc.fineModel == fineModel.name();
        };

        auto transaction = std::find_if(std::begin(descriptors_), std::end(descriptors_), finder);

        if (transaction != std::end(descriptors_))
        {
            return coarseModel.name() + "-" + fineModel.name();
        }
        else
        {
            return {};
        }
    }

    std::unique_ptr<ITransaction> create(std::string transactionName,
                                         IPhysicalModel const& coarseModel,
                                         IPhysicalModel const& fineModel, int const firstLevel) const
    {
        if (transactionName == HybridHybridTransactionStrategy<HybridModel>::stratName)
        {
            auto resourcesManager = dynamic_cast<HybridModel const&>(coarseModel).resourcesManager;

            auto transactionStrategy
                = std::make_unique<HybridHybridTransactionStrategy<HybridModel>>(
                    std::move(resourcesManager), firstLevel);

            return std::make_unique<HybridTransaction<HybridModel>>(std::move(transactionStrategy));
        }




        else if (transactionName == MHDHybridTransactionStrategy<MHDModel, HybridModel>::stratName)
        {
            // caution we move them so don't put a ref
            auto mhdResourcesManager = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;
            auto hybridResourcesManager
                = dynamic_cast<HybridModel const&>(fineModel).resourcesManager;

            auto transactionStrategy
                = std::make_unique<MHDHybridTransactionStrategy<MHDModel, HybridModel>>(
                    std::move(mhdResourcesManager), std::move(hybridResourcesManager), firstLevel);

            return std::make_unique<HybridTransaction<HybridModel>>(std::move(transactionStrategy));
        }

        else if (transactionName == MHDTransaction<MHDModel>::stratName)
        {
            auto mhdResourcesManager = dynamic_cast<MHDModel const&>(coarseModel).resourcesManager;

            return std::make_unique<MHDTransaction<MHDModel>>(std::move(mhdResourcesManager),
                                                              firstLevel);
        }
        else
            return {};
    }


private:
    std::vector<TransactionDescriptor> descriptors_;
};

} // namespace PHARE


#endif
