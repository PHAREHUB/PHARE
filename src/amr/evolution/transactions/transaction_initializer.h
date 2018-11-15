
#ifndef PHARE_TRANSACTION_INITIALIZER_H
#define PHARE_TRANSACTION_INITIALIZER_H


#include "evolution/solvers/solver.h"
#include "evolution/transactions/transaction.h"
#include "physical_models/physical_model.h"


namespace PHARE
{
class TransactionInitializer
{
public:
    static void setup(ITransaction& transaction, IPhysicalModel const& coarseModel,
                      IPhysicalModel const& fineModel, ISolver const& solver)
    {
        auto fromCoarserInfo = transaction.emptyInfoFromCoarser();
        auto fromFinerInfo   = transaction.emptyInfoFromFiner();

        fineModel.fillTransactionInfo(fromFinerInfo);
        coarseModel.fillTransactionInfo(fromCoarserInfo);

        transaction.registerQuantities(std::move(fromCoarserInfo), std::move(fromFinerInfo));


        /*
        solver.fillTransactionInfo(fromFinerInfo);
            */

        /*
                for (auto const& user : users)
                {
                    if (areCompatible(transaction, *user))
                    {
                        transaction.setup(std::move(user->transactionInfo()));
                    }
                }
                */
    }
};



} // namespace PHARE
#endif
