

#include "transaction.h"


namespace PHARE
{
bool areCompatible(ITransaction const& transaction, IPhysicalModel const& model)
{
    return model.name() == transaction.fineModelName()
           || model.name() == transaction.coarseModelName();
}
} // namespace PHARE
