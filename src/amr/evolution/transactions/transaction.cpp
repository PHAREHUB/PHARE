

#include "transaction.h"


namespace PHARE
{
bool areCompatible(ITransaction const& transaction, PhysicalModel const& model)
{
    return model.name() == transaction.fineModelName()
           || model.name() == transaction.coarseModelName();
}
} // namespace PHARE
