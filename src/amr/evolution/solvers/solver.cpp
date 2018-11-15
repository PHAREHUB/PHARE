

#include "evolution/solvers/solver.h"

namespace PHARE
{
bool areCompatible(ITransaction const& transaction, ISolver const& solver)
{
    return solver.modelName() == transaction.fineModelName()
           || solver.modelName() == transaction.coarseModelName();
}


bool areCompatible(IPhysicalModel const& model, ISolver const& solver)
{
    return model.name() == solver.modelName();
}


} // namespace PHARE
