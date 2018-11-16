

#include "evolution/solvers/solver.h"

namespace PHARE
{
bool areCompatible(IMessenger const& messenger, ISolver const& solver)
{
    return solver.modelName() == messenger.fineModelName()
           || solver.modelName() == messenger.coarseModelName();
}


bool areCompatible(IPhysicalModel const& model, ISolver const& solver)
{
    return model.name() == solver.modelName();
}


} // namespace PHARE
