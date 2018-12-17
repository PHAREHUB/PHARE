

#include "messenger.h"


namespace PHARE
{
bool areCompatible(IMessenger const& messenger, IPhysicalModel const& model)
{
    return model.name() == messenger.fineModelName() || model.name() == messenger.coarseModelName();
}
} // namespace PHARE
