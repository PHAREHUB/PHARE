

#include "messenger.h"


namespace PHARE
{
namespace amr_interface
{
    bool areCompatible(IMessenger const& messenger, IPhysicalModel const& model)
    {
        return model.name() == messenger.fineModelName()
               || model.name() == messenger.coarseModelName();
    }
} // namespace amr_interface
} // namespace PHARE
