

#include "solver.h"

namespace PHARE
{
namespace solver
{
    bool areCompatible(amr::IMessenger<IPhysicalModel> const& messenger,
                       ISolver const& solver)
    {
        return solver.modelName() == messenger.fineModelName()
               || solver.modelName() == messenger.coarseModelName();
    }


    /**
     * @brief areCompatible returns true if the model name is equal to the solver modelname
     */
    bool areCompatible(IPhysicalModel const& model, ISolver const& solver)
    {
        return model.name() == solver.modelName();
    }
} // namespace solver
} // namespace PHARE
