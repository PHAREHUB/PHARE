
#ifndef PHARE_MESSENGER_INITIALIZER_H
#define PHARE_MESSENGER_INITIALIZER_H


#include "evolution/messengers/messenger.h"
#include "evolution/solvers/solver.h"
#include "physical_models/physical_model.h"


namespace PHARE
{
class MessengerRegistration
{
public:
    static void registerQuantities(IMessenger& messenger, IPhysicalModel const& coarseModel,
                                   IPhysicalModel const& fineModel, ISolver const& solver)
    {
        auto fromCoarserInfo = messenger.emptyInfoFromCoarser();
        auto fromFinerInfo   = messenger.emptyInfoFromFiner();

        fineModel.fillMessengerInfo(fromFinerInfo);
        coarseModel.fillMessengerInfo(fromCoarserInfo);

        // solver only fills fromFinerInfo since
        // that's on this level it is solving equations
        solver.fillMessengerInfo(fromFinerInfo);

        messenger.registerQuantities(std::move(fromCoarserInfo), std::move(fromFinerInfo));
    }
};



} // namespace PHARE
#endif
