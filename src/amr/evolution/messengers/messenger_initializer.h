
#ifndef PHARE_MESSENGER_INITIALIZER_H
#define PHARE_MESSENGER_INITIALIZER_H


#include "evolution/messengers/messenger.h"
#include "evolution/solvers/solver.h"
#include "physical_models/physical_model.h"


namespace PHARE
{
class MessengerInitializer
{
public:
    static void setup(IMessenger& messenger, IPhysicalModel const& coarseModel,
                      IPhysicalModel const& fineModel, ISolver const& solver)
    {
        auto fromCoarserInfo = messenger.emptyInfoFromCoarser();
        auto fromFinerInfo   = messenger.emptyInfoFromFiner();

        fineModel.fillMessengerInfo(fromFinerInfo);
        coarseModel.fillMessengerInfo(fromCoarserInfo);

        messenger.registerQuantities(std::move(fromCoarserInfo), std::move(fromFinerInfo));
    }
};



} // namespace PHARE
#endif
