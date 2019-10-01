
#ifndef PHARE_MESSENGER_INITIALIZER_H
#define PHARE_MESSENGER_INITIALIZER_H


#include "messengers/messenger.h"
#include "physical_models/physical_model.h"
#include "solvers/solver.h"


namespace PHARE
{
namespace solver
{
    class MessengerRegistration
    {
    public:
        /**
         * @brief registerQuantities asks the fine level model, the coarse level model and the
         * solver, to fill MessengerInfo structures. Two structures are filled, one for the finer
         * and one for the coarser level. Once filled, these two structures are given to the
         * IMessenger which now knows quantities that will need messages.
         */
        template<typename AMR_Types>
        static void registerQuantities(amr::IMessenger<IPhysicalModel<AMR_Types>>& messenger,
                                       IPhysicalModel<AMR_Types> const& coarseModel,
                                       IPhysicalModel<AMR_Types> const& fineModel,
                                       ISolver<AMR_Types> const& solver)
        {
            if (messenger.fineModelName() == fineModel.name()
                && messenger.coarseModelName() == coarseModel.name())
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
            else
            {
                auto error = std::string("Error ") + fineModel.name() + std::string{" or "}
                             + coarseModel.name() + std::string{" incompatible with "}
                             + messenger.name() + std::string{"\n"};
                throw std::runtime_error(error);
            }
        }
    };

} // namespace solver


} // namespace PHARE
#endif
