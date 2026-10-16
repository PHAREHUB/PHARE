#ifndef PHARE_MESSENGER_INITIALIZER_HPP
#define PHARE_MESSENGER_INITIALIZER_HPP


#include "amr/messengers/messenger.hpp"
#include "physical_models/physical_model.hpp"
#include "amr/solvers/solver.hpp"


namespace PHARE
{
namespace solver
{

    // This class is necessary to decouple models and solvers from the Messengers.
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
                                       IPhysicalModel<AMR_Types>& coarseModel,
                                       IPhysicalModel<AMR_Types>& fineModel,
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
                messenger.registerQuantities(coarseModel, fineModel);
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
