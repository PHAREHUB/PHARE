


#ifndef PHARE_MHD_MODEL_H
#define PHARE_MHD_MODEL_H

#include <string>

#include "evolution/messengers/mhd_messenger_info.h"
#include "models/mhd_state.h"
#include "physical_models/physical_model.h"
#include "tools/resources_manager.h"



namespace PHARE
{
namespace amr_interface
{
    template<typename GridLayoutT, typename VecFieldT>
    class MHDModel : public IPhysicalModel
    {
    public:
        static const std::string model_name;
        static constexpr auto dimension = GridLayoutT::dimension;
        using resources_manager_type    = ResourcesManager<GridLayoutT>;


        explicit MHDModel(std::shared_ptr<resources_manager_type> resourcesManager)
            : IPhysicalModel{model_name}
            , resourcesManager{std::move(resourcesManager)}
        {
        }

        virtual void initialize(SAMRAI::hier::Patch& patch) override {}


        virtual void allocate(SAMRAI::hier::Patch& patch, double const allocateTime) override
        {
            resourcesManager->allocate(state.B, patch, allocateTime);
            resourcesManager->allocate(state.V, patch, allocateTime);
        }



        virtual void fillMessengerInfo(std::unique_ptr<IMessengerInfo> const& info) const override
        {
        }

        virtual ~MHDModel() override = default;

        core::MHDState<VecFieldT> state;
        std::shared_ptr<resources_manager_type> resourcesManager;
    };

    template<typename GridLayoutT, typename VecFieldT>
    const std::string MHDModel<GridLayoutT, VecFieldT>::model_name = "MHDModel";


} // namespace amr_interface
} // namespace PHARE

#endif
