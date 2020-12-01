


#ifndef PHARE_MHD_MODEL_H
#define PHARE_MHD_MODEL_H

#include <string>


#include <SAMRAI/hier/PatchLevel.h>

#include "amr/messengers/mhd_messenger_info.h"
#include "core/models/mhd_state.h"
#include "solver/physical_models/physical_model.h"
#include "amr/resources_manager/resources_manager.h"



namespace PHARE
{
namespace solver
{
    template<typename GridLayoutT, typename VecFieldT, typename AMR_Types>
    class MHDModel : public IPhysicalModel<AMR_Types>
    {
    public:
        using patch_t   = typename AMR_Types::patch_t;
        using level_t   = typename AMR_Types::level_t;
        using Interface = IPhysicalModel<AMR_Types>;

        static const std::string model_name;
        static constexpr auto dimension = GridLayoutT::dimension;
        using resources_manager_type    = amr::ResourcesManager<GridLayoutT>;


        explicit MHDModel(std::shared_ptr<resources_manager_type> const& _resourcesManager)
            : IPhysicalModel<AMR_Types>{model_name}
            , resourcesManager{std::move(_resourcesManager)}
        {
        }

        virtual void initialize(level_t& /*level*/) override {}


        virtual void allocate(patch_t& patch, double const allocateTime) override
        {
            resourcesManager->allocate(state.B, patch, allocateTime);
            resourcesManager->allocate(state.V, patch, allocateTime);
        }



        virtual void
        fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override
        {
        }

        auto setOnPatch(patch_t& patch) { return resourcesManager->setOnPatch(patch, *this); }


        virtual ~MHDModel() override = default;

        core::MHDState<VecFieldT> state;
        std::shared_ptr<resources_manager_type> resourcesManager;
    };

    template<typename GridLayoutT, typename VecFieldT, typename AMR_Types>
    const std::string MHDModel<GridLayoutT, VecFieldT, AMR_Types>::model_name = "MHDModel";


} // namespace solver
} // namespace PHARE

#endif
