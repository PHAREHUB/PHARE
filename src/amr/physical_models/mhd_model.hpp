#ifndef PHARE_MHD_MODEL_HPP
#define PHARE_MHD_MODEL_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/models/mhd_state.hpp"

#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/resources_manager.hpp"

#include <SAMRAI/hier/PatchLevel.h>

#include <string>
#include <string_view>


namespace PHARE
{
namespace solver
{
    template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
    class MHDModel : public IPhysicalModel<AMR_Types>
    {
    public:
        using patch_t   = typename AMR_Types::patch_t;
        using level_t   = typename AMR_Types::level_t;
        using Interface = IPhysicalModel<AMR_Types>;

        static constexpr std::string_view model_type_name = "MHDModel";
        static inline std::string const model_name{model_type_name};

        static constexpr auto dimension = GridLayoutT::dimension;
        using resources_manager_type    = amr::ResourcesManager<GridLayoutT, Grid_t>;


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

        NO_DISCARD auto setOnPatch(patch_t& patch)
        {
            return resourcesManager->setOnPatch(patch, *this);
        }


        virtual ~MHDModel() override = default;

        core::MHDState<VecFieldT> state;
        std::shared_ptr<resources_manager_type> resourcesManager;
    };



} // namespace solver
} // namespace PHARE


namespace PHARE::solver
{

template<typename Model>
auto constexpr is_mhd_model(Model* m) -> decltype(m->model_type_name, bool())
{
    return Model::model_type_name == "MHDModel";
}

template<typename... Args>
auto constexpr is_mhd_model(Args...)
{
    return false;
}

template<typename Model>
auto constexpr is_mhd_model_v = is_mhd_model(static_cast<Model*>(nullptr));


} // namespace PHARE::solver

#endif
