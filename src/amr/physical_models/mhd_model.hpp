#ifndef PHARE_MHD_MODEL_HPP
#define PHARE_MHD_MODEL_HPP

#include <string>

#include "core/def/phare_mpi.hpp"


#include <SAMRAI/hier/PatchLevel.h>

#include "amr/messengers/mhd_messenger_info.hpp"
#include "core/models/mhd_state.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/resources_manager.hpp"
#include "core/def.hpp"



namespace PHARE::solver
{
template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
class MHDModel : public IPhysicalModel<AMR_Types>
{
public:
    using patch_t   = typename AMR_Types::patch_t;
    using level_t   = typename AMR_Types::level_t;
    using Interface = IPhysicalModel<AMR_Types>;

    using field_type      = typename VecFieldT::field_type;
    using gridlayout_type = GridLayoutT;

    static inline std::string const model_name = "MHDModel";
    static constexpr auto dimension            = gridlayout_type::dimension;
    using resources_manager_type               = amr::ResourcesManager<gridlayout_type, Grid_t>;

    core::MHDState<VecFieldT> state;
    std::shared_ptr<resources_manager_type> resourcesManager;

    virtual void initialize(level_t& level) override;


    virtual void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
    }



    virtual void
    fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& /*info*/) const override
    {
    }

    NO_DISCARD auto setOnPatch(patch_t& patch)
    {
        return resourcesManager->setOnPatch(patch, *this);
    }

    explicit MHDModel(PHARE::initializer::PHAREDict const& dict,
                      std::shared_ptr<resources_manager_type> const& _resourcesManager)
        : IPhysicalModel<AMR_Types>{model_name}
        , state{dict}
        , resourcesManager{std::move(_resourcesManager)}
    {
    }

    virtual ~MHDModel() override = default;


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(state); }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------
};

template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
void MHDModel<GridLayoutT, VecFieldT, AMR_Types, Grid_t>::initialize(level_t& level)
{
    for (auto& patch : level)
    {
        auto layout = amr::layoutFromPatch<GridLayoutT>(*patch);
        auto _      = this->resourcesManager->setOnPatch(*patch, state);

        state.initialize(layout);
    }
    resourcesManager->registerForRestarts(*this);
}

} // namespace PHARE::solver

#endif
