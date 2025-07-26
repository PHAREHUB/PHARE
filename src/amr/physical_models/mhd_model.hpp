#ifndef PHARE_MHD_MODEL_HPP
#define PHARE_MHD_MODEL_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp"
#include "core/models/mhd_state.hpp"

#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/resources_manager.hpp"

#include <SAMRAI/hier/PatchLevel.h>

#include <string>
#include <string_view>


namespace PHARE::solver
{
template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
class MHDModel : public IPhysicalModel<AMR_Types>
{
public:
    static constexpr auto dimension = GridLayoutT::dimension;

    using type_list = PHARE::core::type_list<GridLayoutT, VecFieldT, AMR_Types, Grid_t>;
    using amr_types = AMR_Types;
    using patch_t   = amr_types::patch_t;
    using level_t   = amr_types::level_t;
    using Interface = IPhysicalModel<AMR_Types>;

    using vecfield_type          = VecFieldT;
    using field_type             = vecfield_type::field_type;
    using state_type             = core::MHDState<vecfield_type>;
    using gridlayout_type        = GridLayoutT;
    using grid_type              = Grid_t;
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, Grid_t>;

    static constexpr std::string_view model_type_name = "MHDModel";
    static inline std::string const model_name{model_type_name};

    state_type state;
    std::shared_ptr<resources_manager_type> resourcesManager;

    void initialize(level_t& level) override;


    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
    }

    void fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const override;

    NO_DISCARD auto setOnPatch(patch_t& patch)
    {
        return resourcesManager->setOnPatch(patch, *this);
    }

    explicit MHDModel(PHARE::initializer::PHAREDict const& dict,
                      std::shared_ptr<resources_manager_type> const& _resourcesManager)
        : IPhysicalModel<AMR_Types>{model_name}
        , state{dict["mhd_state"]}
        , resourcesManager{std::move(_resourcesManager)}
    {
    }

    ~MHDModel() override = default;

    auto get_B() -> auto& { return state.B; }

    auto get_B() const -> auto& { return state.B; }

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

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
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

template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
void MHDModel<GridLayoutT, VecFieldT, AMR_Types, Grid_t>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& MHDInfo = dynamic_cast<amr::MHDMessengerInfo&>(*info);

    MHDInfo.modelDensity     = state.rho.name();
    MHDInfo.modelVelocity    = core::VecFieldNames{state.V};
    MHDInfo.modelMagnetic    = core::VecFieldNames{state.B};
    MHDInfo.modelPressure    = state.P.name();
    MHDInfo.modelMomentum    = core::VecFieldNames{state.rhoV};
    MHDInfo.modelTotalEnergy = state.Etot.name();
    MHDInfo.modelElectric    = core::VecFieldNames{state.E};
    MHDInfo.modelCurrent     = core::VecFieldNames{state.J};

    MHDInfo.initDensity.push_back(MHDInfo.modelDensity);
    MHDInfo.initMomentum.push_back(MHDInfo.modelMomentum);
    MHDInfo.initMagnetic.push_back(MHDInfo.modelMagnetic);
    MHDInfo.initTotalEnergy.push_back(MHDInfo.modelTotalEnergy);

    MHDInfo.ghostDensity.push_back(MHDInfo.modelDensity);
    MHDInfo.ghostVelocity.push_back(MHDInfo.modelVelocity);
    MHDInfo.ghostMagnetic.push_back(MHDInfo.modelMagnetic);
    MHDInfo.ghostPressure.push_back(MHDInfo.modelPressure);
    MHDInfo.ghostMomentum.push_back(MHDInfo.modelMomentum);
    MHDInfo.ghostTotalEnergy.push_back(MHDInfo.modelTotalEnergy);
    MHDInfo.ghostElectric.push_back(MHDInfo.modelElectric);
    MHDInfo.ghostCurrent.push_back(MHDInfo.modelCurrent);
}




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
