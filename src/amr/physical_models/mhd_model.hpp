#ifndef PHARE_MHD_MODEL_HPP
#define PHARE_MHD_MODEL_HPP

#include "core/def.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/models/mhd_state.hpp"
#include "core/utilities/variants.hpp"

#include "amr/messengers/mhd_messenger_info.hpp"
#include "amr/physical_models/physical_model.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include <SAMRAI/hier/PatchLevel.h>

#include <string>
#include <string_view>
#include <variant>


namespace PHARE::solver
{
template<typename Types>
class MHDModel : public IPhysicalModel<typename Types::amr_types>
{
public:
    using GridLayoutT        = Types::GridLayout_t;
    using VecFieldT          = Types::VecField_t;
    using AMR_Types          = Types::amr_types;
    using Grid_t             = Types::Grid_t;
    using ResourcesManager_t = Types::ResourcesManager_t;

    static constexpr auto dimension = GridLayoutT::dimension;

    using amr_types = AMR_Types;
    using patch_t   = amr_types::patch_t;
    using level_t   = amr_types::level_t;
    using Interface = IPhysicalModel<AMR_Types>;

    using physical_quantity_type = core::MHDQuantity;
    using vecfield_type          = VecFieldT;
    using field_type             = vecfield_type::field_type;
    using state_type             = core::MHDState<vecfield_type>;
    using gridlayout_type        = GridLayoutT;
    using grid_type              = Grid_t;
    using resources_manager_type = ResourcesManager_t;
    using Resources              = std::variant<field_type, vecfield_type>;

    static constexpr std::string_view model_type_name = "MHDModel";
    static inline std::string const model_name{model_type_name};

    void initialize(level_t& level) override;

    void registerResources() override { resourcesManager->registerResources(*this); }

    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(*this, patch, allocateTime);
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
        , resourcesManager{_resourcesManager}
    {
        registerTemporaryRequirements();
    }

    ~MHDModel() override = default;


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const { return std::forward_as_tuple(state); }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state); }

    void registerTemporaryRequirements(auto&&... quantities);

    auto& getTmp(auto&& as) { return core::get_from_variants(resources, as); }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    state_type state;
    std::shared_ptr<resources_manager_type> resourcesManager;
    std::vector<Resources> resources;
    std::array<std::size_t, std::variant_size_v<Resources>> n_per_type = {0};
    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
};

template<typename Types>
void MHDModel<Types>::registerTemporaryRequirements(auto&&... quantities)
{
    // todo figure out at runtime
    std::array<std::size_t, std::variant_size_v<Resources>> _n_per_type{2, 2};

    auto const namer = []<auto i>(auto n) -> std::string {
        return "tmp" + std::array<std::string, 3>{"Field", "VecField", "TensorField"}[i]
               + std::to_string(n);
    };

    core::for_N<std::variant_size_v<Resources>>([&](auto i) {
        while (n_per_type[i] < _n_per_type[i])
        {
            auto const name = namer.template operator()<i>(n_per_type[i]);
            if constexpr (i == 0)
                resources.emplace_back(field_type{name, core::MHDQuantity::all_primal_field});
            else if constexpr (i == 1)
                resources.emplace_back(vecfield_type{name, core::MHDQuantity::Vector::V});
            else
                static_assert(core::dependent_false_v<Types>);

            ++n_per_type[i];
        }
    });
}


template<typename Types>
void MHDModel<Types>::initialize(level_t& level)
{
    for (auto& patch : level)
    {
        auto layout = amr::layoutFromPatch<GridLayoutT>(*patch);
        auto _      = this->resourcesManager->setOnPatch(*patch, state);

        state.initialize(layout);
    }
}

template<typename Types>
void MHDModel<Types>::fillMessengerInfo(std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& MHDInfo = dynamic_cast<amr::MHDMessengerInfo&>(*info);

    MHDInfo.modelDensity     = state.rho.name();
    MHDInfo.modelVelocity    = state.V.name();
    MHDInfo.modelMagnetic    = state.B.name();
    MHDInfo.modelPressure    = state.P.name();
    MHDInfo.modelMomentum    = state.rhoV.name();
    MHDInfo.modelTotalEnergy = state.Etot.name();
    MHDInfo.modelElectric    = state.E.name();
    MHDInfo.modelCurrent     = state.J.name();

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
