#ifndef PHARE_MHD_MODEL_HPP
#define PHARE_MHD_MODEL_HPP

#include "core/def.hpp"
#include "core/def/phare_mpi.hpp" // IWYU pragma: keep
#include "core/models/mhd_state.hpp"
#include "core/data/vecfield/vecfield_initializer.hpp"

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
    using resources_manager_type = amr::ResourcesManager<gridlayout_type, Grid_t>;

    static constexpr std::string_view model_type_name = "MHDModel";
    static inline std::string const model_name{model_type_name};

    state_type state;
    // The static background field B0 of the B = B0 + B1 split: a single face-centered vector field
    // (same Yee centering as B1). Held once and shared by every RK stage (B0 does not change
    // between stages, so it is never copied around). B0 is re-evaluated analytically per patch
    // (never interpolated across AMR levels); the numerics linearly interpolate it from this native
    // face-centered storage onto the EMF edges / Riemann faces where they consume it.
    vecfield_type B0{model_name + "_B0", core::MHDQuantity::Vector::B0};
    core::VecFieldInitializer<dimension> B0init_;
    std::shared_ptr<resources_manager_type> resourcesManager;

    // diagnostics buffers
    vecfield_type V_diag_{"diagnostics_V_", core::MHDQuantity::Vector::V};
    field_type P_diag_{"diagnostics_P_", core::MHDQuantity::Scalar::P};
    // total fields reconstructed from the B0 + B1 split for output
    vecfield_type BTotal_diag_{"diagnostics_BTotal_", core::MHDQuantity::Vector::B};
    field_type EtotTotal_diag_{"diagnostics_EtotTotal_", core::MHDQuantity::Scalar::Etot};
    field_type divB_diag_{"diagnostics_divB_", core::MHDQuantity::Scalar::divB};

    // maybe these could have a single allocation shared for hybrid and mhd, as they are strictly
    // temporaries. Right now the hybrid version is in the hybrid_hybrid_messenger_strategy.hpp
    field_type tmpField_{"PHARE_sumField_MHD", core::MHDQuantity::Scalar::ScalarAllPrimal};
    vecfield_type tmpVec_{"PHARE_sumVec_MHD", core::MHDQuantity::Vector::VecAllPrimal};

    void initialize(level_t& level) override;

    // (Re-)evaluate the static background B0 analytically on every patch of a level. B0 is never
    // interpolated across levels: each level (including levels created/changed by a regrid)
    // evaluates B0 at its own coordinates. Called after init/regrid of non-root levels.
    void updateExternalFields(level_t& level)
    {
        for (auto& patch : level)
        {
            auto layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto _      = resourcesManager->setOnPatch(*patch, B0);
            B0init_.initialize(B0, layout);
        }
    }


    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
        resourcesManager->allocate(B0, patch, allocateTime);
        resourcesManager->allocate(V_diag_, patch, allocateTime);
        resourcesManager->allocate(P_diag_, patch, allocateTime);
        resourcesManager->allocate(BTotal_diag_, patch, allocateTime);
        resourcesManager->allocate(EtotTotal_diag_, patch, allocateTime);
        resourcesManager->allocate(divB_diag_, patch, allocateTime);
        resourcesManager->allocate(tmpField_, patch, allocateTime);
        resourcesManager->allocate(tmpVec_, patch, allocateTime);
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
        , B0init_{dict["mhd_state"]["external_magnetic"]["initializer"]}
        , resourcesManager{std::move(_resourcesManager)}
    {
        resourcesManager->registerResources(B0);
        resourcesManager->registerResources(V_diag_);
        resourcesManager->registerResources(P_diag_);
        resourcesManager->registerResources(BTotal_diag_);
        resourcesManager->registerResources(EtotTotal_diag_);
        resourcesManager->registerResources(divB_diag_);
        resourcesManager->registerResources(tmpField_);
        resourcesManager->registerResources(tmpVec_);
    }

    ~MHDModel() override = default;


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable() and B0.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable() and B0.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state, B0);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() { return std::forward_as_tuple(state, B0); }

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
        auto _      = this->resourcesManager->setOnPatch(*patch, state, B0);

        // evaluate the static background B0 first, then the dynamic state, which subtracts B0 from
        // the prescribed total field to form B1.
        B0init_.initialize(B0, layout);
        state.initialize(layout, B0);
    }
}

template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
void MHDModel<GridLayoutT, VecFieldT, AMR_Types, Grid_t>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
{
    auto& MHDInfo = dynamic_cast<amr::MHDMessengerInfo&>(*info);

    MHDInfo.modelDensity     = state.rho.name();
    MHDInfo.modelVelocity    = state.V.name();
    MHDInfo.modelMagnetic    = state.B1.name();
    MHDInfo.modelPressure    = state.P.name();
    MHDInfo.modelMomentum    = state.rhoV.name();
    MHDInfo.modelTotalEnergy = state.Etot1.name();
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
