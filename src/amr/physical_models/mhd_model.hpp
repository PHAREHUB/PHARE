#ifndef PHARE_MHD_MODEL_HPP
#define PHARE_MHD_MODEL_HPP

#include "core/def.hpp"
#include "core/models/quantities/mhd_quantities.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep
#include "core/models/mhd_state.hpp"
#include "core/models/external_field.hpp"
#include "core/models/external_field_updater.hpp"
#include "core/models/external_field_updater_factory.hpp"

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

    using physical_quantity_type      = core::MHDQuantity;
    using vecfield_type               = VecFieldT;
    using field_type                  = vecfield_type::field_type;
    using state_type                  = core::MHDState<vecfield_type>;
    using gridlayout_type             = GridLayoutT;
    using grid_type                   = Grid_t;
    using resources_manager_type      = amr::ResourcesManager<gridlayout_type, Grid_t>;
    using external_field_type         = core::ExternalField<vecfield_type>;
    using external_field_updater_type = core::IExternalFieldUpdater<vecfield_type, gridlayout_type>;
    using external_field_factory_type
        = core::ExternalFieldUpdaterFactory<vecfield_type, gridlayout_type>;

    static constexpr std::string_view model_type_name = "MHDModel";
    static inline std::string const model_name{model_type_name};

    state_type state;
    external_field_type externalField;
    std::shared_ptr<resources_manager_type> resourcesManager;
    std::unique_ptr<external_field_updater_type> externalFieldUpdater;

    // diagnostics buffers
    vecfield_type V_diag_{"diagnostics_V_", core::MHDQuantity::Vector::V};
    field_type P_diag_{"diagnostics_P_", core::MHDQuantity::Scalar::P};

    // maybe these could have a single allocation shared for hybrid and mhd, as they are strictly
    // temporaries. Right now the hybrid version is in the hybrid_hybrid_messenger_strategy.hpp
    field_type tmpField_{"PHARE_sumField_MHD", core::MHDQuantity::Scalar::ScalarAllPrimal};
    vecfield_type tmpVec_{"PHARE_sumVec_MHD", core::MHDQuantity::Vector::VecAllPrimal};

    void initialize(level_t& level) override;


    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
        resourcesManager->allocate(V_diag_, patch, allocateTime);
        resourcesManager->allocate(P_diag_, patch, allocateTime);
        resourcesManager->allocate(tmpField_, patch, allocateTime);
        resourcesManager->allocate(tmpVec_, patch, allocateTime);
        resourcesManager->allocate(externalField, patch, allocateTime);
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
        , externalField{model_name}
        , resourcesManager{_resourcesManager}
        , externalFieldUpdater{dict.contains("external_field")
                                   ? external_field_factory_type::create(dict["external_field"])
                                   : external_field_factory_type::createZero()}
    {
        resourcesManager->registerResources(V_diag_);
        resourcesManager->registerResources(P_diag_);
        resourcesManager->registerResources(tmpField_);
        resourcesManager->registerResources(tmpVec_);
        resourcesManager->registerResources(externalField);
    }

    void initializeExternalField(level_t& level, double time) override
    {
        for (auto const& patch : resourcesManager->enumerate(level, externalField, tmpVec_))
        {
            auto const layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            auto scratch      = core::view_as(tmpVec_, vecfield_type::tensor_t::E, layout);
            (*externalFieldUpdater)(externalField, scratch, layout, time);
        }
    }

    void updateExternalField(level_t& level, double time) override
    {
        if (externalFieldUpdater->isTimeDependent())
            initializeExternalField(level, time);
    }

    ~MHDModel() override = default;


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const { return state.isUsable() and externalField.isUsable(); }

    NO_DISCARD bool isSettable() const { return state.isSettable() and externalField.isSettable(); }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state, externalField);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(state, externalField);
    }

    //-------------------------------------------------------------------------
    //                  ends the ResourcesUser interface
    //-------------------------------------------------------------------------

    std::unordered_map<std::string, std::shared_ptr<core::NdArrayVector<dimension, int>>> tags;
};

template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
void MHDModel<GridLayoutT, VecFieldT, AMR_Types, Grid_t>::initialize(level_t& level)
{
    auto& rm = *(this->resourcesManager);
    for (auto& patch : rm.enumerate(level, *this))
    {
        auto layout = amr::layoutFromPatch<GridLayoutT>(*patch);
        state.initialize(layout);
    }
}

template<typename GridLayoutT, typename VecFieldT, typename AMR_Types, typename Grid_t>
void MHDModel<GridLayoutT, VecFieldT, AMR_Types, Grid_t>::fillMessengerInfo(
    std::unique_ptr<amr::IMessengerInfo> const& info) const
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
