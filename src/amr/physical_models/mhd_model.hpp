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
    // Vector-potential init: when set, B0 = curl(A0) is built from a full 3D vector potential with
    // the discrete curl (div B0 = 0 to machine precision) instead of from the component functions
    // B0init_. Defaults keep the legacy component-wise init.
    bool b0FromPotential_ = false;
    core::VecFieldInitializer<dimension> a0Init_;
    // Time-dependent external field B0(x,t): when enabled the background is re-stamped once per
    // full timestep and the split gains the sources -dB0/dt (on B1) and -dB0/dt . B1 (on the
    // reduced energy). dB0dt holds the analytic time derivative (same face-centering as B0), frozen
    // across the RK stages of a step. The ST initializers evaluate the user f(x,t) at a given time;
    // in potential mode B0 = curl(A0(t)) and dB0dt = curl(dA0/dt(t)) share the discrete curl.
    bool b0TimeDependent_ = false;
    vecfield_type dB0dt{model_name + "_dB0dt", core::MHDQuantity::Vector::B0};
    // Edge-centered scratch holding the vector potential A while building B0 = curl(A) with the
    // discrete curl. Dedicated (not state.E) so the potential restamp never clobbers the electric
    // field: updateExternalFields runs every step for a time-dependent B0, and state.E must survive
    // untouched for that step's constrained transport.
    vecfield_type B0Ascratch_{model_name + "_B0Ascratch", core::MHDQuantity::Vector::E};
    core::SpaceTimeVecFieldInitializer<dimension> B0initST_;
    core::SpaceTimeVecFieldInitializer<dimension> dB0dtInitST_;
    core::SpaceTimeVecFieldInitializer<dimension> a0InitST_;
    core::SpaceTimeVecFieldInitializer<dimension> da0dtInitST_;
    std::shared_ptr<resources_manager_type> resourcesManager;

    // diagnostics buffers
    vecfield_type V_diag_{"diagnostics_V_", core::MHDQuantity::Vector::V};
    field_type P_diag_{"diagnostics_P_", core::MHDQuantity::Scalar::P};
    // The total-field / total-energy / divB output buffers of the B0 + B1 split live on the
    // diagnostic ModelView (diagnostic_model_view.hpp), which owns the reconstruction; the model
    // holds none of them here (they would be dead and would collide on the SAMRAI variable name).

    // maybe these could have a single allocation shared for hybrid and mhd, as they are strictly
    // temporaries. Right now the hybrid version is in the hybrid_hybrid_messenger_strategy.hpp
    field_type tmpField_{"PHARE_sumField_MHD", core::MHDQuantity::Scalar::ScalarAllPrimal};
    vecfield_type tmpVec_{"PHARE_sumVec_MHD", core::MHDQuantity::Vector::VecAllPrimal};

    void initialize(level_t& level) override;

    // (Re-)evaluate the static background B0 analytically on every patch of a level. B0 is never
    // interpolated across levels: each level (including levels created/changed by a regrid)
    // evaluates B0 at its own coordinates. Called after init/regrid of non-root levels.
    void updateExternalFields(level_t& level, double time = 0.)
    {
        for (auto& patch : level)
        {
            auto layout = amr::layoutFromPatch<GridLayoutT>(*patch);
            if (b0TimeDependent_)
            {
                if (b0FromPotential_)
                {
                    // B0 = curl(A0(t)) and dB0/dt = curl(dA0/dt(t)) via the same discrete curl,
                    // both using the dedicated B0Ascratch_ vecfield as the A scratch.
                    auto _ = resourcesManager->setOnPatch(*patch, B0, dB0dt, B0Ascratch_);
                    core::initBFromPotential(a0InitST_, B0, B0Ascratch_, layout, time);
                    core::initBFromPotential(da0dtInitST_, dB0dt, B0Ascratch_, layout, time);
                }
                else
                {
                    auto _ = resourcesManager->setOnPatch(*patch, B0, dB0dt);
                    B0initST_.initialize(B0, layout, time);
                    dB0dtInitST_.initialize(dB0dt, layout, time);
                }
            }
            else if (b0FromPotential_)
            {
                // B0 = curl(A0), built with the discrete curl using B0Ascratch_ as the A scratch.
                auto _ = resourcesManager->setOnPatch(*patch, B0, B0Ascratch_);
                core::initBFromPotential(a0Init_, B0, B0Ascratch_, layout);
            }
            else
            {
                auto _ = resourcesManager->setOnPatch(*patch, B0);
                B0init_.initialize(B0, layout);
            }
        }
    }


    void allocate(patch_t& patch, double const allocateTime) override
    {
        resourcesManager->allocate(state, patch, allocateTime);
        resourcesManager->allocate(B0, patch, allocateTime);
        resourcesManager->allocate(dB0dt, patch, allocateTime);
        resourcesManager->allocate(B0Ascratch_, patch, allocateTime);
        resourcesManager->allocate(V_diag_, patch, allocateTime);
        resourcesManager->allocate(P_diag_, patch, allocateTime);
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
        , b0FromPotential_{cppdict::get_value(dict["mhd_state"], "b0_init_mode",
                                              std::string{"components"})
                           == "potential"}
        , b0TimeDependent_{cppdict::get_value(dict["mhd_state"]["external_magnetic"],
                                              "time_dependent", false)}
        , resourcesManager{std::move(_resourcesManager)}
    {
        // Static component B0: the initializer components are InitFunctions. In time-dependent mode
        // they are SpaceTimeFunctions instead (read into B0initST_ below), so B0init_ must not be
        // built from them then (to<InitFunction> would throw "invalid type").
        if (not b0TimeDependent_)
            B0init_ = core::VecFieldInitializer<dimension>{
                dict["mhd_state"]["external_magnetic"]["initializer"]};

        // Vector-potential init: read the vector potential only in "potential" mode so dicts that
        // predate this feature need not provide the potential key.
        if (b0FromPotential_ and not b0TimeDependent_)
            a0Init_ = core::VecFieldInitializer<dimension>{
                dict["mhd_state"]["external_magnetic"]["initializer"]["potential"]};

        // Time-dependent B0: read the space+time initializers and the analytic derivative for the
        // active mode only, so static dicts (no derivative subtree) are unaffected.
        if (b0TimeDependent_)
        {
            auto const& extB = dict["mhd_state"]["external_magnetic"];
            if (b0FromPotential_)
            {
                a0InitST_ = core::SpaceTimeVecFieldInitializer<dimension>{
                    extB["initializer"]["potential"]};
                da0dtInitST_ = core::SpaceTimeVecFieldInitializer<dimension>{
                    extB["derivative"]["potential"]};
            }
            else
            {
                B0initST_ = core::SpaceTimeVecFieldInitializer<dimension>{extB["initializer"]};
                dB0dtInitST_ = core::SpaceTimeVecFieldInitializer<dimension>{extB["derivative"]};
            }
        }

        resourcesManager->registerResources(B0);
        resourcesManager->registerResources(dB0dt);
        resourcesManager->registerResources(B0Ascratch_);
        resourcesManager->registerResources(V_diag_);
        resourcesManager->registerResources(P_diag_);
        resourcesManager->registerResources(tmpField_);
        resourcesManager->registerResources(tmpVec_);
    }

    ~MHDModel() override = default;


    //-------------------------------------------------------------------------
    //                  start the ResourcesUser interface
    //-------------------------------------------------------------------------

    NO_DISCARD bool isUsable() const
    {
        return state.isUsable() and B0.isUsable() and dB0dt.isUsable();
    }

    NO_DISCARD bool isSettable() const
    {
        return state.isSettable() and B0.isSettable() and dB0dt.isSettable();
    }

    NO_DISCARD auto getCompileTimeResourcesViewList() const
    {
        return std::forward_as_tuple(state, B0, dB0dt);
    }

    NO_DISCARD auto getCompileTimeResourcesViewList()
    {
        return std::forward_as_tuple(state, B0, dB0dt);
    }

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
        auto _ = this->resourcesManager->setOnPatch(*patch, state, B0, dB0dt, B0Ascratch_);

        // evaluate the background B0 first, then the dynamic state, which subtracts B0 from the
        // prescribed total field to form B1. The initial stamp is at time 0: the Python total field
        // was assembled with B0(.,0), so B1 = total - B0(.,0). The solver re-stamps B0/dB0dt to the
        // new time at the end of each timestep.
        if (b0TimeDependent_)
        {
            if (b0FromPotential_)
            {
                core::initBFromPotential(a0InitST_, B0, B0Ascratch_, layout, 0.);
                core::initBFromPotential(da0dtInitST_, dB0dt, B0Ascratch_, layout, 0.);
            }
            else
            {
                B0initST_.initialize(B0, layout, 0.);
                dB0dtInitST_.initialize(dB0dt, layout, 0.);
            }
        }
        else
        {
            if (b0FromPotential_)
                // B0 = curl(A0) via the discrete curl, using B0Ascratch_ as the A scratch.
                core::initBFromPotential(a0Init_, B0, B0Ascratch_, layout);
            else
                B0init_.initialize(B0, layout);
            // dB0/dt is unused for a static B0; zero it so the allocated field is well-defined.
            for (auto const& component : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& dc = dB0dt(component);
                layout.evalOnGhostBox(dc, [&](auto&... args) mutable { dc(args...) = 0.0; });
            }
        }
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
