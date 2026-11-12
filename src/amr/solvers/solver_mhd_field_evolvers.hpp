#ifndef PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP

#include "core/numerics/time_integrator_utils.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/constrained_transport/dissipative_electric_field.hpp"
#include "core/numerics/constrained_transport/upwind_constrained_transport.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"

#include "amr/resources_manager/amr_utils.hpp"

#include "solver_field_evolvers.hpp"

namespace PHARE::solver
{

template<typename Model>
class ToConservativeTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::ToConservativeConverter<GridLayout>;

public:
    explicit ToConservativeTransformer(level_t& level, auto& model)
        : level{level}
        , model{model}
    {
    }

    void operator()(auto& state, double const gamma)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout, gamma}(state.rho, state.V, state.B, state.P, state.rhoV, state.Etot);
        }
    }

    level_t& level;
    Model& model;
};
template<typename Model>
ToConservativeTransformer(typename Model::amr_types::level_t&, Model&)
    -> ToConservativeTransformer<Model>;



template<typename Model>
class ToPrimitiveTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::ToPrimitiveConverter<GridLayout>;

public:
    explicit ToPrimitiveTransformer(level_t& level, auto& model)
        : level{level}
        , model{model}
    {
    }

    void operator()(auto& state, double const gamma)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(gamma, state.rho, state.rhoV, state.B, state.Etot, state.V, state.P);
        }
    }

    level_t& level;
    Model& model;
};
template<typename Model>
ToPrimitiveTransformer(typename Model::amr_types::level_t&, Model&)
    -> ToPrimitiveTransformer<Model>;




template<typename Model, typename FVMethod>
class FVMethodTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = FVMethod;

public:
    using info_type   = core_type::Info_t;
    using Equations_t = core_type::Equations_t;

    template<typename T>
    using Rec = core_type::template Rec<T>;

    constexpr static auto Hall = core_type::Hall;

    explicit FVMethodTransformer(level_t& level, auto& model, info_type const& info)
        : level{level}
        , model{model}
        , info{info}
    {
    }


    template<typename MHDState>
    void operator()(auto& ct_state, auto& dissipative_electric_state, MHDState& state, auto& fluxes)
        requires(not core::has_tiled_field_type_c<typename MHDState::vecfield_type>)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, ct_state, dissipative_electric_state, state, fluxes))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type finite_volume_method{info, layout};
            finite_volume_method(ct_state, dissipative_electric_state, state, fluxes);
        }
    }

    // fluxes and ct coefficients are computed per tile over a box grown by the reconstruction
    // ghosts, so no inner ghost sync of the outputs
    template<typename MHDState>
    void operator()(auto& ct_state, auto& dissipative_electric_state, MHDState& state, auto& fluxes)
        requires(core::has_tiled_field_type_c<typename MHDState::vecfield_type>)
    {
        auto accessor = amr::make_model_level_accessor(level, model, ct_state,
                                                       dissipative_electric_state, state, fluxes);
        tiled_level_exec(
            accessor, [](auto& view) -> auto& { return std::get<2>(view.args).rho; },
            [&](auto& view, auto const tile_idx) {
                auto& [ct_v, dissipative_v, state_v, fluxes_v] = view.args;
                auto const& layout = core::tile_layout(state_v.rho, tile_idx);
                auto ct_t          = core::tile_at(ct_v, tile_idx);
                auto dissipative_t = core::tile_at(dissipative_v, tile_idx);
                auto state_t       = core::tile_at(state_v, tile_idx);
                auto fluxes_t      = core::tile_at(fluxes_v, tile_idx);
                core_type{info, layout}(ct_t, dissipative_t, state_t, fluxes_t);
            });
    }

    level_t& level;
    Model& model;
    info_type const& info;
};



template<typename Model>
class FiniteVolumeEulerTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::FiniteVolumeEuler<GridLayout>;

public:
    explicit FiniteVolumeEulerTransformer(level_t& level, auto& model)
        : level{level}
        , model{model}
    {
    }

    void operator()(Model::state_type& state, Model::state_type& statenew, auto& fluxes,
                    double const dt)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, state, statenew, fluxes))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(state, statenew, fluxes, dt);
        }
    }


    level_t& level;
    Model& model;
};
template<typename Model>
FiniteVolumeEulerTransformer(typename Model::amr_types::level_t&, Model&)
    -> FiniteVolumeEulerTransformer<Model>;




template<typename Model>
class DissipativeElectricFieldTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::DissipativeElectricField<GridLayout>;

public:
    using info_type = core_type::Info_t;

    explicit DissipativeElectricFieldTransformer(level_t& level, Model& model,
                                                 info_type const& info)
        : level{level}
        , model{model}
        , info{info}
    {
    }

    template<typename DissipativeState, typename MHDState>
    void operator()(DissipativeState& dissipative_electric_state, MHDState& mhd_state)
        requires(not core::has_tiled_field_type_c<typename MHDState::vecfield_type>)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, dissipative_electric_state, mhd_state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{info, layout}(dissipative_electric_state, mhd_state);
        }
    }

    template<typename DissipativeState, typename MHDState>
    void operator()(DissipativeState& dissipative_electric_state, MHDState& mhd_state)
        requires(core::has_tiled_field_type_c<typename MHDState::vecfield_type>)
    {
        auto accessor
            = amr::make_model_level_accessor(level, model, dissipative_electric_state, mhd_state);
        tiled_level_transform(
            accessor, [](auto& view) -> auto& { return std::get<0>(view.args).E(); },
            [&](auto& view, auto const tile_idx) {
                auto& [dissipative_v, mhd_v] = view.args;
                auto const& layout           = core::tile_layout(mhd_v.B, tile_idx);
                auto dissipative_t           = core::tile_at(dissipative_v, tile_idx);
                auto const mhd_t             = core::tile_at(mhd_v, tile_idx);
                using Layout_t               = std::decay_t<decltype(layout)>;
                core::DissipativeElectricField<Layout_t>{info, layout}(dissipative_t, mhd_t);
            });
    }

    level_t& level;
    Model& model;
    info_type const& info;
};




template<typename GridLayout, typename Model, template<typename> typename Reconstruction, auto Hall>
class ConstrainedTransportTransformer
{
    using level_t   = Model::amr_types::level_t;
    using core_type = core::UpwindConstrainedTransport<GridLayout, Reconstruction, Hall>;

public:
    using info_type = core_type::Info_t;

    explicit ConstrainedTransportTransformer(level_t& level, auto& model, info_type const& info)
        : level{level}
        , model{model}
        , info{info}
    {
    }

    template<typename MHDState>
    void operator()(auto& ct_state, auto& dissipative_electric_state, MHDState& mhd_state)
        requires(not core::has_tiled_field_type_c<typename MHDState::vecfield_type>)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, ct_state, dissipative_electric_state, mhd_state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type constrained_transport_{info, layout};
            constrained_transport_(ct_state, dissipative_electric_state, mhd_state);
        }
    }

    template<typename MHDState>
    void operator()(auto& ct_state, auto& dissipative_electric_state, MHDState& mhd_state)
        requires(core::has_tiled_field_type_c<typename MHDState::vecfield_type>)
    {
        auto accessor = amr::make_model_level_accessor(level, model, ct_state,
                                                       dissipative_electric_state, mhd_state);
        tiled_level_transform(
            accessor, [](auto& view) -> auto& { return std::get<2>(view.args).E; },
            [&](auto& view, auto const tile_idx) {
                auto& [ct_v, dissipative_v, mhd_v] = view.args;
                auto const& layout                 = core::tile_layout(mhd_v.E, tile_idx);
                auto ct_t                          = core::tile_at(ct_v, tile_idx);
                auto dissipative_t                 = core::tile_at(dissipative_v, tile_idx);
                auto mhd_t                         = core::tile_at(mhd_v, tile_idx);
                core_type{info, layout}(ct_t, dissipative_t, mhd_t);
            });
    }


    level_t& level;
    Model& model;
    info_type const info;
};




template<typename Model>
class RKUtilsTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;
    using core_type  = core::RKUtils<GridLayout>;

public:
    void operator()(Model::state_type& res, auto... pairs)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, res, pairs.state...))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(res, pairs...);
        }
    }


    level_t& level;
    Model& model;
};


template<typename Model>
struct Dispatchers : FieldEvolverDispatchers<Model>
{
    using GridLayout = Model::gridlayout_type;

    // Ampere_t and Faraday_t inherited from FieldEvolverDispatchers<Model>

    using ToPrimitiveConverter_t    = ToPrimitiveTransformer<Model>;
    using ToConservativeConverter_t = ToConservativeTransformer<Model>;

    template<typename FVMethodStrategy>
    using FVMethod_t = FVMethodTransformer<Model, FVMethodStrategy>;

    using FiniteVolumeEuler_t = FiniteVolumeEulerTransformer<Model>;

    using DissipativeElectricField_t = DissipativeElectricFieldTransformer<Model>;

    template<template<typename> typename Reconstruction, auto Hall>
    using ConstrainedTransport_t
        = ConstrainedTransportTransformer<GridLayout, Model, Reconstruction, Hall>;

    using RKUtils_t = RKUtilsTransformer<Model>;
};



}; // namespace PHARE::solver

#endif // PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
