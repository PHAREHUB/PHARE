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
        amr::onLevel(amr::make_model_level_accessor(level, model, state), [&](auto& view) {
            auto& [state_v] = view.args;
            core_type{view.layout, gamma}(state_v.rho, state_v.V, state_v.B, state_v.P,
                                          state_v.rhoV, state_v.Etot);
        });
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
        amr::onLevel(amr::make_model_level_accessor(level, model, state), [&](auto& view) {
            auto& [state_v] = view.args;
            core_type{view.layout}(gamma, state_v.rho, state_v.rhoV, state_v.B, state_v.Etot,
                                   state_v.V, state_v.P);
        });
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


    void operator()(auto& ct_state, auto& dissipative_electric_state, auto& state, auto& fluxes)
    {
        amr::onLevel(amr::make_model_level_accessor(level, model, ct_state,
                                                    dissipative_electric_state, state, fluxes),
                     [&](auto& view) {
                         auto& [ct_state_v, dissipative_v, state_v, fluxes_v] = view.args;
                         core_type finite_volume_method{info, view.layout};
                         finite_volume_method(ct_state_v, dissipative_v, state_v, fluxes_v);
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
        amr::onLevel(amr::make_model_level_accessor(level, model, state, statenew, fluxes),
                     [&](auto& view) {
                         auto& [state_v, statenew_v, fluxes_v] = view.args;
                         core_type{view.layout}(state_v, statenew_v, fluxes_v, dt);
                     });
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

    void operator()(auto& dissipative_electric_state, auto& mhd_state)
    {
        amr::onLevel(
            amr::make_model_level_accessor(level, model, dissipative_electric_state, mhd_state),
            [&](auto& view) {
                auto& [dissipative_v, mhd_state_v] = view.args;
                core_type{info, view.layout}(dissipative_v, mhd_state_v);
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

    void operator()(auto& ct_state, auto& dissipative_electric_state, auto& mhd_state)
    {
        amr::onLevel(amr::make_model_level_accessor(level, model, ct_state,
                                                    dissipative_electric_state, mhd_state),
                     [&](auto& view) {
                         auto& [ct_state_v, dissipative_v, mhd_state_v] = view.args;
                         core_type constrained_transport_{info, view.layout};
                         constrained_transport_(ct_state_v, dissipative_v, mhd_state_v);
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
        amr::onLevel(
            amr::make_model_level_accessor(level, model, res, pairs.state...), [&](auto& view) {
                std::apply(
                    [&](auto& res_v, auto&... states_v) {
                        core_type{view.layout}(res_v, decltype(pairs){pairs.weight, states_v}...);
                    },
                    view.args);
            });
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
