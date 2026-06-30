#ifndef PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP


#include "core/numerics/time_integrator_utils.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
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

    void operator()(auto& state, double const gamma, double const newTime)
    {
        TimeSetter setTime{level, model, newTime};

        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout, gamma}(state.rho, state.V, state.B, state.P, state.rhoV, state.Etot);
        }

        setTime(state.rho, state.V, state.P, state.rhoV, state.Etot);
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

    void operator()(auto& state, double const gamma, double const newTime)
    {
        TimeSetter setTime{level, model, newTime};

        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(gamma, state.rho, state.rhoV, state.B, state.Etot, state.V, state.P);
        }

        setTime(state.rho, state.rhoV, state.Etot, state.V, state.P);
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
    using info_type    = core_type::Info_t;
    using Equations_t  = core_type::Equations_t;

    template<typename T>
    using Rec = core_type::template Rec<T>;

    constexpr static auto Hall = core_type::Hall;

    explicit FVMethodTransformer(level_t& level, auto& model, info_type const& info)
        : level{level}
        , model{model}
        , info{info}
    {
    }


    void operator()(auto& fvm_state, auto& ct_state, auto& state, auto& fluxes, double const newTime)
    {
        TimeSetter setTime{level, model, newTime};

        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, fvm_state, ct_state, state, fluxes))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type finite_volume_method{info, layout};
            finite_volume_method(fvm_state, ct_state, state, fluxes);
        }

        setTime(state.rho, state.V, state.P, state.J);
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

    void operator()(double const newTime, Model::state_type& state, Model::state_type& statenew,
                    auto& fluxes, double const dt)
    {
        TimeSetter setTime{level, model, newTime};

        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, state, statenew, fluxes))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(state, statenew, fluxes, dt);
        }

        setTime(state.rho, state.rhoV, state.Etot);
    }


    level_t& level;
    Model& model;
};
template<typename Model>
FiniteVolumeEulerTransformer(typename Model::amr_types::level_t&, Model&)
    -> FiniteVolumeEulerTransformer<Model>;




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

    void operator()(auto& ct_state, auto& mhd_state)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, ct_state, mhd_state))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type constrained_transport_{info, layout};
            constrained_transport_(ct_state, mhd_state);
        }
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
    void operator()(double const newTime, Model::state_type& res, auto... pairs)
    {
        TimeSetter setTime{level, model, newTime};

        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, res, pairs.state...))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(res, pairs...);
        }

        setTime(res.rho, res.rhoV, res.Etot);
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

    template<template<typename> typename Reconstruction, auto Hall>
    using ConstrainedTransport_t
        = ConstrainedTransportTransformer<GridLayout, Model, Reconstruction, Hall>;

    using RKUtils_t = RKUtilsTransformer<Model>;
};



}; // namespace PHARE::solver

#endif // PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
