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
            core_type{layout, gamma}(state.rho, state.V, state.B1, state.P, state.rhoV,
                                     state.Etot1);
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
            core_type{layout}(gamma, state.rho, state.rhoV, state.B1, state.Etot1, state.V,
                              state.P);
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


    void operator()(auto& ct_state, auto& dissipative_electric_state, auto& state, auto& fluxes)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, ct_state, dissipative_electric_state, state, fluxes,
                                        model.externalField))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type finite_volume_method{info, layout};
            finite_volume_method(ct_state, dissipative_electric_state, state,
                                 model.externalField.B0, fluxes);
        }
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

    void operator()(auto& dissipative_electric_state, auto& mhd_state)
    {
        auto& rm = *model.resourcesManager;
        for (auto& patch :
             rm.enumerate(level, dissipative_electric_state, mhd_state, model.externalField))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{info, layout}(dissipative_electric_state, mhd_state, model.externalField.B0);
        }
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
        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, ct_state, dissipative_electric_state, mhd_state,
                                        model.externalField))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type constrained_transport_{info, layout};
            constrained_transport_(ct_state, dissipative_electric_state, mhd_state,
                                   model.externalField.B0);
        }
    }


    level_t& level;
    Model& model;
    info_type const info;
};




template<typename Model>
class ExternalFieldSourceTransformer
{
    using GridLayout = Model::gridlayout_type;
    using level_t    = Model::amr_types::level_t;

public:
    explicit ExternalFieldSourceTransformer(level_t& level, Model& model)
        : level{level}
        , model{model}
    {
    }

    void operator()(Model::state_type& state, auto& sources)
    {
        auto& rm = *model.resourcesManager;

        if (not model.externalFieldUpdater->isTimeDependent())
        {
            for (auto& _ : rm.enumerate(level, sources))
                sources.zero();
            return;
        }

        for (auto& patch : rm.enumerate(level, state, sources, model.externalField))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            auto const& dB0dt = model.externalField.dB0dt;

            for (auto const c : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& s        = sources.B1_source(c);
                auto const& dc = dB0dt(c);
                layout.evalOnBox(s, [&](auto&... args) { s(args...) = -dc(args...); });
            }

            layout.evalOnBox(sources.Etot_source, [&](auto&... args) {
                sources.Etot_source(args...) = -dotOnCellCenter_(dB0dt, state.B1, {args...});
            });
        }
    }

private:
    static auto dotOnCellCenter_(auto const& A, auto const& B,
                                 core::MeshIndex<GridLayout::dimension> index)
    {
        using implT = GridLayout::implT;

        auto const ax
            = GridLayout::template project<implT::faceXToCellCenter>(A(core::Component::X), index);
        auto const ay
            = GridLayout::template project<implT::faceYToCellCenter>(A(core::Component::Y), index);
        auto const az
            = GridLayout::template project<implT::faceZToCellCenter>(A(core::Component::Z), index);
        auto const bx
            = GridLayout::template project<implT::faceXToCellCenter>(B(core::Component::X), index);
        auto const by
            = GridLayout::template project<implT::faceYToCellCenter>(B(core::Component::Y), index);
        auto const bz
            = GridLayout::template project<implT::faceZToCellCenter>(B(core::Component::Z), index);

        return ax * bx + ay * by + az * bz;
    }

    level_t& level;
    Model& model;
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

    using ExternalFieldSource_t = ExternalFieldSourceTransformer<Model>;

    template<template<typename> typename Reconstruction, auto Hall>
    using ConstrainedTransport_t
        = ConstrainedTransportTransformer<GridLayout, Model, Reconstruction, Hall>;

    using RKUtils_t = RKUtilsTransformer<Model>;
};



}; // namespace PHARE::solver

#endif // PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
