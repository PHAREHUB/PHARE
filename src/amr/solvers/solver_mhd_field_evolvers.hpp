#ifndef PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
#define PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP


#include "core/numerics/time_integrator_utils.hpp"
#include "core/numerics/finite_volume_euler/finite_volume_euler.hpp"
#include "core/numerics/constrained_transport/upwind_constrained_transport.hpp"
#include "core/numerics/primite_conservative_converter/to_primitive_converter.hpp"
#include "core/numerics/primite_conservative_converter/to_conservative_converter.hpp"
#include "core/data/vecfield/vecfield_component.hpp"
#include "core/utilities/index/index.hpp"

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
            core_type{layout, gamma}(state.rho, state.V, state.B1, state.P, state.rhoV,
                                     state.Etot1);
        }

        setTime(state.rho, state.V, state.P, state.rhoV, state.Etot1);
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
            core_type{layout}(gamma, state.rho, state.rhoV, state.B1, state.Etot1, state.V,
                              state.P);
        }

        setTime(state.rho, state.rhoV, state.Etot1, state.V, state.P);
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

    constexpr static auto Hall             = core_type::Hall;
    constexpr static auto Resistivity      = core_type::Resistivity;
    constexpr static auto HyperResistivity = core_type::HyperResistivity;

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
        for (auto& patch : rm.enumerate(level, fvm_state, ct_state, state, fluxes, model.B0))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type finite_volume_method{info, layout};
            finite_volume_method(fvm_state, ct_state, state, model.B0, fluxes);
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

        setTime(state.rho, state.rhoV, state.Etot1);
    }


    level_t& level;
    Model& model;
};
template<typename Model>
FiniteVolumeEulerTransformer(typename Model::amr_types::level_t&, Model&)
    -> FiniteVolumeEulerTransformer<Model>;




template<typename GridLayout, typename Model, template<typename> typename Reconstruction, auto Hall,
         auto Resistivity, auto HyperResistivity>
class ConstrainedTransportTransformer
{
    using level_t   = Model::amr_types::level_t;
    using core_type = core::UpwindConstrainedTransport<GridLayout, Reconstruction, Hall,
                                                       Resistivity, HyperResistivity>;

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
        for (auto& patch : rm.enumerate(level, ct_state, mhd_state, model.B0))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type constrained_transport_{info, layout};
            constrained_transport_(ct_state, mhd_state, model.B0);
        }
    }


    level_t& level;
    Model& model;
    info_type const info;
};








// Computes the body-source terms of the time-dependent B = B0(x,t) + B1 split (see theory/mhd.md)
// and writes them into a per-stage MHDSources buffer (does NOT touch the state). With B0 frozen
// across the RK stages of a step:
//   B1_source   = -dB0/dt                 (face-centered, added to the B1 induction equation)
//   Etot_source = -dB0/dt . B1_stage      (cell-centered, added to the reduced-energy equation)
// The energy source uses the current stage's B1. When B0 is static (b0TimeDependent_ == false) the
// buffer is zero-filled so the Butcher accumulation/apply is an exact no-op.
template<typename Model>
class ExternalFieldSourceTransformer
{
    using GridLayout                = Model::gridlayout_type;
    using level_t                   = Model::amr_types::level_t;
    static constexpr auto dimension = GridLayout::dimension;

public:
    explicit ExternalFieldSourceTransformer(level_t& level, auto& model)
        : level{level}
        , model{model}
    {
    }

    void operator()(typename Model::state_type& state, auto& sources)
    {
        auto& rm = *model.resourcesManager;

        if (not model.b0TimeDependent_)
        {
            for (auto& patch : rm.enumerate(level, sources))
            {
                auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
                for (auto const& c : {core::Component::X, core::Component::Y, core::Component::Z})
                {
                    auto& s = sources.B1_source(c);
                    layout.evalOnBox(s, [&](auto const&... args) mutable { s(args...) = 0.0; });
                }
                layout.evalOnBox(sources.Etot_source, [&](auto const&... args) mutable {
                    sources.Etot_source(args...) = 0.0;
                });
            }
            return;
        }

        for (auto& patch : rm.enumerate(level, state, sources, model.dB0dt))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            // B1_source = -dB0/dt, co-located with B1 (same face centering).
            for (auto const& c : {core::Component::X, core::Component::Y, core::Component::Z})
            {
                auto& s        = sources.B1_source(c);
                auto const& dc = model.dB0dt(c);
                layout.evalOnBox(s, [&](auto&... args) mutable { s(args...) = -dc(args...); });
            }
            // Etot_source = -(dB0/dt . B1); both face-centered fields projected to the cell center
            // exactly as the B1^2/2 magnetic energy is (ToConservativeConverter), for consistency.
            layout.evalOnBox(sources.Etot_source, [&](auto&... args) mutable {
                energySource_(model.dB0dt, state.B1, sources.Etot_source, {args...});
            });
        }
    }

private:
    template<typename VecField, typename Field>
    static void energySource_(VecField const& dB0dt, VecField const& B1, Field& Etot_source,
                              core::MeshIndex<dimension> index)
    {
        auto const& dx = dB0dt(core::Component::X);
        auto const& dy = dB0dt(core::Component::Y);
        auto const& dz = dB0dt(core::Component::Z);
        auto const& bx = B1(core::Component::X);
        auto const& by = B1(core::Component::Y);
        auto const& bz = B1(core::Component::Z);

        auto const db0x = GridLayout::template project<GridLayout::faceXToCellCenter>(dx, index);
        auto const db0y = GridLayout::template project<GridLayout::faceYToCellCenter>(dy, index);
        auto const db0z = GridLayout::template project<GridLayout::faceZToCellCenter>(dz, index);
        auto const b1x  = GridLayout::template project<GridLayout::faceXToCellCenter>(bx, index);
        auto const b1y  = GridLayout::template project<GridLayout::faceYToCellCenter>(by, index);
        auto const b1z  = GridLayout::template project<GridLayout::faceZToCellCenter>(bz, index);

        Etot_source(index) = -(db0x * b1x + db0y * b1y + db0z * b1z);
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
    void operator()(double const newTime, Model::state_type& res, auto... pairs)
    {
        TimeSetter setTime{level, model, newTime};

        auto& rm = *model.resourcesManager;
        for (auto& patch : rm.enumerate(level, res, pairs.state...))
        {
            auto const layout = amr::layoutFromPatch<GridLayout>(*patch);
            core_type{layout}(res, pairs...);
        }

        setTime(res.rho, res.rhoV, res.Etot1);
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

    using ExternalFieldSource_t = ExternalFieldSourceTransformer<Model>;

    template<template<typename> typename Reconstruction, auto Hall, auto Resistivity,
             auto HyperResistivity>
    using ConstrainedTransport_t
        = ConstrainedTransportTransformer<GridLayout, Model, Reconstruction, Hall, Resistivity,
                                          HyperResistivity>;

    using RKUtils_t = RKUtilsTransformer<Model>;
};



}; // namespace PHARE::solver

#endif // PHARE_AMR_SOLVERS_SOLVER_MHD_FIELD_EVOLVERS_HPP
