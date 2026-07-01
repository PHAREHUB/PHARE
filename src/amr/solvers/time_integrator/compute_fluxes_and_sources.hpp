#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_COMPUTE_FLUXES_AND_SOURCES_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_COMPUTE_FLUXES_AND_SOURCES_HPP

#include "initializer/data_provider.hpp"
#include "core/numerics/godunov_fluxes/godunov_utils.hpp"
#include "amr/solvers/solver_mhd_field_evolvers.hpp"

namespace PHARE::solver
{
template<typename FVMethodStrategy, typename MHDModel>
class ComputeFluxesAndSources
{
    using level_t = MHDModel::level_t;
    // using Layout        = MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<MHDModel>;

    using Ampere_t = Dispatchers_t::Ampere_t;

    using FVMethod_t     = Dispatchers_t::template FVMethod_t<FVMethodStrategy>;
    using FVMethodInfo_t = FVMethod_t::info_type;

    constexpr static auto Hall             = FVMethod_t::Hall;
    constexpr static auto Resistivity      = FVMethod_t::Resistivity;
    constexpr static auto HyperResistivity = FVMethod_t::HyperResistivity;

    template<typename T>
    using Rec = FVMethod_t::template Rec<T>;

    using ConstrainedTransport_t
        = Dispatchers_t::template ConstrainedTransport_t<Rec, Hall, Resistivity, HyperResistivity>;
    using ConstrainedTransportInfo_t = ConstrainedTransport_t::info_type;

    using ToPrimitiveConverter_t    = Dispatchers_t::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = Dispatchers_t::ToConservativeConverter_t;
    using ExternalFieldSource_t     = Dispatchers_t::ExternalFieldSource_t;

    using VecField    = MHDModel::vecfield_type;
    using Equations_t = FVMethod_t::Equations_t;


public:
    ComputeFluxesAndSources(PHARE::initializer::PHAREDict const& dict)
        : fVMethodInfo_{FVMethodInfo_t::FROM(dict["fv_method"])}
        , constrainedTransportInfo_{ConstrainedTransportInfo_t::FROM(dict["constrained_transport"])}
        , to_primitive_gamma_{dict["to_primitive"]["heat_capacity_ratio"]}
        , to_conservative_gamma_{dict["to_conservative"]["heat_capacity_ratio"]}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& sources, auto& bc,
                    level_t& level, double const newTime)
    {
        ToPrimitiveConverter_t{level, model}(state, to_primitive_gamma_, newTime);

        if constexpr (Hall || Resistivity || HyperResistivity)
        {
            // NOTE: J is built as curl(B1) only, so the static background contributes curl(B0) = 0
            // to the current. This is exact for a curl-free B0 (the usual background). For a
            // current-carrying B0 (curl(B0) != 0) the Hall/resistive/hyper-resistive terms would
            // miss curl(B0); to support that, thread model.B0 here and compute J = curl(B1 + B0).
            Ampere_t{level, model}(state.B1, state.J);
            TimeSetter{level, model, newTime}(state.B1, state.J);
        }

        FVMethod_t{level, model, fVMethodInfo_}(fvm_, ct_, state, fluxes, newTime);

        // unecessary if we decide to store both primitive and conservative variables
        ToConservativeConverter_t{level, model}(state, to_conservative_gamma_, newTime);

        ConstrainedTransport_t{level, model, constrainedTransportInfo_}(ct_, state);

        // Body sources of the time-dependent B0(x,t) split, computed from THIS stage's B1 into the
        // per-stage `sources` buffer (zero when B0 is static). The integrator accumulates them in a
        // Butcher buffer with the same RK weights as the fluxes, and the Butcher-flux Euler applies
        // them. `state.B1` is unchanged by the prim/cons conversions above.
        ExternalFieldSource_t{level, model}(state, sources);
    }

    void registerResources(MHDModel& model)
    {
        model.resourcesManager->registerResources(fvm_);
        model.resourcesManager->registerResources(ct_);
    }

    void allocate(MHDModel& model, auto& patch, double const allocateTime) const
    {
        model.resourcesManager->allocate(fvm_, patch, allocateTime);
        model.resourcesManager->allocate(ct_, patch, allocateTime);
    }

private:
    FVMethodInfo_t fVMethodInfo_;
    ConstrainedTransportInfo_t constrainedTransportInfo_;

    // Ampere_t ampere_;
    core::GodunovState<VecField, Equations_t> fvm_{};
    core::UpwindConstrainedTransportState<VecField, Hall, Resistivity> ct_{};
    // ToPrimitiveConverter_t to_primitive_;
    // ToConservativeConverter_t to_conservative_;
    double to_primitive_gamma_;
    double to_conservative_gamma_;
};
} // namespace PHARE::solver

#endif
