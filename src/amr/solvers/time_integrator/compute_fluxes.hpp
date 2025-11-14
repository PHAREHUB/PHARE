#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_COMPUTE_FLUXES_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_COMPUTE_FLUXES_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class ComputeFluxes
{
    using level_t       = typename MHDModel::level_t;
    using Layout        = typename MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<Layout>;

    using Ampere_t   = Dispatchers_t::Ampere_t;
    using FVMethod_t = Dispatchers_t::template FVMethod_t<FVMethodStrategy>;

    constexpr static auto Hall             = FVMethod_t::Hall;
    constexpr static auto Resistivity      = FVMethod_t::Resistivity;
    constexpr static auto HyperResistivity = FVMethod_t::HyperResistivity;

    using ConstrainedTransport_t
        = Dispatchers_t::template ConstrainedTransport_t<Resistivity, HyperResistivity>;
    using ToPrimitiveConverter_t    = Dispatchers_t::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = Dispatchers_t::ToConservativeConverter_t;


public:
    ComputeFluxes(PHARE::initializer::PHAREDict const& dict)
        : fvm_{dict["fv_method"]}
        , ct_{dict["fv_method"]}
        , to_primitive_{dict["to_primitive"]}
        , to_conservative_{dict["to_conservative"]}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, level_t& level,
                    double const newTime)
    {
        to_primitive_(level, model, newTime, state);

        if constexpr (Hall || Resistivity || HyperResistivity)
        {
            ampere_(level, model, newTime, state);

            bc.fillCurrentGhosts(state.J, level, newTime);
        }

        fvm_(level, model, newTime, state, fluxes);

        // unecessary if we decide to store both primitive and conservative variables
        to_conservative_(level, model, newTime, state);

        bc.fillMagneticFluxesXGhosts(fluxes.B_fx, level, newTime);

        if constexpr (MHDModel::dimension >= 2)
        {
            bc.fillMagneticFluxesYGhosts(fluxes.B_fy, level, newTime);

            if constexpr (MHDModel::dimension == 3)
            {
                bc.fillMagneticFluxesZGhosts(fluxes.B_fz, level, newTime);
            }
        }

        ct_(level, model, state, fluxes);

        bc.fillElectricGhosts(state.E, level, newTime);
    }

private:
    Ampere_t ampere_;
    FVMethod_t fvm_;
    ConstrainedTransport_t ct_;
    ToPrimitiveConverter_t to_primitive_;
    ToConservativeConverter_t to_conservative_;
};
} // namespace PHARE::solver

#endif
