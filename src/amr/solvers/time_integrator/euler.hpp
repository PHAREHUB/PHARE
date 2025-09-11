#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class Euler
{
    using level_t       = typename MHDModel::level_t;
    using Layout        = typename MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<Layout>;

    using Ampere_t               = Dispatchers_t::Ampere_t;
    using FVMethod_t             = Dispatchers_t::template FVMethod_t<FVMethodStrategy>;
    using FiniteVolumeEuler_t    = Dispatchers_t::FiniteVolumeEuler_t;
    using ConstrainedTransport_t = Dispatchers_t::ConstrainedTransport_t;
    using Faraday_t              = Dispatchers_t::Faraday_t;

    using ToPrimitiveConverter_t    = Dispatchers_t::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = Dispatchers_t::ToConservativeConverter_t;

    constexpr static auto Hall             = FVMethod_t::Hall;
    constexpr static auto Resistivity      = FVMethod_t::Resistivity;
    constexpr static auto HyperResistivity = FVMethod_t::HyperResistivity;

public:
    Euler(PHARE::initializer::PHAREDict const& dict)
        : fvm_{dict["fv_method"]}
        , to_primitive_{dict["to_primitive"]}
        , to_conservative_{dict["to_conservative"]}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& statenew, auto& fluxes, auto& bc,
                    level_t& level, double const currentTime, double const newTime)
    {
        double const dt = newTime - currentTime;

        to_primitive_(level, model, newTime, state);

        bc.fillMomentsGhosts(state, level, newTime);

        if constexpr (Hall || Resistivity || HyperResistivity)
        {
            ampere_(level, model, newTime, state);

            bc.fillCurrentGhosts(state.J, level, newTime);
        }

        fvm_(level, model, newTime, state, fluxes);

        // unecessary if we decide to store both primitive and conservative variables
        to_conservative_(level, model, newTime, state);

        fv_euler_(level, model, newTime, state, statenew, fluxes, dt);

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

        faraday_(level, model, state, statenew, dt);
    }

private:
    Ampere_t ampere_;
    FVMethod_t fvm_;
    FiniteVolumeEuler_t fv_euler_;
    ConstrainedTransport_t ct_;
    Faraday_t faraday_;
    ToPrimitiveConverter_t to_primitive_;
    ToConservativeConverter_t to_conservative_;
};
} // namespace PHARE::solver

#endif
