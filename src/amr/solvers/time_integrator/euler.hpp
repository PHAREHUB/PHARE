#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class Euler
{
    using Layout        = typename MHDModel::gridlayout_type;
    using Dispatchers_t = Dispatchers<Layout>;

    using FVMethod_t             = Dispatchers_t::template FVMethod_t<FVMethodStrategy>;
    using FiniteVolumeEuler_t    = Dispatchers_t::FiniteVolumeEuler_t;
    using ConstrainedTransport_t = Dispatchers_t::ConstrainedTransport_t;
    using Faraday_t              = Dispatchers_t::Faraday_t;

    using ToPrimitiveConverter_t    = Dispatchers_t::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = Dispatchers_t::ToConservativeConverter_t;

public:
    Euler(PHARE::initializer::PHAREDict const& dict)
        : fvm_{dict["fv_method"]}
        , fv_euler_{dict["fv_euler"]}
        , to_primitive_{dict["to_primitive"]}
        , to_conservative_{dict["to_conservative"]}
    {
    }

    void operator()(auto layouts, auto& state, auto& statenew, auto& fluxes, auto& bc, auto& level,
                    double const currentTime, double const newTime)
    {
        double const dt = newTime - currentTime;

        to_primitive_(layouts, state);

        bc.fillMomentsGhosts(state, level, newTime);

        fvm_(layouts, state, fluxes);

        // unecessary if we decide to store both primitive and conservative variables
        to_conservative_(layouts, state);

        fv_euler_(layouts, state, statenew, dt, fluxes);

        auto& B_fx = fluxes.B_fx;

        bc.fillMagneticFluxGhosts(B_fx, level, newTime);

        if constexpr (MHDModel::dimension >= 2)
        {
            auto& B_fy = fluxes.B_fy;

            bc.fillMagneticFluxGhosts(B_fy, level, newTime);

            if constexpr (MHDModel::dimension == 3)
            {
                auto& B_fz = fluxes.B_fz;

                bc.fillMagneticFluxGhosts(B_fz, level, newTime);
            }
        }

        ct_(layouts, state, fluxes);

        bc.fillElectricGhosts(state.E, level, newTime);

        faraday_(layouts, state, statenew, dt);
    }

private:
    FVMethod_t fvm_;
    FiniteVolumeEuler_t fv_euler_;
    ConstrainedTransport_t ct_{};
    Faraday_t faraday_{};
    ToPrimitiveConverter_t to_primitive_;
    ToConservativeConverter_t to_conservative_;
};
} // namespace PHARE::solver

#endif
