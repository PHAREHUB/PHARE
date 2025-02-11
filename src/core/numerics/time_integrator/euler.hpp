#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/solver_mhd_model_view.hpp"

namespace PHARE::core
{
template<template<typename> typename FVMethod, typename MHDModel>
class Euler
{
    using ModelView_t = solver::MHDModelView<MHDModel>;

    using FVMethod_t             = ModelView_t::template FVMethod_t<FVMethod>;
    using FiniteVolumeEuler_t    = ModelView_t::FiniteVolumeEuler_t;
    using ConstrainedTransport_t = ModelView_t::ConstrainedTransport_t;
    using Faraday_t              = ModelView_t::Faraday_t;

    using ToPrimitiveConverter_t    = ModelView_t::ToPrimitiveConverter_t;
    using ToConservativeConverter_t = ModelView_t::ToConservativeConverter_t;

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

        std::apply(
            [&](auto&&... fluxArgs) {
                fvm_(layouts, state, std::forward<decltype(fluxArgs)>(fluxArgs)...);

                // unecessary if we decide to store both primitive and conservative variables
                to_conservative_(layouts, state);

                fv_euler_(layouts, state, statenew, dt,
                          std::forward<decltype(fluxArgs)>(fluxArgs)...);
            },
            fluxes);

        auto& B_fx = std::get<2>(fluxes);

        bc.fillMagneticFluxGhosts(B_fx, level, newTime);

        if constexpr (MHDModel::dimension >= 2)
        {
            auto& B_fy = std::get<6>(fluxes);

            bc.fillMagneticFluxGhosts(B_fy, level, newTime);

            if constexpr (MHDModel::dimension == 3)
            {
                auto& B_fz = std::get<10>(fluxes);

                bc.fillMagneticFluxGhosts(B_fz, level, newTime);
            }
        }

        std::apply(
            [&](auto&&... fluxArgs) {
                ct_(layouts, state, std::forward<decltype(fluxArgs)>(fluxArgs)...);
            },
            fluxes);

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
} // namespace PHARE::core

#endif
