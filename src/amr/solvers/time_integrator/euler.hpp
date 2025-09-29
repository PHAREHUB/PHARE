#ifndef PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP
#define PHARE_CORE_NUMERICS_TIME_INTEGRATOR_EULER_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/compute_fluxes.hpp"
#include "amr/solvers/time_integrator/euler_using_computed_flux.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class Euler
{
    using level_t = typename MHDModel::level_t;

    using ComputeFluxes_t          = ComputeFluxes<FVMethodStrategy, MHDModel>;
    using EulerUsingComputedFlux_t = EulerUsingComputedFlux<MHDModel>;

public:
    Euler(PHARE::initializer::PHAREDict const& dict)
        : compute_fluxes_{dict}
    {
    }

    void operator()(MHDModel& model, auto& state, auto& statenew, auto& fluxes, auto& bc,
                    level_t& level, double const currentTime, double const newTime,
                    double dt = std::nan(""))
    {
        if (std::isnan(dt))
            dt = newTime - currentTime;

        compute_fluxes_(model, state, fluxes, bc, level, newTime);

        euler_using_computed_flux_(model, state, statenew, state.E, fluxes, bc, level, newTime, dt);
    }

private:
    ComputeFluxes_t compute_fluxes_;
    EulerUsingComputedFlux_t euler_using_computed_flux_;
};
} // namespace PHARE::solver

#endif
