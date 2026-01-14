#ifndef PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/base_mhd_timestepper.hpp"
#include "amr/solvers/time_integrator/euler.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class EulerIntegrator : public BaseMHDTimestepper<MHDModel>
{
    using Super = BaseMHDTimestepper<MHDModel>;

public:
    EulerIntegrator(PHARE::initializer::PHAREDict const& dict)
        : Super{dict}
        , euler_{dict}
    {
    }

    // Butcher fluxes are used to accumulate fluxes over multiple stages, the corresponding buffer
    // should only contain the fluxes over one time step. The accumulation over all substeps is
    // delegated to the solver.
    void operator()(MHDModel& model, auto& state, auto& fluxes, auto& bc, auto& level,
                    double const currentTime, double const newTime)
    {
        this->resetButcherFluxes_(model, level);

        euler_(model, state, state, fluxes, bc, level, currentTime, newTime);

        this->accumulateButcherFluxes_(model, state.E, fluxes, level);
    }

    using Super::allocate;
    using Super::exposeFluxes;
    using Super::fillMessengerInfo;
    using Super::getCompileTimeResourcesViewList;
    using Super::registerResources;

private:
    Euler<FVMethodStrategy, MHDModel> euler_;
};
} // namespace PHARE::solver

#endif
