#ifndef PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "amr/solvers/time_integrator/euler.hpp"

namespace PHARE::solver
{
template<template<typename> typename FVMethodStrategy, typename MHDModel>
class EulerIntegrator
{
public:
    EulerIntegrator(PHARE::initializer::PHAREDict const& dict)
        : euler_{dict}
    {
    }

    void operator()(auto layouts, auto& state, auto& fluxes, auto& bc, auto& level,
                    double const currentTime, double const newTime)
    {
        euler_(layouts, state, state, fluxes, bc, level, currentTime, newTime);
    }

private:
    Euler<FVMethodStrategy, MHDModel> euler_;
};
} // namespace PHARE::solver

#endif
