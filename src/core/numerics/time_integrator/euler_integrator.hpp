#ifndef PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP
#define PHARE_CORE_NUMERICS_EULER_INTEGRATOR_HPP

#include "initializer/data_provider.hpp"
#include "core/numerics/time_integrator/euler.hpp"

namespace PHARE::core
{
template<template<typename> typename FVMethod, typename MHDModel>
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
    Euler<FVMethod, MHDModel> euler_;
};
} // namespace PHARE::core

#endif
