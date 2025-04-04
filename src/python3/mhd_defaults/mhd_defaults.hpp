#ifndef PHARE_MHD_DEFAULTS_HPP
#define PHARE_MHD_DEFAULTS_HPP

#include "initializer/data_provider.hpp"

namespace PHARE
{
template<template<typename> typename FVmethod, typename MHDModel>
struct DefaultTimeIntegrator
{
    DefaultTimeIntegrator(PHARE::initializer::PHAREDict const& /*dict*/) {}

    void operator()(MHDModel& /*model*/, MHDModel::state_type& /*state*/, auto& /*fluxes*/,
                    auto& /*fromCoarser*/, auto& /*level*/, double const /*currentTime*/,
                    double const /*newTime*/)
    {
    }

    void registerResources(MHDModel& /*model*/) {}

    void allocate(MHDModel& /*model*/, auto& /*patch*/, double const /*allocateTime*/) const {}

    void fillMessengerInfo(auto& /*info*/) const {}
};

template<typename GridLayout, typename SlopeLimiter>
struct DefaultReconstruction
{
};

template<typename GridLayout, bool HallFlag>
struct DefaultRiemannSolver
{
};

template<bool Hall, bool Resistivity, bool HyperResistivity>
struct DefaultEquations
{
};
} // namespace PHARE

#endif
