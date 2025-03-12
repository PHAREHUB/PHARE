#ifndef PHARE_MHD_DEFAULTS_HPP
#define PHARE_MHD_DEFAULTS_HPP

namespace PHARE
{
template<template<typename> typename FVmethod, typename MHDModel>
struct DefaultTimeIntegrator
{
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
