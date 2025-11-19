#ifndef PHARE_SIMULATOR_OPTIONS_HPP
#define PHARE_SIMULATOR_OPTIONS_HPP

#include "core/utilities/meta/meta_utilities.hpp"
#include "python3/mhd_defaults/default_mhd_time_stepper.hpp"

#include <cstddef>

namespace PHARE
{

template<template<typename> typename MHDTimeStepper = DefaultMHDTimeStepper_t>
struct SimOpts
{
    std::size_t dimension    = 1;
    std::size_t interp_order = 1;

    std::size_t nbRefinedPart = core::defaultNbrRefinedParts(dimension, interp_order);

    template<typename Model>
    using MHDTimeStepper_t = MHDTimeStepper<Model>;
};



} // namespace PHARE

#endif // PHARE_SIMULATOR_OPTIONS_HPP
