#ifndef PHARE_DEFAULT_MHD_TIME_STEPPER_HPP
#define PHARE_DEFAULT_MHD_TIME_STEPPER_HPP

#include "python3/mhd_defaults/mhd_defaults.hpp"
#include "python3/mhd_resolver.hpp"

namespace PHARE
{
template<typename Model>
struct DefaultMHDTimeStepper
{
    using type
        = MHDResolver<DefaultTimeIntegrator, DefaultReconstruction, void, DefaultRiemannSolver,
                      DefaultEquations, false, false, false>::TimeIntegrator_t<Model>;
};

template<typename Model>
using DefaultMHDTimeStepper_t = typename DefaultMHDTimeStepper<Model>::type;
} // namespace PHARE

#endif // DEFAULT_MHD_TIME_STEPPER_HPP
