#ifndef DEFAULT_MHD_TIME_STEPPER_HPP
#define DEFAULT_MHD_TIME_STEPPER_HPP

#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"
#include "core/numerics/MHD_equations/MHD_equations.hpp"
#include "core/numerics/reconstructions/wenoz.hpp"
#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "python3/MHDResolver.hpp"

namespace PHARE
{
template<typename Model>
struct DefaultMHDTimeStepper
{
    using type
        = MHDResolver<solver::TVDRK3Integrator, core::WENOZReconstruction, void, core::Rusanov,
                      core::MHDEquations, true, false, false>::TimeIntegrator_t<Model>;
};
} // namespace PHARE

#endif // DEFAULT_MHD_TIME_STEPPER_HPP
