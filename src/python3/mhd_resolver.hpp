#ifndef PHARE_MHD_RESOLVER_HPP
#define PHARE_MHD_RESOLVER_HPP

#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"
#include "phare_simulator_options.hpp"

#include "amr/solvers/time_integrator/euler_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk2_integrator.hpp"
#include "amr/solvers/time_integrator/tvdrk3_integrator.hpp"
#include "amr/solvers/time_integrator/ssprk4_5_integrator.hpp"

#include "core/numerics/reconstructions/constant.hpp"
#include "core/numerics/reconstructions/linear.hpp"
#include "core/numerics/reconstructions/weno3.hpp"
#include "core/numerics/reconstructions/wenoz.hpp"
#include "core/numerics/reconstructions/mp5.hpp"

#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/numerics/slope_limiters/van_leer.hpp"

#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "core/numerics/riemann_solvers/hll.hpp"
#include "core/numerics/riemann_solvers/hlld.hpp"

#include "core/numerics/MHD_equations/MHD_equations.hpp"
#include "python3/mhd_defaults/mhd_defaults.hpp"

namespace PHARE
{

// Selectors

template<MHDOpts::TimeIntegratorType T, typename MHDModel>
struct TimeIntegratorSelector;

template<MHDOpts::ReconstructionType T>
struct ReconstructionSelector;

template<MHDOpts::ReconstructionType R, MHDOpts::SlopeLimiterType S>
struct SlopeLimiterSelector;

template<MHDOpts::RiemannSolverType T>
struct RiemannSolverSelector;

template<typename MHDModel>
struct TimeIntegratorSelector<MHDOpts::TimeIntegratorType::Default, MHDModel>
{
    template<template<typename> typename FVmethod>
    using type = DefaultTimeIntegrator<FVmethod, MHDModel>;
};

template<typename MHDModel>
struct TimeIntegratorSelector<MHDOpts::TimeIntegratorType::Euler, MHDModel>
{
    template<template<typename> typename FVmethod>
    using type = solver::EulerIntegrator<FVmethod, MHDModel>;
};

template<typename MHDModel>
struct TimeIntegratorSelector<MHDOpts::TimeIntegratorType::TVDRK2, MHDModel>
{
    template<template<typename> typename FVmethod>
    using type = solver::TVDRK2Integrator<FVmethod, MHDModel>;
};

template<typename MHDModel>
struct TimeIntegratorSelector<MHDOpts::TimeIntegratorType::TVDRK3, MHDModel>
{
    template<template<typename> typename FVmethod>
    using type = solver::TVDRK3Integrator<FVmethod, MHDModel>;
};

template<typename MHDModel>
struct TimeIntegratorSelector<MHDOpts::TimeIntegratorType::SSPRK4_5, MHDModel>
{
    template<template<typename> typename FVmethod>
    using type = solver::SSPRK4_5Integrator<FVmethod, MHDModel>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::Default>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = DefaultReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::Constant>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::ConstantReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::Linear>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::LinearReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::WENO3>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::WENO3Reconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::WENOZ>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::WENOZReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::MP5>
{
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::MP5Reconstruction<GridLayout, SlopeLimiter>;
};

template<MHDOpts::ReconstructionType R, MHDOpts::SlopeLimiterType S>
struct SlopeLimiterSelector
{
    using type = void;
};

template<>
struct SlopeLimiterSelector<MHDOpts::ReconstructionType::Linear, MHDOpts::SlopeLimiterType::VanLeer>
{
    using type = core::VanLeerLimiter;
};

template<>
struct SlopeLimiterSelector<MHDOpts::ReconstructionType::Linear, MHDOpts::SlopeLimiterType::MinMod>
{
    using type = core::MinModLimiter;
};

template<>
struct RiemannSolverSelector<MHDOpts::RiemannSolverType::Default>
{
    template<bool Hall>
    using type = DefaultRiemannSolver<Hall>;
};

template<>
struct RiemannSolverSelector<MHDOpts::RiemannSolverType::Rusanov>
{
    template<bool Hall>
    using type = core::Rusanov<Hall>;
};

template<>
struct RiemannSolverSelector<MHDOpts::RiemannSolverType::HLL>
{
    template<bool Hall>
    using type = core::HLL<Hall>;
};

template<>
struct RiemannSolverSelector<MHDOpts::RiemannSolverType::HLLD>
{
    template<bool Hall>
    using type = core::HLLD<Hall>;
};

template<auto opts, typename MHDModel>
struct MHDResolver
{
    // Get the types from opts

    static constexpr bool Hall             = opts.Hall;
    static constexpr bool Resistivity      = opts.Resistivity;
    static constexpr bool HyperResistivity = opts.HyperResistivity;

    using SlopeLimiter
        = SlopeLimiterSelector<opts.reconstruction_type, opts.slope_limiter_type>::type;

    template<bool HallFlag>
    using RiemannSolver = RiemannSolverSelector<opts.riemann_solver_type>::template type<HallFlag>;

    template<typename Layout, typename Limiter>
    using Reconstruction
        = ReconstructionSelector<opts.reconstruction_type>::template type<Layout, Limiter>;

    template<template<typename> typename FVMethod>
    using MHDTimeStepper
        = TimeIntegratorSelector<opts.time_integrator_type, MHDModel>::template type<FVMethod>;

    // Resolution

    using Equations_t = core::MHDEquations<Hall, Resistivity, HyperResistivity>;

    using RiemannSolver_t = RiemannSolver<Hall>;

    template<typename Layout>
    using Reconstruction_t = Reconstruction<Layout, SlopeLimiter>;

    template<typename Layout>
    using FVMethodStrategy
        = core::Godunov<Layout, MHDModel, Reconstruction_t, RiemannSolver_t, Equations_t>;

    using MHDTimeStepper_t = MHDTimeStepper<FVMethodStrategy>;
};
} // namespace PHARE

#endif // PHARE_MHD_RESOLVER_HPP
