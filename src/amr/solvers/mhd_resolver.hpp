#ifndef PHARE_AMR_SOLVERS_MHD_RESOLVER_HPP
#define PHARE_AMR_SOLVERS_MHD_RESOLVER_HPP

#include "amr/solvers/time_integrator/time_integrator.hpp"

#include "core/models/options/mhd_options_def.hpp"

#include "core/numerics/MHD_equations/MHD_equations.hpp"
#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"
#include "core/numerics/reconstructions/constant.hpp"
#include "core/numerics/reconstructions/linear.hpp"
#include "core/numerics/reconstructions/mp5.hpp"
#include "core/numerics/reconstructions/weno3.hpp"
#include "core/numerics/reconstructions/wenoz.hpp"
#include "core/numerics/riemann_solvers/hll.hpp"
#include "core/numerics/riemann_solvers/hlld.hpp"
#include "core/numerics/riemann_solvers/rusanov.hpp"
#include "core/numerics/slope_limiters/min_mod.hpp"
#include "core/numerics/slope_limiters/van_leer.hpp"


namespace PHARE::solver
{

// Selectors
template<MHDOpts::ReconstructionType T>
struct ReconstructionSelector;

template<MHDOpts::ReconstructionType R, MHDOpts::SlopeLimiterType S>
struct SlopeLimiterSelector;

template<MHDOpts::RiemannSolverType T>
struct RiemannSolverSelector;

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

// SlopeLimiterSelector is only declared above, never defined: every (reconstruction, limiter) pair
// we support must be listed explicitly below, and any pair that is not listed fails to compile
// rather than silently resolving to something. That is how a half-configured MHD build is caught --
// e.g. reconstruction set but limiter left at MHDOff has no specialization, so it does not build.
// Only Linear actually consults a limiter; the others still need an entry for None, resolving to
// void, to say "this combination is valid, the limiter is simply unused".
template<>
struct SlopeLimiterSelector<MHDOpts::ReconstructionType::Constant, MHDOpts::SlopeLimiterType::None>
{
    using type = void;
};

template<>
struct SlopeLimiterSelector<MHDOpts::ReconstructionType::WENO3, MHDOpts::SlopeLimiterType::None>
{
    using type = void;
};

template<>
struct SlopeLimiterSelector<MHDOpts::ReconstructionType::WENOZ, MHDOpts::SlopeLimiterType::None>
{
    using type = void;
};

template<>
struct SlopeLimiterSelector<MHDOpts::ReconstructionType::MP5, MHDOpts::SlopeLimiterType::None>
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

    static constexpr bool Hall = opts.Hall;

    using SlopeLimiter
        = SlopeLimiterSelector<opts.reconstruction_type, opts.slope_limiter_type>::type;

    template<bool HallFlag>
    using RiemannSolver = RiemannSolverSelector<opts.riemann_solver_type>::template type<HallFlag>;

    template<typename Layout, typename Limiter>
    using Reconstruction
        = ReconstructionSelector<opts.reconstruction_type>::template type<Layout, Limiter>;

    // Resolution

    using GridLayout = MHDModel::gridlayout_type;

    using Equations_t = core::MHDEquations<Hall>;

    using RiemannSolver_t = RiemannSolver<Hall>;

    template<typename Layout>
    using Reconstruction_t = Reconstruction<Layout, SlopeLimiter>;

    using FVMethodStrategy
        = core::Godunov<GridLayout, Reconstruction_t, RiemannSolver_t, Equations_t>;

    static constexpr bool is_mhd
        = (opts.reconstruction_type != MHDOpts::ReconstructionType::MHDOff);

    using MHDTimeStepper_t = solver::TimeIntegrator<FVMethodStrategy, MHDModel>;
};

} // namespace PHARE::solver

#endif // PHARE_AMR_SOLVERS_MHD_RESOLVER_HPP
