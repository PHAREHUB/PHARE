#ifndef PHARE_MHD_RESOLVER_HPP
#define PHARE_MHD_RESOLVER_HPP

#include "core/numerics/godunov_fluxes/godunov_fluxes.hpp"
#include "core/numerics/reconstructions/reconstruction_nghosts.hpp"
#include "phare_simulator_options.hpp"

#include "amr/solvers/time_integrator/time_integrator.hpp"

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

template<MHDOpts::ReconstructionType T>
struct ReconstructionSelector;

template<MHDOpts::ReconstructionType R, MHDOpts::SlopeLimiterType S>
struct SlopeLimiterSelector;

template<MHDOpts::RiemannSolverType T>
struct RiemannSolverSelector;

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::Default>
{
    static constexpr std::uint32_t nghosts
        = MHDOpts::reconstruction_nghosts_v<MHDOpts::ReconstructionType::Default>;
    template<typename GridLayout, typename SlopeLimiter>
    using type = DefaultReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::Constant>
{
    static constexpr std::uint32_t nghosts
        = MHDOpts::reconstruction_nghosts_v<MHDOpts::ReconstructionType::Constant>;
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::ConstantReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::Linear>
{
    static constexpr std::uint32_t nghosts
        = MHDOpts::reconstruction_nghosts_v<MHDOpts::ReconstructionType::Linear>;
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::LinearReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::WENO3>
{
    static constexpr std::uint32_t nghosts
        = MHDOpts::reconstruction_nghosts_v<MHDOpts::ReconstructionType::WENO3>;
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::WENO3Reconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::WENOZ>
{
    static constexpr std::uint32_t nghosts
        = MHDOpts::reconstruction_nghosts_v<MHDOpts::ReconstructionType::WENOZ>;
    template<typename GridLayout, typename SlopeLimiter>
    using type = core::WENOZReconstruction<GridLayout, SlopeLimiter>;
};

template<>
struct ReconstructionSelector<MHDOpts::ReconstructionType::MP5>
{
    static constexpr std::uint32_t nghosts
        = MHDOpts::reconstruction_nghosts_v<MHDOpts::ReconstructionType::MP5>;
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

// Avoid instantiating TimeIntegrator for non-MHD opts (Default reconstruction has stub types
// that would fail inside Godunov). Partial specialization ensures only the selected branch compiles.
template<bool IsMHD, typename FVMethod, typename MHDModel>
struct MHDTimestepperSelector
{
    using type = solver::TimeIntegrator<FVMethod, MHDModel>;
};

template<typename FVMethod, typename MHDModel>
struct MHDTimestepperSelector<false, FVMethod, MHDModel>
{
    using type = DefaultTimeIntegrator<FVMethod, MHDModel>;
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
    using VecField   = MHDModel::vecfield_type;

    using Equations_t = core::MHDEquations<Hall>;

    using RiemannSolver_t = RiemannSolver<Hall>;

    template<typename Layout>
    using Reconstruction_t = Reconstruction<Layout, SlopeLimiter>;

    using FVMethodStrategy
        = core::Godunov<GridLayout, Reconstruction_t, RiemannSolver_t, Equations_t>;

    static constexpr bool is_mhd
        = (opts.reconstruction_type != MHDOpts::ReconstructionType::Default);

    using MHDTimeStepper_t =
        MHDTimestepperSelector<is_mhd, FVMethodStrategy, MHDModel>::type;
};
} // namespace PHARE

#endif // PHARE_MHD_RESOLVER_HPP
