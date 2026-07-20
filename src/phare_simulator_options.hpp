#ifndef PHARE_SIMULATOR_OPTIONS_HPP
#define PHARE_SIMULATOR_OPTIONS_HPP

#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/meta/enum.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <utility>

namespace PHARE
{

// if mhd is off, use default empty objects
namespace MHDOpts
{
    enum class TimeIntegratorType : uint8_t { Default, Euler, TVDRK2, TVDRK3, SSPRK4_5, count };
    enum class ReconstructionType : uint8_t { Default, Constant, Linear, WENO3, WENOZ, MP5, count };
    enum class SlopeLimiterType : uint8_t { None, VanLeer, MinMod, count };
    enum class RiemannSolverType : uint8_t { Default, Rusanov, HLL, HLLD, count };
} // namespace MHDOpts

namespace core
{
    template<>
    struct EnumTraits<MHDOpts::TimeIntegratorType>
    {
        static constexpr std::string_view label = "time integrator type";
        static constexpr std::array<std::pair<std::string_view, MHDOpts::TimeIntegratorType>, 4>
            names{{
                {"euler", MHDOpts::TimeIntegratorType::Euler},
                {"tvdrk2", MHDOpts::TimeIntegratorType::TVDRK2},
                {"tvdrk3", MHDOpts::TimeIntegratorType::TVDRK3},
                {"ssprk4_5", MHDOpts::TimeIntegratorType::SSPRK4_5},
            }};
    };

    template<>
    struct EnumTraits<MHDOpts::ReconstructionType>
    {
        static constexpr std::string_view label = "reconstruction type";
        static constexpr std::array<std::pair<std::string_view, MHDOpts::ReconstructionType>, 5>
            names{{
                {"constant", MHDOpts::ReconstructionType::Constant},
                {"linear", MHDOpts::ReconstructionType::Linear},
                {"weno3", MHDOpts::ReconstructionType::WENO3},
                {"wenoz", MHDOpts::ReconstructionType::WENOZ},
                {"mp5", MHDOpts::ReconstructionType::MP5},
            }};
    };

    template<>
    struct EnumTraits<MHDOpts::SlopeLimiterType>
    {
        static constexpr std::string_view label = "slope limiter type";
        static constexpr std::array<std::pair<std::string_view, MHDOpts::SlopeLimiterType>, 3>
            names{{
                {"none", MHDOpts::SlopeLimiterType::None},
                {"vanleer", MHDOpts::SlopeLimiterType::VanLeer},
                {"minmod", MHDOpts::SlopeLimiterType::MinMod},
            }};
    };

    template<>
    struct EnumTraits<MHDOpts::RiemannSolverType>
    {
        static constexpr std::string_view label = "riemann solver type";
        static constexpr std::array<std::pair<std::string_view, MHDOpts::RiemannSolverType>, 3>
            names{{
                {"rusanov", MHDOpts::RiemannSolverType::Rusanov},
                {"hll", MHDOpts::RiemannSolverType::HLL},
                {"hlld", MHDOpts::RiemannSolverType::HLLD},
            }};
    };
} // namespace core

struct SimOpts
{
    std::size_t dimension    = 1;
    std::size_t interp_order = 1;

    std::size_t nbRefinedPart = core::defaultNbrRefinedParts(dimension, interp_order);

    MHDOpts::ReconstructionType reconstruction_type = MHDOpts::ReconstructionType::Default;
    MHDOpts::SlopeLimiterType slope_limiter_type    = MHDOpts::SlopeLimiterType::None;
    MHDOpts::RiemannSolverType riemann_solver_type  = MHDOpts::RiemannSolverType::Default;
    bool Hall                                       = false;
};


} // namespace PHARE

#endif // PHARE_SIMULATOR_OPTIONS_HPP
