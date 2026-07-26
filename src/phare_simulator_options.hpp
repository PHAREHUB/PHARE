#ifndef PHARE_SIMULATOR_OPTIONS_HPP
#define PHARE_SIMULATOR_OPTIONS_HPP

#include "core/utilities/meta/meta_utilities.hpp"
#include "core/utilities/meta/enum.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

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
        static constexpr std::array names{
            enumEntry("euler", MHDOpts::TimeIntegratorType::Euler),
            enumEntry("tvdrk2", MHDOpts::TimeIntegratorType::TVDRK2),
            enumEntry("tvdrk3", MHDOpts::TimeIntegratorType::TVDRK3),
            enumEntry("ssprk4_5", MHDOpts::TimeIntegratorType::SSPRK4_5),
        };
    };

    template<>
    struct EnumTraits<MHDOpts::ReconstructionType>
    {
        static constexpr std::string_view label = "reconstruction type";
        static constexpr std::array names{
            enumEntry("constant", MHDOpts::ReconstructionType::Constant),
            enumEntry("linear", MHDOpts::ReconstructionType::Linear),
            enumEntry("weno3", MHDOpts::ReconstructionType::WENO3),
            enumEntry("wenoz", MHDOpts::ReconstructionType::WENOZ),
            enumEntry("mp5", MHDOpts::ReconstructionType::MP5),
        };
    };

    template<>
    struct EnumTraits<MHDOpts::SlopeLimiterType>
    {
        static constexpr std::string_view label = "slope limiter type";
        static constexpr std::array names{
            enumEntry("none", MHDOpts::SlopeLimiterType::None),
            enumEntry("vanleer", MHDOpts::SlopeLimiterType::VanLeer),
            enumEntry("minmod", MHDOpts::SlopeLimiterType::MinMod),
        };
    };

    template<>
    struct EnumTraits<MHDOpts::RiemannSolverType>
    {
        static constexpr std::string_view label = "riemann solver type";
        static constexpr std::array names{
            enumEntry("rusanov", MHDOpts::RiemannSolverType::Rusanov),
            enumEntry("hll", MHDOpts::RiemannSolverType::HLL),
            enumEntry("hlld", MHDOpts::RiemannSolverType::HLLD),
        };
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
