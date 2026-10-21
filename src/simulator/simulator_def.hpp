#ifndef PHARE_SIMULATOR_SIMULATOR_DEF_HPP
#define PHARE_SIMULATOR_SIMULATOR_DEF_HPP

#include "core/models/options/mhd_options_def.hpp"

#include "initializer/data_provider.hpp"

#include <array>
#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

namespace PHARE
{

class ISimulator
{
public:
    virtual double startTime()   = 0;
    virtual double endTime()     = 0;
    virtual double currentTime() = 0;
    virtual double timeStep()    = 0;

    virtual void initialize() = 0;
    virtual double advance()  = 0;

    virtual std::vector<int> const& domainBox() const    = 0;
    virtual std::vector<double> const& cellWidth() const = 0;
    virtual std::size_t interporder() const              = 0;

    virtual std::string to_str() = 0;

    virtual ~ISimulator() {}

    virtual bool dump_diagnostics(double timestamp, double timestep)
    {
        return false; // overriding optional
    }
    virtual bool dump_restarts(double timestamp, double timestep)
    {
        return false; // overriding optional
    }
};

constexpr std::size_t defaultNbrRefinedParts(std::size_t dim, std::size_t interp);

struct SimOpts
{
    std::size_t dimension     = 1;
    std::size_t interp_order  = 1;
    std::size_t nbRefinedPart = defaultNbrRefinedParts(dimension, interp_order);

    MHDOpts::ReconstructionType reconstruction_type = MHDOpts::ReconstructionType::MHDOff;
    MHDOpts::SlopeLimiterType slope_limiter_type    = MHDOpts::SlopeLimiterType::MHDOff;
    MHDOpts::RiemannSolverType riemann_solver_type  = MHDOpts::RiemannSolverType::MHDOff;
    bool Hall                                       = false;

    // derived — two independent axes, never one as the other's negation
    bool hybrid_enabled = interp_order > 0;
    bool mhd_enabled    = reconstruction_type != MHDOpts::ReconstructionType::MHDOff;

    // all four MHD axes must agree on off-ness; a member fn keeps SimOpts an aggregate
    constexpr bool mhd_axes_consistent() const
    {
        bool const mhd_is_off = !mhd_enabled;
        return mhd_is_off == (slope_limiter_type == MHDOpts::SlopeLimiterType::MHDOff)
               && mhd_is_off == (riemann_solver_type == MHDOpts::RiemannSolverType::MHDOff);
    }

    static SimOpts FROM(initializer::PHAREDict const& dict);
};


SimOpts inline SimOpts::FROM(initializer::PHAREDict const& dict)
{
    SimOpts static const defaults{};

    auto const enum_from = [&](std::string const& key, auto const default_value) {
        using Enum = decltype(default_value);
        return static_cast<Enum>(cppdict::get_value(dict, key, static_cast<int>(default_value)));
    };

    auto const dimension = static_cast<std::size_t>(
        cppdict::get_value(dict, "dimension", static_cast<int>(defaults.dimension)));
    auto const interp_order = static_cast<std::size_t>(
        cppdict::get_value(dict, "interp_order", static_cast<int>(defaults.interp_order)));

    return {
        dimension,
        interp_order,
        static_cast<std::size_t>(
            cppdict::get_value(dict, "refined_particle_nbr",
                               static_cast<int>(defaultNbrRefinedParts(dimension, interp_order)))),
        enum_from("reconstruction", defaults.reconstruction_type),
        enum_from("limiter", defaults.slope_limiter_type),
        enum_from("riemann", defaults.riemann_solver_type),
        cppdict::get_value(dict, "hall", defaults.Hall),
    };
}


constexpr bool operator==(SimOpts const& lhs, SimOpts const& rhs)
{
    return lhs.dimension == rhs.dimension                         //
           and lhs.interp_order == rhs.interp_order               //
           and lhs.nbRefinedPart == rhs.nbRefinedPart             //
           and lhs.reconstruction_type == rhs.reconstruction_type //
           and lhs.slope_limiter_type == rhs.slope_limiter_type   //
           and lhs.riemann_solver_type == rhs.riemann_solver_type //
           and lhs.Hall == rhs.Hall;
}


constexpr auto inline possibleSimulators()
{
    auto constexpr opts = [](std::size_t dim, std::size_t interp, std::size_t nbRefinedPart) {
        return SimOpts{.dimension = dim, .interp_order = interp, .nbRefinedPart = nbRefinedPart};
    };

    return std::array{
        opts(1, 1, 2),
        opts(1, 1, 3),
        opts(1, 2, 2),
        opts(1, 2, 3),
        opts(1, 2, 4),
        opts(1, 3, 2),
        opts(1, 3, 3),
        opts(1, 3, 4),
        opts(1, 3, 5),
        opts(2, 1, 4),
        opts(2, 1, 5),
        opts(2, 1, 8),
        opts(2, 1, 9),
        opts(2, 2, 4),
        opts(2, 2, 5),
        opts(2, 2, 8),
        opts(2, 2, 9),
        opts(2, 2, 16),
        opts(2, 3, 4),
        opts(2, 3, 5),
        opts(2, 3, 8),
        opts(2, 3, 9),
        opts(2, 3, 25),

        // TODO add in the rest of 3d nbrParticles permutations
        opts(3, 1, 6),
        opts(3, 1, 12),
        opts(3, 2, 6),
        opts(3, 2, 12),
        opts(3, 3, 6),
        opts(3, 3, 12),
    };
}


constexpr std::size_t defaultNbrRefinedParts(std::size_t dim, std::size_t interp)
{
    for (auto const& opts : possibleSimulators())
        if (opts.dimension == dim and opts.interp_order == interp)
            return opts.nbRefinedPart;

    assert(false); // is static_assert in constexpr call

    return 0;
}


constexpr auto inline phare_exe_default_simulators()
{
    // feel free to change as you wish
    auto constexpr opts
        = [](std::size_t dim, std::size_t interp, std::size_t nbRefinedPart, auto&&... args) {
              return SimOpts{dim, interp, nbRefinedPart, args...};
          };

    return std::array{
        opts(1, 1, 2), opts(1, 2, 2), opts(1, 3, 2),

        opts(2, 1, 4), opts(2, 2, 4), opts(2, 3, 4),

        opts(3, 1, 6), opts(3, 2, 6), opts(3, 3, 6),
    };
}



} // namespace PHARE

#endif /*PHARE_SIMULATOR_SIMULATOR_DEF_HPP*/
