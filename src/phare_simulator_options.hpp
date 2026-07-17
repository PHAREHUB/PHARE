#ifndef PHARE_SIMULATOR_OPTIONS_HPP
#define PHARE_SIMULATOR_OPTIONS_HPP

#include "core/def/phare_config.hpp"
#include "core/data/particles/particle_array.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstddef>

namespace PHARE
{

namespace MHDOpts
{
    enum class TimeIntegratorType : uint8_t { Default, Euler, TVDRK2, TVDRK3, SSPRK4_5, count };
    enum class ReconstructionType : uint8_t { Default, Constant, Linear, WENO3, WENOZ, MP5, count };
    enum class SlopeLimiterType : uint8_t { None, VanLeer, MinMod, count };
    enum class RiemannSolverType : uint8_t { Default, Rusanov, HLL, HLLD, count };
}; // namespace MHDOpts


struct SimOpts
{
    std::size_t dimension     = 1;
    std::size_t interp_order  = 1;
    std::size_t nbRefinedPart = core::defaultNbrRefinedParts(dimension, interp_order);

    core::LayoutMode layout_mode = core::LayoutMode::AoSMapped;
    AllocatorMode alloc_mode     = AllocatorMode::CPU;

    MHDOpts::TimeIntegratorType time_integrator_type = MHDOpts::TimeIntegratorType::Default;
    MHDOpts::ReconstructionType reconstruction_type  = MHDOpts::ReconstructionType::Default;
    MHDOpts::SlopeLimiterType slope_limiter_type     = MHDOpts::SlopeLimiterType::None;
    MHDOpts::RiemannSolverType riemann_solver_type   = MHDOpts::RiemannSolverType::Default;
    bool Hall                                        = false;
    bool Resistivity                                 = false;
    bool HyperResistivity                            = false;

    template<auto>
    struct Particles;
};

template<auto sopts>
struct SimOpts::Particles
{
    constexpr static core::ParticleArrayOptions opts{sopts.dimension, sopts.layout_mode,
                                                     core::StorageMode::VECTOR, sopts.alloc_mode};
    using value_type = core::ParticleArray<opts>;
};

} // namespace PHARE

#endif // PHARE_SIMULATOR_OPTIONS_HPP
