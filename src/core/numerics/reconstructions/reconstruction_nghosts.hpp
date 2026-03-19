#ifndef PHARE_RECONSTRUCTION_NGHOSTS_HPP
#define PHARE_RECONSTRUCTION_NGHOSTS_HPP

#include "phare_simulator_options.hpp"
#include <cstdint>

namespace PHARE::MHDOpts
{

// Compile-time mapping from ReconstructionType enum to stencil width (nghosts).
// These values correspond to the nghosts member in each reconstruction class.
template<ReconstructionType R>
struct ReconstructionNghosts;

template<>
struct ReconstructionNghosts<ReconstructionType::Default>
{
    // Default is for Hybrid-only simulations, no MHD reconstruction needed.
    // nghosts=0 signals to use Hybrid ghost width instead.
    static constexpr std::uint32_t value = 0;
};

template<>
struct ReconstructionNghosts<ReconstructionType::Constant>
{
    static constexpr std::uint32_t value = 1;
};

template<>
struct ReconstructionNghosts<ReconstructionType::Linear>
{
    static constexpr std::uint32_t value = 2;
};

template<>
struct ReconstructionNghosts<ReconstructionType::WENO3>
{
    static constexpr std::uint32_t value = 2;
};

template<>
struct ReconstructionNghosts<ReconstructionType::WENOZ>
{
    static constexpr std::uint32_t value = 3;
};

template<>
struct ReconstructionNghosts<ReconstructionType::MP5>
{
    static constexpr std::uint32_t value = 3;
};

template<ReconstructionType R>
inline constexpr std::uint32_t reconstruction_nghosts_v = ReconstructionNghosts<R>::value;

} // namespace PHARE::MHDOpts

#endif // PHARE_RECONSTRUCTION_NGHOSTS_HPP
