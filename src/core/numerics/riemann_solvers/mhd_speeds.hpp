#ifndef PHARE_CORE_NUMERICS_MHD_SPEEDS_HPP
#define PHARE_CORE_NUMERICS_MHD_SPEEDS_HPP

#include <cmath>

namespace PHARE::core
{
auto compute_fast_magnetosonic_(auto gamma, auto const& rho, auto const& B, auto const& BdotB,
                                auto const& P)
{
    // squared sound, Alfven and directional Alfven speeds, without taking their roots
    auto const invRho = 1. / rho;
    auto const c02    = gamma * P * invRho;
    auto const cA2    = BdotB * invRho;
    auto const cAdir2 = B * B * invRho;

    return std::sqrt((c02 + cA2) * 0.5
                     + std::sqrt((c02 + cA2) * (c02 + cA2) - 4.0 * c02 * cAdir2) * 0.5);
}

auto compute_whistler_(auto const& invMeshSize, auto const& rho, auto const& BdotB)
{
    auto const vw = std::sqrt(1 + 0.25 * invMeshSize * invMeshSize) + 0.5 * invMeshSize;
    return std::sqrt(BdotB) * vw / rho;
}
} // namespace PHARE::core


#endif
