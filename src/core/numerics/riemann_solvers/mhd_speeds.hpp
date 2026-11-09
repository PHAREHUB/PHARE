#ifndef PHARE_CORE_NUMERICS_MHD_SPEEDS_HPP
#define PHARE_CORE_NUMERICS_MHD_SPEEDS_HPP

#include <algorithm>
#include <cmath>
#include <numbers>

namespace PHARE::core
{
auto compute_fast_magnetosonic_(auto gamma, auto const& rho, auto const& B, auto const& BdotB,
                                auto const& P)
{
    auto const Sound     = std::sqrt((gamma * P) / rho);
    auto const AlfvenDir = std::sqrt(B * B / rho); // directionnal alfven
    auto const Alfven    = std::sqrt(BdotB / rho);

    auto const c02    = Sound * Sound;
    auto const cA2    = Alfven * Alfven;
    auto const cAdir2 = AlfvenDir * AlfvenDir;

    return std::sqrt((c02 + cA2) * 0.5
                     + std::sqrt(std::max(0.0, (c02 + cA2) * (c02 + cA2) - 4.0 * c02 * cAdir2))
                           * 0.5);
}

auto compute_whistler_(auto const& invMeshSize, auto const& rho, auto const& BdotB)
{
    auto const sqrtRho = std::sqrt(rho);
    auto const vA      = std::sqrt(BdotB) / sqrtRho;
    auto const x       = std::numbers::pi * invMeshSize / (2.0 * sqrtRho);
    return vA * (std::sqrt(1.0 + x * x) + x);
}
} // namespace PHARE::core


#endif
