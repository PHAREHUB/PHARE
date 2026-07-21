#ifndef PHARE_CORE_NUMERICS_MHD_SPEEDS_HPP
#define PHARE_CORE_NUMERICS_MHD_SPEEDS_HPP

#include <cmath>

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
                     + std::sqrt((c02 + cA2) * (c02 + cA2) - 4.0 * c02 * cAdir2) * 0.5);
}

// Same fast-magnetosonic speed as compute_fast_magnetosonic_ but parameterized by an
// already-known sound speed instead of (gamma, P). Useful where a Thermo object exposes
// soundSpeed() but not gamma — keeps the EOS abstraction intact. @c Bn is the
// boundary-normal magnetic component, @c BdotB the total |B|^2 (both co-located with rho).
auto compute_fast_magnetosonic_from_cs_(auto const& cs, auto const& rho, auto const& Bn,
                                        auto const& BdotB)
{
    auto const c02    = cs * cs;
    auto const cA2    = BdotB / rho;
    auto const cAdir2 = Bn * Bn / rho;

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
