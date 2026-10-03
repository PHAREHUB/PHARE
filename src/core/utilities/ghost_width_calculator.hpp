#ifndef PHARE_CORE_UTILITIES_GHOST_WIDTH_CALCULATOR_HPP
#define PHARE_CORE_UTILITIES_GHOST_WIDTH_CALCULATOR_HPP

#include <cstdint>

namespace PHARE::core
{

// ============================================================================
// Utility: Round up to nearest even number
// ============================================================================

constexpr inline std::uint32_t roundUpToEven(std::uint32_t n)
{
    return (n % 2 == 0) ? n : n + 1;
}


// ============================================================================
// Ghost Width Computation Functions
// ============================================================================

/**
 * @brief Compute ghost width for Hybrid PIC model based on interpolation order.
 *
 * Ghost cells are needed for:
 * - Particle-mesh interpolation: (interp_order + 1) / 2
 * - One extra layer for particles that may leave cells
 * - Rounded to even for Toth & Roe (2002) magnetic refinement formulas
 */
template<std::uint32_t interp_order>
constexpr std::uint32_t nbrGhostsFromInterpOrder()
{
    if constexpr (interp_order == 1)
        return 2;
    else if constexpr (interp_order == 2)
        return 4;
    else if constexpr (interp_order == 3)
        return 4;
    else
        return roundUpToEven((interp_order + 1) / 2 + 1);
}


/**
 * @brief Compute ghost width for MHD model based on reconstruction stencil.
 *
 * Ghost cells are needed for:
 * - Reconstruction stencil width
 * - One layer because ampere computes J on the ghost box shrinked by one, so J is only valid on
 *   ghost_width - 1 layers while the reconstruction reads it out to the stencil width
 * - Rounded to even for Toth & Roe (2002) magnetic refinement formulas
 *
 * There is no hyper-resistivity layer: the dissipative electric field takes its Laplacian on J at
 * its native edge location, and the flux contributions project it to the face, so they stay within
 * the reconstruction reach. See DissipativeElectricField and Godunov::transverse_on_face_.
 */
template<std::uint32_t reconstruction_nghosts>
constexpr std::uint32_t nbrGhostsFromReconstruction()
{
    return roundUpToEven(reconstruction_nghosts + 1);
}


/**
 * @brief For particles, ghost width depends on interpolation order.
 *
 * This is the same as the Hybrid field ghost width.
 */
template<std::uint32_t interp_order>
constexpr std::uint32_t particleGhostWidth()
{
    return nbrGhostsFromInterpOrder<interp_order>();
}


} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_GHOST_WIDTH_CALCULATOR_HPP
