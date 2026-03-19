#ifndef PHARE_CORE_UTILITIES_GHOST_WIDTH_CALCULATOR_HPP
#define PHARE_CORE_UTILITIES_GHOST_WIDTH_CALCULATOR_HPP

#include <cstdint>
#include <type_traits>

namespace PHARE::core
{

// ============================================================================
// Utility: Round up to nearest even number
// ============================================================================

constexpr std::uint32_t roundUpToEven(std::uint32_t n)
{
    return (n % 2 == 0) ? n : n + 1;
}


// ============================================================================
// Model Type Tags
// ============================================================================

struct HybridModelTag
{
};
struct MHDModelTag
{
};
struct MultiModelTag
{
};


// ============================================================================
// Ghost Width Requirement Traits
// ============================================================================

// Base trait for ghost width requirements
template<typename ModelTag, typename Config>
struct GhostWidthRequirement;


// Hybrid model: Based on interpolation order
template<typename Config>
struct GhostWidthRequirement<HybridModelTag, Config>
{
    static constexpr std::uint32_t interp_order = Config::interp_order;

    static constexpr std::uint32_t compute()
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
};


// MHD model: Based on reconstruction stencil
template<typename Config>
struct GhostWidthRequirement<MHDModelTag, Config>
{
    static constexpr std::uint32_t reconstruction_stencil = Config::reconstruction_nghosts;

    static constexpr std::uint32_t compute()
    {
        // Reconstruction stencil + one layer for J on the full ghost box + one more
        // for the J Laplacian used by hyper-resistivity, then round up for magnetic
        // refinement requirements.
        return roundUpToEven(reconstruction_stencil + 2);
    }
};


// Multi-model: Take maximum of all model requirements
template<typename Config>
struct GhostWidthRequirement<MultiModelTag, Config>
{
    static constexpr std::uint32_t hybrid_ghosts
        = GhostWidthRequirement<HybridModelTag, Config>::compute();

    static constexpr std::uint32_t mhd_ghosts
        = GhostWidthRequirement<MHDModelTag, Config>::compute();

    static constexpr std::uint32_t compute()
    {
        // Both are already even, so max is also even
        return (hybrid_ghosts > mhd_ghosts) ? hybrid_ghosts : mhd_ghosts;
    }
};


// ============================================================================
// Configuration Types
// ============================================================================

// Hybrid PIC configuration
template<std::uint32_t InterpOrder>
struct HybridConfig
{
    using model_tag                                       = HybridModelTag;
    static constexpr std::uint32_t interp_order           = InterpOrder;
    static constexpr std::uint32_t reconstruction_nghosts = 0; // Not used for Hybrid
};


// MHD configuration
template<std::uint32_t ReconstructionGhosts, std::uint32_t InterpOrder = 1>
struct MHDConfig
{
    using model_tag                             = MHDModelTag;
    static constexpr std::uint32_t interp_order = InterpOrder; // For GridLayout compatibility
    static constexpr std::uint32_t reconstruction_nghosts = ReconstructionGhosts;
};


// Multi-model configuration
template<std::uint32_t InterpOrder, std::uint32_t ReconstructionGhosts>
struct MultiModelConfig
{
    using model_tag                                       = MultiModelTag;
    static constexpr std::uint32_t interp_order           = InterpOrder;
    static constexpr std::uint32_t reconstruction_nghosts = ReconstructionGhosts;
};


// ============================================================================
// Main Calculator
// ============================================================================

template<typename Config>
struct GhostWidthCalculator
{
    using model_tag = typename Config::model_tag;

    static constexpr std::uint32_t compute()
    {
        return GhostWidthRequirement<model_tag, Config>::compute();
    }

    // Convenience aliases
    static constexpr std::uint32_t value         = compute();
    static constexpr std::uint32_t primal_ghosts = compute();
    static constexpr std::uint32_t dual_ghosts   = compute();
};


// ============================================================================
// Convenience Functions for Backward Compatibility
// ============================================================================

// Legacy function matching GridLayout's current interface
template<std::uint32_t interp_order>
constexpr std::uint32_t nbrGhostsFromInterpOrder()
{
    using Config = HybridConfig<interp_order>;
    return GhostWidthCalculator<Config>::value;
}


// For reconstruction-based calculations
template<std::uint32_t reconstruction_nghosts>
constexpr std::uint32_t nbrGhostsFromReconstruction()
{
    using Config = MHDConfig<reconstruction_nghosts>;
    return GhostWidthCalculator<Config>::value;
}


} // namespace PHARE::core

#endif // PHARE_CORE_UTILITIES_GHOST_WIDTH_CALCULATOR_HPP
