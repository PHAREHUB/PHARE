#ifndef PHARE_OPTIONS_HYBRID_OPTIONS_HPP
#define PHARE_OPTIONS_HYBRID_OPTIONS_HPP


#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp"

#include <array>


namespace PHARE
{

// SimOpts says if the build is fully kinetic, option structs of tests may not
template<auto opts>
bool constexpr isFullyKinetic()
{
    if constexpr (requires { opts.fully_kinetic; })
        return opts.fully_kinetic;
    else
        return false;
}

template<auto opts>
struct HybridFieldOptions
{
    using Quantity = core::HybridQuantity;
    using Scalar   = Quantity::Scalar;
    using Vector   = Quantity::Vector;
    using Tensor   = Quantity::Tensor;

    auto static constexpr dimension     = opts.dimension;
    auto static constexpr interp_order  = opts.interp_order;
    auto static constexpr fully_kinetic = isFullyKinetic<opts>();
    // fully kinetic: the level ghost band holds every particle whose current or charge reaches
    // the first ghost layer of J and rho outside a fine level, which the current filter reads:
    // one layer, the half support and half a cell (segment end, displacement guard), ceil(3/2 +
    // half support) cells. The fields reach the push and the deposit of its deepest particle,
    // floor(band + 1/2 + half support) nodes, rounded up to even: a ratio 2 fine ghost box must
    // be whole coarse cells (magnetic_refine_patch_strategy.hpp). The single source of both
    // widths for every particle and field datum.
    auto static constexpr field_ghost_width
        = (isFullyKinetic<opts>() ? std::array{4, 6, 6} : std::array{2, 4, 4})[interp_order - 1];
    auto static constexpr particle_ghost_width
        = (isFullyKinetic<opts>() ? std::array{3, 3, 4} : std::array{1, 2, 2})[interp_order - 1];
};

template<HybridFieldOptions opts>
struct HybridOptions
{
    using GridLayoutImpl = core::hybrid::GridLayoutImplYee<opts>;
    using FieldOptions   = decltype(opts);

    auto static constexpr field_options        = opts;
    auto static constexpr dimension            = opts.dimension;
    auto static constexpr interp_order         = opts.interp_order;
    auto static constexpr fully_kinetic        = opts.fully_kinetic; // isFullyKinetic<options>
    auto static constexpr field_ghost_width    = opts.field_ghost_width;
    auto static constexpr particle_ghost_width = opts.particle_ghost_width;
};




} // namespace PHARE

#endif // PHARE_OPTIONS_HYBRID_OPTIONS_HPP
