#ifndef PHARE_OPTIONS_HYBRID_OPTIONS_HPP
#define PHARE_OPTIONS_HYBRID_OPTIONS_HPP


#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp"

#include <array>
#include <cstddef>

namespace PHARE
{

template<auto opts>
struct HybridFieldOptions
{
    auto static constexpr dimension    = opts.dimension;
    auto static constexpr interp_order = opts.interp_order;
    std::size_t ghost_width            = std::array{2, 4, 4}[interp_order - 1];

    using Quantity = core::HybridQuantity;
    using Scalar   = Quantity::Scalar;
};

template<HybridFieldOptions opts>
struct HybridOptions
{
    auto static constexpr field_options = opts;
    auto static constexpr dimension     = opts.dimension;
    auto static constexpr interp_order  = opts.interp_order;
    std::size_t ghost_width             = opts.ghost_width;
    using GridLayoutImpl                = core::hybrid::GridLayoutImplYee<opts>;
    using FieldOptions                  = decltype(field_options);
};




} // namespace PHARE

#endif // PHARE_OPTIONS_HYBRID_OPTIONS_HPP
