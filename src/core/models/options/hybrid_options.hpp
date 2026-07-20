#ifndef PHARE_OPTIONS_HYBRID_OPTIONS_HPP
#define PHARE_OPTIONS_HYBRID_OPTIONS_HPP


#include "core/models/quantities/hybrid_quantities.hpp"
#include "core/data/grid/impl/yee/gridlayout_hybrid_yee.hpp"

#include <array>


namespace PHARE
{

template<auto opts>
struct HybridFieldOptions
{
    using Quantity = core::HybridQuantity;
    using Scalar   = Quantity::Scalar;
    using Vector   = Quantity::Vector;
    using Tensor   = Quantity::Tensor;

    auto static constexpr dimension            = opts.dimension;
    auto static constexpr interp_order         = opts.interp_order;
    auto static constexpr field_ghost_width    = std::array{2, 4, 4}[interp_order - 1];
    auto static constexpr particle_ghost_width = std::array{1, 2, 2}[interp_order - 1];
};

template<HybridFieldOptions opts>
struct HybridOptions
{
    using GridLayoutImpl = core::hybrid::GridLayoutImplYee<opts>;
    using FieldOptions   = decltype(opts);

    auto static constexpr field_options        = opts;
    auto static constexpr dimension            = opts.dimension;
    auto static constexpr interp_order         = opts.interp_order;
    auto static constexpr field_ghost_width    = opts.field_ghost_width;
    auto static constexpr particle_ghost_width = opts.particle_ghost_width;
};




} // namespace PHARE

#endif // PHARE_OPTIONS_HYBRID_OPTIONS_HPP
