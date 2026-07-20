#ifndef PHARE_OPTIONS_MHD_OPTIONS_HPP
#define PHARE_OPTIONS_MHD_OPTIONS_HPP

#include "core/models/quantities/mhd_quantities.hpp"
#include "core/utilities/ghost_width_calculator.hpp"
#include "core/data/grid/impl/yee/gridlayout_mhd_yee.hpp"
#include "core/numerics/reconstructions/reconstruction_nghosts.hpp"


namespace PHARE
{

template<auto opts_>
struct MHDFieldOptions
{
    using Quantity = core::MHDQuantity;
    using Scalar   = Quantity::Scalar;
    using Vector   = Quantity::Vector;
    using Tensor   = Quantity::Tensor;

    auto static constexpr opts      = opts_;
    auto static constexpr dimension = opts.dimension;

    auto static constexpr reconstruction_nghosts
        = MHDOpts::reconstruction_nghosts_v<opts.reconstruction_type>;

    auto static constexpr field_ghost_width
        = core::nbrGhostsFromReconstruction<reconstruction_nghosts>();
};


template<MHDFieldOptions opts>
struct MHDOptions
{
    using GridLayoutImpl = core::mhd::GridLayoutImplYee<opts>;
    using FieldOptions   = decltype(opts);

    auto static constexpr field_options     = opts;
    auto static constexpr dimension         = opts.dimension;
    auto static constexpr field_ghost_width = opts.field_ghost_width;
};


} // namespace PHARE

#endif // PHARE_OPTIONS_MHD_OPTIONS_HPP
