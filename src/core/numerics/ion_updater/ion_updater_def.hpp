#ifndef PHARE_ION_UPDATER_DEF_HPP
#define PHARE_ION_UPDATER_DEF_HPP


#include "core/utilities/box/box.hpp"
#include "core/hybrid/hybrid_quantities.hpp"
#include "core/data/tiles/tile_set_overlaps.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <cstdint>

namespace PHARE::core
{
enum class UpdaterMode : std::uint16_t { domain_only = 0, all };

template<typename GridLayout>
struct UpdaterSelectionBoxing;

template<typename Selector_t, typename GridLayout>
struct UpdaterCellMapSelectionBoxing;

template<typename Updater_t, typename GridLayout>
constexpr auto* cellmapselection_boxing_impl()
{
    using Selector_t = Updater_t::Pusher::ParticleSelector;
    return static_cast<UpdaterCellMapSelectionBoxing<Selector_t, GridLayout>*>(0);
}

template<typename Updater_t, typename GridLayout>
constexpr auto* selection_boxing_impl()
{
    using enum LayoutMode;
    if constexpr (Updater_t::ParticleArray_t::layout_mode == AoSMapped)
        return cellmapselection_boxing_impl<Updater_t, GridLayout>();
    else
        return static_cast<UpdaterSelectionBoxing<GridLayout>*>(0);
}



template<typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;

    UpdaterSelectionBoxing(GridLayout_t const& layout_, std::vector<Box_t> const& nonLevelGhostBox_)
        : layout{layout_}
        , nonLevelGhostBox{nonLevelGhostBox_}
    {
    }

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);
};


template<typename Selector_t, typename GridLayout>
struct UpdaterCellMapSelectionBoxing : public UpdaterSelectionBoxing<GridLayout>
{
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;
    using Super                          = UpdaterSelectionBoxing<GridLayout>;

    UpdaterCellMapSelectionBoxing(GridLayout_t const& layout_,
                                  std::vector<Box_t> const& nonLevelGhostBox_)
        : Super{layout_, nonLevelGhostBox_}
    {
    }

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inDomainBox = [domainBox = Super::domainBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return core::isIn(cell, domainBox); });
    };

    Selector_t const inGhostBox = [ghostBox = Super::ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return isIn(cell, ghostBox); });
    };

    Selector_t const inNonLevelGhostBox
        = [nonLevelGhostBox = Super::nonLevelGhostBox](auto& particleRange) {
              return particleRange.array().partition(
                  particleRange, [&](auto const& cell) { return isIn(cell, nonLevelGhostBox); });
          };

    Selector_t const inGhostLayer
        = [ghostBox = Super::ghostBox, domainBox = Super::domainBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(cell, ghostBox) and !isIn(cell, domainBox);
              });
          };

    Selector_t const outsideGhostBox = [ghostBox = Super::ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return !isIn(cell, ghostBox); });
    };
};



template<typename GridLayout>
struct UpdaterTileSetSelectionBoxing : public UpdaterSelectionBoxing<GridLayout>
{
    auto constexpr static dimension      = GridLayout::dimension;
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;
    using Super                          = UpdaterSelectionBoxing<GridLayout>;
    using BoxSpanSet_t                   = TileBoxSpanSet<dimension>;
    using scalar_t                       = HybridQuantity::Scalar;


    UpdaterTileSetSelectionBoxing(GridLayout_t const& layout_,
                                  std::vector<Box_t> const& nonLevelGhostBox_)
        : Super{layout_, nonLevelGhostBox_}
    {
    }

    static std::size_t constexpr field_quantities_count()
    {
        if constexpr (dimension == 1)
            return 2;
        if constexpr (dimension == 2)
            return 4;
        if constexpr (dimension == 3)
            return 7; // no ALL DUAL!
    }

    auto static field_quantities(auto const& layout)
    {
        // for all primal/dual permutations -- 1 = primal / 0 = dual
        if constexpr (dimension == 1)
            return std::array<BoxSpanSet_t, field_quantities_count()>{
                make_nd_span_set_for_qty(layout, scalar_t::By), // 0
                make_nd_span_set_for_qty(layout, scalar_t::Bx), // 1
            };
        if constexpr (dimension == 2)
            return std::array<BoxSpanSet_t, field_quantities_count()>{
                make_nd_span_set_for_qty(layout, scalar_t::Bz), // 0 0
                make_nd_span_set_for_qty(layout, scalar_t::By), // 0 1
                make_nd_span_set_for_qty(layout, scalar_t::By), // 1 0
                make_nd_span_set_for_qty(layout, scalar_t::Vx), // 1 1
            };
        if constexpr (dimension == 3)
            return std::array<BoxSpanSet_t, field_quantities_count()>{
                make_nd_span_set_for_qty(layout, scalar_t::Bz), // 0 0 1
                make_nd_span_set_for_qty(layout, scalar_t::By), // 0 1 0
                make_nd_span_set_for_qty(layout, scalar_t::Bx), // 1 0 0
                make_nd_span_set_for_qty(layout, scalar_t::Ex), // 0 1 1
                make_nd_span_set_for_qty(layout, scalar_t::Ey), // 1 0 1
                make_nd_span_set_for_qty(layout, scalar_t::Ez), // 1 1 0
                make_nd_span_set_for_qty(layout, scalar_t::Vx), // 1 1 1
            };
    }


    // std::array<BoxSpanSet_t, field_quantities_count()> fieldings =
    // field_quantities(Super::layout); BoxSpanSet_t particles{make_nd_span_set_from(Super::layout,
    // [](auto const& layout) {
    //     return box_from_zero_to_upper_minus_one(
    //         *grow(layout.AMRBox(), GridLayout_t::nbrParticleGhosts()).shape().as_unsigned());
    // })};
};


} // namespace PHARE::core


// ###############################################

#if PHARE_HAVE_MKN_GPU || defined(_MKN_WITH_MKN_KUL_)
#include "mkn/kul/os.hpp"

namespace PHARE::core::detail
{
auto static const timings_dir_str
    = get_env_as("PHARE_ASYNC_TIMES", std::string{".phare/async/multi_updater"});

static bool ion_updater_io_setup = []() {
    mkn::kul::Dir timings{timings_dir_str};
    timings.mk();
    return true;
}();

} // namespace PHARE::core::detail

#endif // PHARE_HAVE_MKN_GPU

// ###############################################


#endif // ION_UPDATER_DEF_HPP
