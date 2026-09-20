#ifndef PHARE_ION_UPDATER_DEF_HPP
#define PHARE_ION_UPDATER_DEF_HPP

#include "core/utilities/box/box.hpp"

#include <cstdint>

namespace PHARE::core
{
enum class UpdaterMode : std::uint16_t { domain_only = 0, all };

template<typename GridLayout>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
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
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
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


} // namespace PHARE::core

#endif // ION_UPDATER_DEF_HPP
