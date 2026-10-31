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

    UpdaterSelectionBoxing(GridLayout_t const& layout_, std::vector<Box_t> const& nonLevelGhostBox_,
                           std::vector<Box_t> const& nonOwnedBoxes_ = {})
        : layout{layout_}
        , nonLevelGhostBox{nonLevelGhostBox_}
        , nonOwnedBoxes{nonOwnedBoxes_}
    {
    }

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    // domain cells owned by an overlapping same-level patch, usually empty
    std::vector<Box_t> const nonOwnedBoxes;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);

    static bool isOwned(auto const& cell, Box_t const& domainBox,
                        std::vector<Box_t> const& nonOwnedBoxes)
    {
        return isIn(cell, domainBox)
               and (nonOwnedBoxes.empty() or !isIn(Point{cell}, nonOwnedBoxes));
    }

    // a particle belongs to this patch if and only if its cell is in the domain and owned
    bool isOwned(auto const& cell) const { return isOwned(cell, domainBox, nonOwnedBoxes); }
};


template<typename Selector_t, typename GridLayout>
struct UpdaterCellMapSelectionBoxing : public UpdaterSelectionBoxing<GridLayout>
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout_t::dimension>;
    using Super                          = UpdaterSelectionBoxing<GridLayout>;

    UpdaterCellMapSelectionBoxing(GridLayout_t const& layout_,
                                  std::vector<Box_t> const& nonLevelGhostBox_,
                                  std::vector<Box_t> const& nonOwnedBoxes_ = {})
        : Super{layout_, nonLevelGhostBox_, nonOwnedBoxes_}
    {
    }

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inOwnedDomain = [domainBox     = Super::domainBox,
                                      nonOwnedBoxes = Super::nonOwnedBoxes](auto& particleRange) {
        return particleRange.array().partition(particleRange, [&](auto const& cell) {
            return Super::isOwned(cell, domainBox, nonOwnedBoxes);
        });
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
