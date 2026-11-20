#ifndef PHARE_CORE_NUMERICS_ION_UPDATER_DEF_HPP
#define PHARE_CORE_NUMERICS_ION_UPDATER_DEF_HPP


#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/range/range.hpp"

#include <vector>
#include <functional>


namespace PHARE::core
{
enum class UpdaterMode { in_place = 0, copy };


// this is to detach how we partition particles from the updater directly
template<typename GridLayout, typename ParticleArray>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::options.particle_ghost_width;
    using GridLayout_t                   = GridLayout;
    using Box_t                          = Box<int, GridLayout::dimension>;
    using ParticleRange                  = IndexRange<ParticleArray>;
    using Selector_t                     = std::function<ParticleRange(ParticleRange&)>;

    GridLayout_t const layout;
    std::vector<Box_t> const nonLevelGhostBox;
    // domain cells owned by an overlapping same-level patch, usually empty
    std::vector<Box_t> const nonOwnedBoxes{};
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

    bool isInNonLevelGhostBox(auto const& cell) const
    {
        return isIn(Point{cell}, nonLevelGhostBox);
    }

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inOwnedDomain
        = [domainBox = domainBox, nonOwnedBoxes = nonOwnedBoxes](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isOwned(cell, domainBox, nonOwnedBoxes);
              });
          };

    Selector_t const inGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return isIn(cell, ghostBox); });
    };

    Selector_t const inNonLevelGhostBox
        = [nonLevelGhostBox = nonLevelGhostBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(Point{cell}, nonLevelGhostBox);
              });
          };

    // ghost cells outside the level: cells of a same-level neighbour's domain are excluded, a
    // level ghost particle entering one is that neighbour's (it holds its own copy)
    Selector_t const inLevelGhostLayer
        = [ghostBox = ghostBox, nonLevelGhostBox = nonLevelGhostBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& cell) {
                  return isIn(cell, ghostBox) and !isIn(cell, nonLevelGhostBox);
              });
          };

    Selector_t const outsideGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& cell) { return !isIn(cell, ghostBox); });
    };


} // namespace PHARE::core


#endif // PHARE_CORE_NUMERICS_ION_UPDATER_DEF_HPP
