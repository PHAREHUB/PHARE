#ifndef PHARE_CORE_NUMERICS_ION_UPDATER_DEF_HPP
#define PHARE_CORE_NUMERICS_ION_UPDATER_DEF_HPP


#include "core/utilities/box/box.hpp"
#include "core/utilities/point/point.hpp"
#include "core/utilities/range/range.hpp"

#include <functional>

namespace PHARE::core
{


enum class UpdaterMode { domain_only = 1, all = 2 };

// this is to detach how we partition particles from the updater directly
template<typename GridLayout, typename ParticleArray_t>
struct UpdaterSelectionBoxing
{
    auto constexpr static partGhostWidth = GridLayout::nbrParticleGhosts();
    using GridLayout_t                   = GridLayout;
    using Box_t                          = GridLayout_t::AMRBox_t;
    using ParticleRange                  = IndexRange<ParticleArray_t>;
    using Selector_t                     = std::function<ParticleRange(ParticleRange&)>;


    GridLayout_t const layout;
    std::vector<Box_t> const patchGhostBox;
    Box_t const domainBox = layout.AMRBox();
    Box_t const ghostBox  = grow(domainBox, partGhostWidth);

    bool isInGhostBox(auto& particle) const { return isIn(particle, ghostBox); };
    bool isInDomainBox(auto& particle) const { return isIn(particle, domainBox); };
    bool isInPatchGhostBox(auto& particle) const { return isIn(particle, patchGhostBox); };
    bool isInLevelGhostBox(auto& particle) const { return !isInNonLevelGhostBox(particle); };
    bool isInNonLevelGhostBox(auto& particle) const
    {
        return isInDomainBox(particle) or isIn(particle, patchGhostBox);
    };

    Selector_t const noop = [](auto& particleRange) { return particleRange; };

    // lambda copy captures to detach from above references in case of class copy construct
    Selector_t const inDomainBox = [domainBox = domainBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& particle) { return core::isIn(particle, domainBox); });
    };

    Selector_t const inGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& particle) { return isIn(particle, ghostBox); });
    };

    Selector_t const inNonLevelGhostBox
        = [domainBox = domainBox, patchGhostBox = patchGhostBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& particle) {
                  return isIn(particle, domainBox) or isIn(particle, patchGhostBox);
              });
          };

    Selector_t const inGhostLayer
        = [ghostBox = ghostBox, domainBox = domainBox](auto& particleRange) {
              return particleRange.array().partition(particleRange, [&](auto const& particle) {
                  return isIn(particle, ghostBox) and !isIn(particle, domainBox);
              });
          };

    Selector_t const outsideGhostBox = [ghostBox = ghostBox](auto& particleRange) {
        return particleRange.array().partition(
            particleRange, [&](auto const& particle) { return !isIn(particle, ghostBox); });
    };
};

} // namespace PHARE::core


#endif // PHARE_CORE_NUMERICS_ION_UPDATER_DEF_HPP
