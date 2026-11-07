#ifndef PHARE_CORE_NUMERICS_MOMENTS_LEVEL_GHOST_OWNERSHIP_HPP
#define PHARE_CORE_NUMERICS_MOMENTS_LEVEL_GHOST_OWNERSHIP_HPP

#include "core/utilities/box/box.hpp"

#include <vector>
#include <cstddef>

namespace PHARE::core
{
// Level ghost particles are refined from the coarser level into the whole particle ghost box of
// every patch, so the particles of a level ghost cell within the ghost boxes of several same-level
// patches have a copy on each of them. All copies are pushed, but only one patch deposits its
// copy, so that the border sums add it once: the first of those patches in lexicographic order of
// (lower, upper) corners. This order is translation invariant, so periodic images agree on it,
// and it does not depend on the MPI decomposition.
template<std::size_t dim>
struct LevelGhostOwnership
{
    using Box_t = Box<int, dim>;

    Box_t ghostBox;           // patch box grown by the particle ghost width
    std::vector<Box_t> ceded; // parts of ghostBox in the ghost boxes of patches before this one

    bool owns(auto const& cell) const
    {
        return isIn(Point{cell}, ghostBox) and !isIn(Point{cell}, ceded);
    }
};

// Two ghost boxes of width w overlap iff the patch boxes are at most 2w - 1 cells apart along every
// direction, i.e. iff one patch box grown by 2w overlaps the other: the width of the neighbour
// search that finds every patch whose ghost box overlaps that of a given patch.
constexpr int levelGhostNeighbourWidth(int const particleGhostWidth)
{
    return 2 * particleGhostWidth;
}

// `neighbours`: same-level boxes, periodic images included, in the index space of `patch`; at least
// all those whose ghost box overlaps that of `patch`, others are ignored.
template<std::size_t dim>
auto makeLevelGhostOwnership(Box<int, dim> const& patch, auto const& neighbours,
                             int const particleGhostWidth)
{
    LevelGhostOwnership<dim> ownership{grow(patch, particleGhostWidth), {}};
    for (auto const& neighbour : neighbours)
        if (lexicographicLess(neighbour, patch))
            if (auto const overlap = grow(neighbour, particleGhostWidth) * ownership.ghostBox)
                ownership.ceded.push_back(*overlap);
    return ownership;
}

} // namespace PHARE::core

#endif
