// Level ghost ownership is a partition of the level ghost region.
//
// Every patch holds a copy of the level ghost particles in its particle ghost box (width w) and
// deposits those in the cells it owns. A level ghost cell (outside every same-level box) must be
// owned by exactly one patch, among those whose ghost box contains it, or its particles are
// deposited twice or never once the border sums add up the patches' deposits.
//
// The neighbours each patch sees are those SAMRAI's same-level connector of width
// levelGhostNeighbourWidth(w) returns: boxes overlapping the patch box grown by that width,
// periodic images included. The check is a brute force over all patches and their images.

#include "core/numerics/moments/level_ghost_ownership.hpp"
#include "core/utilities/box/box_spans.hpp"

#include "gtest/gtest.h"

#include <array>
#include <vector>
#include <string>
#include <cstddef>
#include <sstream>

using namespace PHARE::core;

namespace
{
template<std::size_t dim>
using Box_t = Box<int, dim>;

template<std::size_t dim>
struct Image
{
    std::size_t real;      // index of the patch it is an image of
    Point<int, dim> shift; // periodic shift, 0 for the patch itself
    Box_t<dim> box;
};

// `period[d] == 0`: not periodic along d
template<std::size_t dim>
auto makeImages(std::vector<Box_t<dim>> const& patches, std::array<int, dim> const& period)
{
    std::vector<Image<dim>> images;
    std::size_t nShifts = 1;
    for (std::size_t d = 0; d < dim; ++d)
        nShifts *= 3;
    for (std::size_t ip = 0; ip < patches.size(); ++ip)
        for (std::size_t i = 0; i < nShifts; ++i)
        {
            Point<int, dim> shift;
            bool valid = true;
            for (std::size_t d = 0, k = i; d < dim; ++d, k /= 3)
            {
                int const s = static_cast<int>(k % 3) - 1;
                valid       = valid and (s == 0 or period[d] != 0);
                shift[d]    = s * period[d];
            }
            if (valid)
                images.push_back({ip, shift, PHARE::core::shift(patches[ip], shift)});
        }
    return images;
}

// returns a description of the first cell not owned exactly once, empty if none
template<std::size_t dim>
std::string checkPartition(std::vector<Box_t<dim>> const& patches,
                           std::array<int, dim> const& period, int const w,
                           int const neighbourWidth)
{
    auto const images = makeImages(patches, period);

    std::vector<Box_t<dim>> allBoxes;
    for (auto const& image : images)
        allBoxes.push_back(image.box);

    // what each patch builds from the neighbours a connector of `neighbourWidth` gives it
    std::vector<LevelGhostOwnership<dim>> ownership;
    for (auto const& patch : patches)
    {
        auto const reach = grow(patch, neighbourWidth);
        std::vector<Box_t<dim>> neighbours;
        for (auto const& box : allBoxes)
            if (!(box == patch) and (box * reach))
                neighbours.push_back(box);
        ownership.push_back(makeLevelGhostOwnership(patch, neighbours, w));
    }

    std::string failure;
    for (auto const& patch : patches)
        for_each_point(grow(patch, w), [&](auto const& cell) {
            if (!failure.empty() or isIn(cell, allBoxes))
                return; // not a level ghost cell
            std::size_t owners = 0;
            for (auto const& image : images)
                if (isIn(cell, grow(image.box, w)))
                    owners += ownership[image.real].owns(cell - image.shift);
            if (owners != 1)
            {
                std::stringstream ss;
                ss << "cell " << cell << " has " << owners << " owners";
                failure = ss.str();
            }
        });
    return failure;
}


// 2^dim patches of `size` cells on a grid, with the gap `gaps[d]` between the two along d; the
// upper patches along direction 0 are shifted by `stagger` along the other directions
template<std::size_t dim>
auto makeGrid(std::array<int, dim> const& gaps, int const size, int const stagger)
{
    std::vector<Box_t<dim>> patches;
    for (std::size_t i = 0; i < (std::size_t{1} << dim); ++i)
    {
        Box_t<dim> box;
        for (std::size_t d = 0; d < dim; ++d)
        {
            bool const upper = (i >> d) & 1u;
            box.lower[d]     = upper ? size + gaps[d] : 0;
            if (d > 0 and ((i & 1u) == 1u))
                box.lower[d] += stagger;
            box.upper[d] = box.lower[d] + size - 1;
        }
        patches.push_back(box);
    }
    return patches;
}

template<std::size_t dim>
void checkAllGrids(int const w, int const neighbourWidth, bool const expectPartition)
{
    std::size_t nbrFailures = 0, nbrLayouts = 0;
    int const maxGap    = 2 * w + 1;
    std::size_t nCombos = 1;
    for (std::size_t d = 0; d < dim; ++d)
        nCombos *= maxGap + 1;

    for (int const stagger : {0, 1, 3})
        for (std::size_t i = 0; i < nCombos; ++i)
        {
            std::array<int, dim> gaps;
            for (std::size_t d = 0, k = i; d < dim; ++d, k /= (maxGap + 1))
                gaps[d] = static_cast<int>(k % (maxGap + 1));

            auto const all = makeGrid<dim>(gaps, 5, stagger);

            // all patches, and every pair (along one direction, or diagonal)
            std::vector<std::vector<Box_t<dim>>> layouts{all};
            for (std::size_t a = 0; a < all.size(); ++a)
                for (std::size_t b = a + 1; b < all.size(); ++b)
                    layouts.push_back({all[a], all[b]});

            for (auto const& patches : layouts)
            {
                ++nbrLayouts;
                auto const failure = checkPartition<dim>(patches, {}, w, neighbourWidth);
                if (!failure.empty())
                {
                    ++nbrFailures;
                    if (expectPartition)
                    {
                        ADD_FAILURE() << "w=" << w << " stagger=" << stagger
                                      << " gap[0]=" << gaps[0] << ": " << failure;
                    }
                }
            }
        }
    EXPECT_GT(nbrLayouts, 0u);
    if (!expectPartition)
    {
        EXPECT_GT(nbrFailures, 0u) << "a narrower neighbour search should miss some neighbour";
    }
}

} // namespace


// particle ghost width w = 1 (interp 1), 2 (interp 2 and 3)
class LevelGhostOwnershipTest : public ::testing::TestWithParam<int>
{
};

TEST_P(LevelGhostOwnershipTest, isAPartitionOfTheLevelGhostRegion1D)
{
    auto const w = GetParam();
    checkAllGrids<1>(w, levelGhostNeighbourWidth(w), true);
}

TEST_P(LevelGhostOwnershipTest, isAPartitionOfTheLevelGhostRegion2D)
{
    auto const w = GetParam();
    checkAllGrids<2>(w, levelGhostNeighbourWidth(w), true);
}

TEST_P(LevelGhostOwnershipTest, isAPartitionOfTheLevelGhostRegion3D)
{
    auto const w = GetParam();
    checkAllGrids<3>(w, levelGhostNeighbourWidth(w), true);
}

// the test detects a neighbour search one cell too narrow
TEST_P(LevelGhostOwnershipTest, narrowerNeighbourSearchIsDetected)
{
    auto const w = GetParam();
    checkAllGrids<2>(w, levelGhostNeighbourWidth(w) - 1, false);
}

// patches on both sides of a periodic boundary, gap g across it, gap h inside the domain
TEST_P(LevelGhostOwnershipTest, isAPartitionAcrossPeriodicBoundaries)
{
    auto const w = GetParam();
    for (int g = 0; g <= 2 * w + 1; ++g)
        for (int h = 0; h <= 2 * w + 1; ++h)
        {
            int const L = 5 + h + 5 + g; // period along x
            // 2D, periodic along x only: two patches along x, a third one touching the upper x
            // boundary, above both at gap g along y (its image is diagonal to the first one)
            std::vector<Box_t<2>> patches{Box_t<2>{{0, 0}, {4, 4}},
                                          Box_t<2>{{5 + h, 0}, {9 + h, 4}},
                                          Box_t<2>{{L - 5, 5 + g}, {L - 1, 9 + g}}};
            auto const failure = checkPartition<2>(patches, {L, 0}, w, levelGhostNeighbourWidth(w));
            EXPECT_TRUE(failure.empty())
                << "w=" << w << " g=" << g << " h=" << h << ": " << failure;

            // 1D, fully periodic
            std::vector<Box_t<1>> line{Box_t<1>{{0}, {4}}, Box_t<1>{{5 + h}, {9 + h}}};
            auto const failure1D = checkPartition<1>(line, {L}, w, levelGhostNeighbourWidth(w));
            EXPECT_TRUE(failure1D.empty())
                << "1D w=" << w << " g=" << g << " h=" << h << ": " << failure1D;
        }
}

INSTANTIATE_TEST_SUITE_P(ParticleGhostWidth, LevelGhostOwnershipTest, ::testing::Values(1, 2));


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
