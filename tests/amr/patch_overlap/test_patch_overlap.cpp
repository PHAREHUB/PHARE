// Same-level patches may overlap. These tests cover the pieces that make overlapping
// patches behave as a partition of the level: ownership of shared cells, the region used
// to exchange leaving particles, and the border sum between distinct patches with equal boxes.

#include "phare_core.hpp"
#include "phare_mpi.hpp" // IWYU pragma: keep

#include "core/models/quantities/hybrid_quantities.hpp"
#include "amr/resources_manager/amr_utils.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/data/field/field_variable_fill_pattern.hpp"
#include "amr/data/tensorfield/tensor_field_geometry.hpp"
#include "amr/data/particles/particles_variable_fill_pattern.hpp"

#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/pdat/CellOverlap.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gtest/gtest.h"

#include <set>
#include <limits>
#include <random>
#include <numeric>
#include <algorithm>


using namespace PHARE;

template<std::size_t dim, std::size_t interp>
using GridLayout_t = core::PHARE_Types<SimOpts{dim, interp}>::Hybrid::GridLayout_t;

template<std::size_t dim>
using Box_t = core::Box<int, dim>;

template<std::size_t dim>
using Shift_t = core::Point<int, dim>;


// ---------------------------------------------------------------------------------------------
// particle exchange region: without overlap, equal to the clipped cell overlap region

template<std::size_t dim>
SAMRAI::hier::Box samraiBox(Box_t<dim> const& box)
{
    return amr::samrai_box_from(box);
}

template<std::size_t dim>
std::set<std::array<int, dim>> cellsOf(SAMRAI::hier::BoxContainer const& boxes)
{
    std::set<std::array<int, dim>> cells;
    for (auto const& box : boxes)
        for (auto const& cell : amr::phare_box_from<dim>(box))
            cells.insert(cell.toArray());
    return cells;
}

// reference region: the cell overlap (destination ghost box * source box, minus the
// destination box), grown by the particle ghost width and clipped to the destination box
template<typename GridLayout>
auto clippedRegion(SAMRAI::pdat::CellGeometry const& dst, SAMRAI::pdat::CellGeometry const& src,
                   SAMRAI::hier::Box const& src_mask, SAMRAI::hier::Box const& fill_box,
                   SAMRAI::hier::Transformation const& transformation)
{
    auto constexpr dim = GridLayout::dimension;
    auto const overlap = static_cast<SAMRAI::hier::BoxGeometry const&>(dst).calculateOverlap(
        src, src_mask, fill_box, false, transformation);
    auto const& cellOverlap = dynamic_cast<SAMRAI::pdat::CellOverlap const&>(*overlap);

    SAMRAI::hier::BoxContainer boxes;
    for (auto const& box : cellOverlap.getDestinationBoxContainer())
        if (auto const region
            = core::grow(amr::phare_box_from<dim>(box), GridLayout::options.particle_ghost_width)
              * amr::phare_box_from<dim>(dst.getBox()))
            boxes.pushBack(amr::samrai_box_from(*region));
    return boxes;
}

template<std::size_t dim, std::size_t interp>
void checkExchangeRegion()
{
    using GridLayout          = GridLayout_t<dim, interp>;
    auto constexpr ghostWidth = core::ghostWidthForParticles<interp>();
    SAMRAI::tbox::Dimension const sdim{dim};
    SAMRAI::hier::IntVector const ghosts{sdim, ghostWidth};
    amr::ParticleDomainFromGhostFillPattern<GridLayout> pattern;

    // destination box, sources around it at every position of a small window,
    // directly or through a periodic shift
    Box_t<dim> const dstBox{core::ConstArray<int, dim>(0), core::ConstArray<int, dim>(3)};
    auto const dst = samraiBox(dstBox);
    SAMRAI::pdat::CellGeometry const dstGeometry{dst, ghosts};
    int constexpr period = 40;

    std::size_t nbrDisjoint = 0, nbrOverlapping = 0;
    std::array<int, dim> lower;
    auto const visit = [&](auto const& self, std::size_t d) -> void {
        if (d == dim)
        {
            for (int size : {1, 3})
                for (int offset : {0, period})
                {
                    Box_t<dim> srcBox;
                    for (std::size_t i = 0; i < dim; ++i)
                    {
                        srcBox.lower[i] = lower[i];
                        srcBox.upper[i] = lower[i] + size - 1;
                    }
                    Shift_t<dim> shift{};
                    SAMRAI::hier::IntVector offsetVector{sdim, 0};
                    for (std::size_t i = 0; i < dim; ++i)
                        shift[i] = offsetVector[i] = offset;
                    auto const src = samraiBox(core::shift(srcBox, shift * -1));
                    SAMRAI::hier::Transformation const transformation{offsetVector};
                    SAMRAI::pdat::CellGeometry const srcGeometry{src, ghosts};

                    // what RefineSchedule passes for a same-level schedule
                    auto fill_box = dst;
                    fill_box.grow(ghosts);
                    auto transformedSrc = src;
                    transformation.transform(transformedSrc);
                    auto src_mask = fill_box * transformedSrc;
                    if (src_mask.empty())
                        continue;
                    transformation.inverseTransform(src_mask);

                    auto const overlap = pattern.calculateOverlap(
                        dstGeometry, srcGeometry, dst, src_mask, fill_box, true, transformation);
                    auto const& boxes = dynamic_cast<amr::ParticlesDomainOverlap const&>(*overlap)
                                            .getDestinationBoxContainer();
                    ASSERT_LE(boxes.size(), 1u);

                    auto const region = cellsOf<dim>(boxes);
                    auto const shared = srcBox * dstBox;
                    if (!shared)
                    {
                        ++nbrDisjoint;
                        auto const reference = cellsOf<dim>(clippedRegion<GridLayout>(
                            dstGeometry, srcGeometry, src_mask, fill_box, transformation));
                        ASSERT_EQ(region, reference) << "source " << srcBox;
                    }
                    else
                    {
                        // the shared cells may hold particles that left the source
                        ++nbrOverlapping;
                        for (auto const& cell : *shared)
                            ASSERT_TRUE(region.count(cell.toArray())) << "source " << srcBox;
                    }
                }
            return;
        }
        for (int l = -5; l <= 6; ++l)
        {
            lower[d] = l;
            self(self, d + 1);
        }
    };
    visit(visit, 0);

    EXPECT_GT(nbrDisjoint, 0u);
    EXPECT_GT(nbrOverlapping, 0u);
}

TEST(PatchOverlapExchangeRegion, matchesClippedRegionWithoutOverlap1D)
{
    checkExchangeRegion<1, 1>();
    checkExchangeRegion<1, 2>();
    checkExchangeRegion<1, 3>();
}
TEST(PatchOverlapExchangeRegion, matchesClippedRegionWithoutOverlap2D)
{
    checkExchangeRegion<2, 1>();
    checkExchangeRegion<2, 2>();
    checkExchangeRegion<2, 3>();
}
TEST(PatchOverlapExchangeRegion, matchesClippedRegionWithoutOverlap3D)
{
    checkExchangeRegion<3, 1>();
    checkExchangeRegion<3, 2>();
    checkExchangeRegion<3, 3>();
}



// ---------------------------------------------------------------------------------------------
// border sum: a patch skips itself, not another patch with the same box

template<std::size_t dim>
auto makeLayout(SAMRAI::hier::Box const& box)
{
    using GridLayout = GridLayout_t<dim, 1>;
    std::array<double, dim> dl;
    std::array<std::uint32_t, dim> nbrCells;
    core::Point<double, dim> origin;
    for (std::size_t i = 0; i < dim; ++i)
    {
        dl[i]       = 0.1;
        nbrCells[i] = box.numberCells(i);
        origin[i]   = 0;
    }
    return GridLayout{dl, nbrCells, origin};
}

template<std::size_t dim>
void checkBorderSumSkipsOnlySelf()
{
    using GridLayout = GridLayout_t<dim, 1>;
    using Scalar     = core::HybridQuantity::Scalar;
    using Field_g    = amr::FieldGeometry<GridLayout, Scalar>;
    using Vector_g   = amr::TensorFieldGeometry<1, GridLayout, core::HybridQuantity>;

    SAMRAI::tbox::Dimension const sdim{dim};
    auto const extent
        = samraiBox(Box_t<dim>{core::ConstArray<int, dim>(0), core::ConstArray<int, dim>(7)});
    auto const layout = makeLayout<dim>(extent);

    SAMRAI::hier::Box const patch{extent, SAMRAI::hier::LocalId{0}, 0};
    SAMRAI::hier::Box const twin{extent, SAMRAI::hier::LocalId{1}, 0};
    SAMRAI::hier::Box const image{extent, SAMRAI::hier::LocalId{0}, 0, SAMRAI::hier::PeriodicId{1}};
    SAMRAI::hier::Transformation const noShift{SAMRAI::hier::IntVector::getZero(sdim)};
    auto fill_box = extent;
    fill_box.grow(SAMRAI::hier::IntVector{sdim, 2});

    auto const isEmpty = [&](auto const& pattern, auto const& dst, auto const& src) {
        // through the virtual interface, as RefineSchedule calls it
        using Geometry = SAMRAI::hier::BoxGeometry const&;
        return pattern
            .calculateOverlap(static_cast<Geometry>(dst), static_cast<Geometry>(src), patch, extent,
                              fill_box, true, noShift)
            ->isOverlapEmpty();
    };

    amr::FieldGhostInterpOverlapFillPattern<GridLayout> scalarPattern;
    Field_g const patchField{patch, layout, Scalar::rho};
    EXPECT_TRUE(isEmpty(scalarPattern, patchField, Field_g{patch, layout, Scalar::rho}));
    EXPECT_FALSE(isEmpty(scalarPattern, patchField, Field_g{twin, layout, Scalar::rho}));
    EXPECT_FALSE(isEmpty(scalarPattern, patchField, Field_g{image, layout, Scalar::rho}));

    auto constexpr V = core::HybridQuantity::Vector::V;
    amr::TensorFieldGhostInterpOverlapFillPattern<GridLayout> vectorPattern;
    Vector_g const patchVector{patch, layout, V};
    EXPECT_TRUE(isEmpty(vectorPattern, patchVector, Vector_g{patch, layout, V}));
    EXPECT_FALSE(isEmpty(vectorPattern, patchVector, Vector_g{twin, layout, V}));
}

TEST(PatchOverlapBorderSum, skipsSelfButNotIdenticalBox1D)
{
    checkBorderSumSkipsOnlySelf<1>();
}
TEST(PatchOverlapBorderSum, skipsSelfButNotIdenticalBox2D)
{
    checkBorderSumSkipsOnlySelf<2>();
}
TEST(PatchOverlapBorderSum, skipsSelfButNotIdenticalBox3D)
{
    checkBorderSumSkipsOnlySelf<3>();
}



// ---------------------------------------------------------------------------------------------
// ownership: every covered cell is owned by exactly one box, the one with the smallest id

template<std::size_t dim>
struct Configuration
{
    std::vector<Box_t<dim>> boxes;
    std::vector<int> ids;
};

template<std::size_t dim>
Box_t<dim> randomBox(std::mt19937& gen, int const domainSize)
{
    std::uniform_int_distribution<int> lowerDist(0, domainSize - 1);
    std::uniform_int_distribution<int> sizeDist(1, domainSize / 2);
    Box_t<dim> box;
    for (std::size_t i = 0; i < dim; ++i)
    {
        box.lower[i] = lowerDist(gen);
        box.upper[i] = std::min(domainSize - 1, box.lower[i] + sizeDist(gen) - 1);
    }
    return box;
}

template<std::size_t dim>
Configuration<dim> randomConfiguration(std::mt19937& gen, int const domainSize)
{
    Configuration<dim> config;
    std::uniform_int_distribution<int> nbrDist(2, 6);
    auto const nbrBoxes = nbrDist(gen);
    for (int i = 0; i < nbrBoxes; ++i)
        config.boxes.push_back(randomBox<dim>(gen, domainSize));

    // geometries that random boxes rarely produce
    std::uniform_int_distribution<int> pick(0, nbrBoxes - 1);
    auto const some = config.boxes[pick(gen)];
    config.boxes.push_back(some); // identical twin
    if (auto inner = some; some.upper[0] > some.lower[0])
    {
        inner.lower[0] += 1; // nested
        config.boxes.push_back(inner);
    }
    auto corner = some; // corner contact only
    for (std::size_t i = 0; i < dim; ++i)
        corner.lower[i] = corner.upper[i] = some.upper[i] + 1;
    config.boxes.push_back(corner);

    config.ids.resize(config.boxes.size());
    std::iota(config.ids.begin(), config.ids.end(), 0);
    std::shuffle(config.ids.begin(), config.ids.end(), gen);
    for (auto& id : config.ids)
        id = 7 * id + 3; // ids need not be contiguous
    return config;
}

template<std::size_t dim>
auto neighborsOf(Configuration<dim> const& config, std::size_t const self,
                 std::vector<Shift_t<dim>> const& shifts)
{
    std::vector<std::pair<Box_t<dim>, int>> neighbors;
    for (std::size_t j = 0; j < config.boxes.size(); ++j)
        for (auto const& shift : shifts)
            if (j != self or shift != Shift_t<dim>{}) // a patch is not its own neighbor
                neighbors.emplace_back(core::shift(config.boxes[j], shift), config.ids[j]);
    return neighbors;
}

template<std::size_t dim>
void checkEveryCoveredCellHasOneOwner(Configuration<dim> const& config, int const domainSize,
                                      std::vector<Shift_t<dim>> const& shifts)
{
    std::vector<std::vector<Box_t<dim>>> foreign;
    for (std::size_t i = 0; i < config.boxes.size(); ++i)
        foreign.push_back(
            amr::makeForeignBoxes(config.boxes[i], config.ids[i], neighborsOf(config, i, shifts)));

    Box_t<dim> domain;
    for (std::size_t i = 0; i < dim; ++i)
    {
        domain.lower[i] = 0;
        domain.upper[i] = domainSize; // corner boxes may stick out by one cell
    }

    for (auto const& cell : domain)
    {
        std::size_t nbrCovering = 0, nbrOwners = 0;
        int smallestId = std::numeric_limits<int>::max(), ownerId = -1;
        for (std::size_t i = 0; i < config.boxes.size(); ++i)
            if (core::isIn(cell, config.boxes[i]))
            {
                ++nbrCovering;
                smallestId = std::min(smallestId, config.ids[i]);
                if (!core::isIn(cell, foreign[i]))
                {
                    ++nbrOwners;
                    ownerId = config.ids[i];
                }
            }
        ASSERT_EQ(nbrOwners, nbrCovering > 0 ? 1u : 0u) << "cell " << cell;
        if (nbrCovering > 0)
        {
            ASSERT_EQ(ownerId, smallestId);
        }
    }
}

template<std::size_t dim>
void checkOwnership(int const domainSize, int const nbrConfigurations)
{
    std::mt19937 gen(1234 + dim);

    std::vector<Shift_t<dim>> const noShift{Shift_t<dim>{}};

    // periodic images: same id, shifted by one period in each direction
    auto periodic = noShift;
    for (std::size_t i = 0; i < dim; ++i)
        for (int sign : {-1, 1})
        {
            Shift_t<dim> shift{};
            shift[i] = sign * (domainSize + 1);
            periodic.push_back(shift);
        }

    for (int c = 0; c < nbrConfigurations; ++c)
    {
        auto const config = randomConfiguration<dim>(gen, domainSize);
        checkEveryCoveredCellHasOneOwner(config, domainSize, noShift);
        checkEveryCoveredCellHasOneOwner(config, domainSize, periodic);
        if (::testing::Test::HasFatalFailure())
            return;
    }
}

TEST(PatchOverlapOwnership, everyCoveredCellHasExactlyOneOwner1D)
{
    checkOwnership<1>(32, 1000);
}
TEST(PatchOverlapOwnership, everyCoveredCellHasExactlyOneOwner2D)
{
    checkOwnership<2>(16, 1000);
}
TEST(PatchOverlapOwnership, everyCoveredCellHasExactlyOneOwner3D)
{
    checkOwnership<3>(10, 1000);
}

TEST(PatchOverlapOwnership, smallerIdOwnsIdenticalBoxes)
{
    Box_t<2> const box{{0, 0}, {7, 7}};
    std::vector<std::pair<Box_t<2>, int>> const smaller{{box, 1}}, larger{{box, 2}};
    EXPECT_TRUE(amr::makeForeignBoxes(box, 1, larger).empty());
    auto const foreign = amr::makeForeignBoxes(box, 2, smaller);
    ASSERT_EQ(foreign.size(), 1u);
    EXPECT_EQ(foreign[0], box);
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();

    int testResult = RUN_ALL_TESTS();

    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
