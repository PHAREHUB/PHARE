// Same-level patches may overlap. These tests cover the pieces that make overlapping
// patches behave as a partition of the level: the region used to exchange leaving particles.

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
