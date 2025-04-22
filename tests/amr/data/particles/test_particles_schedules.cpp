
#include "phare_core.hpp"
#include "phare/phare.hpp"
#include "core/utilities/meta/meta_utilities.hpp"
#include <amr/utilities/box/amr_box.hpp>

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"

#include "gtest/gtest.h"

#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>
#include <core/utilities/box/box.hpp>
#include <core/utilities/types.hpp>

static constexpr std::size_t interp = 1;
static constexpr std::size_t ppc    = 100;

template<std::size_t _dim>
struct TestParam
{
    auto constexpr static dim = _dim;
    using PhareTypes          = PHARE::core::PHARE_Types<dim, interp>;
    using GridLayout_t        = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using Box_t               = PHARE::core::Box<int, dim>;
    using ParticleArray_t     = PHARE::core::ParticleArray<dim>;

    using Hierarchy_t
        = AfullHybridBasicHierarchy<dim, interp, defaultNbrRefinedParts<dim, interp>()>;
};

template<typename TestParam>
struct ParticleScheduleHierarchyTest : public ::testing::Test
{
    auto constexpr static dim = TestParam::dim;
    using Hierarchy_t         = TestParam::Hierarchy_t;
    using GridLayout_t        = TestParam::GridLayout_t;
    using ResourceManager_t   = Hierarchy_t::ResourcesManagerT;
    using DomainGhostPartRefinerPool
        = RefinerPool<ResourceManager_t, RefinerType::ExteriorGhostParticles>;

    std::string configFile = "test_particles_schedules_inputs/" + std::to_string(dim) + "d_L0.txt";
    Hierarchy_t hierarchy{configFile};
};

using ParticlesDatas = testing::Types<TestParam<1>, TestParam<2> /*,TestParam<3>*/>;


TYPED_TEST_SUITE(ParticleScheduleHierarchyTest, ParticlesDatas);
namespace PHARE::core
{


TYPED_TEST(ParticleScheduleHierarchyTest, testing_inject_ghost_layer)
{
    using GridLayout_t                = TypeParam::GridLayout_t;
    auto constexpr static dim         = TypeParam::dim;
    auto constexpr static ghost_cells = GridLayout_t::nbrParticleGhosts();

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    for (auto& patch : *lvl0)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
        {
            pop.domainParticles().clear();
            EXPECT_EQ(pop.domainParticles().size(), 0);

            GridLayout_t const layout{phare_box_from<dim>(patch->getBox())};
            auto const ghostBox = grow(layout.AMRBox(), ghost_cells);
            for (auto const& box : ghostBox.remove(layout.AMRBox()))
                core::add_particles_in(pop.patchGhostParticles(), box, ppc);
        }
        rm.setTime(ions, *patch, 1);
    }

    // domainGhostPartRefiners.fill(0, 0);

    this->hierarchy.messenger->fillIonGhostParticles(ions, *lvl0, 0);

    for (auto& patch : *lvl0)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);

        SAMRAI::hier::HierarchyNeighbors const hier_nbrs{
            *this->hierarchy.basicHierarchy->hierarchy(), patch->getPatchLevelNumber(),
            patch->getPatchLevelNumber()};

        auto domainSamBox    = patch->getBox();
        auto const domainBox = phare_box_from<dim>(domainSamBox);

        auto const neighbors = core::generate(
            [](auto const& el) { return phare_box_from<dim>(el); },
            hier_nbrs.getSameLevelNeighbors(domainSamBox, patch->getPatchLevelNumber()));
        auto const ncells = core::sum_from(
            neighbors, [&](auto& el) { return (*(grow(el, ghost_cells) * domainBox)).size(); });

        for (auto& pop : ions)
        {
            EXPECT_EQ(pop.domainParticles().size(), ncells * ppc);

            for (auto const& p : pop.domainParticles())
            {
                EXPECT_TRUE(isIn(p, domainBox));
            }
        }
    }
}



} // namespace PHARE::core


int main(int argc, char** argv)
{
    PHARE::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
