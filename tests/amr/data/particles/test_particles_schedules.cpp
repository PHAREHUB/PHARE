
#include "phare_core.hpp"


#include <core/utilities/types.hpp>
#include <core/utilities/box/box.hpp>
#include <amr/utilities/box/amr_box.hpp>

#include "tests/core/data/gridlayout/test_gridlayout.hpp"
#include "tests/core/data/particles/test_particles_fixtures.hpp"

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"



#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>

#include "gtest/gtest.h"

namespace PHARE::amr
{

static constexpr std::size_t ppc = 100;

template<auto opts>
struct TestParam
{
    auto constexpr static dim = opts.dimension;
    using PhareTypes          = PHARE::core::PHARE_Types<opts>;
    using GridLayout_t        = TestGridLayout<typename PhareTypes::GridLayout_t>;
    using Hierarchy_t         = AfullHybridBasicHierarchy<opts>;
};

template<typename TestParam_>
struct ParticleScheduleHierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    auto constexpr static dim = TestParam::dim;
    using Hierarchy_t         = TestParam::Hierarchy_t;
    using GridLayout_t        = TestParam::GridLayout_t;
    using ResourceManager_t   = Hierarchy_t::ResourcesManagerT;

    std::string configFile = "test_particles_schedules_inputs/" + std::to_string(dim) + "d_L0.txt";
    Hierarchy_t hierarchy{configFile};
};

// clang-format off
using ParticlesDatas = testing::Types<
   TestParam<SimOpts{}>

PHARE_WITH_MKN_GPU(
  ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSTS}>
  ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSPCTS}>
)

>;
// clang-format on


TYPED_TEST_SUITE(ParticleScheduleHierarchyTest, ParticlesDatas);


TYPED_TEST(ParticleScheduleHierarchyTest, testing_inject_ghost_layer)
{
    using ParticleArray_t             = TypeParam::TestParam::PhareTypes::ParticleArray_t;
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
            if constexpr (core::is_tiled(ParticleArray_t::layout_mode))
                for (auto const& tile : pop.domainParticles()())
                {
                    EXPECT_EQ(tile().size(), 0);
                }

            GridLayout_t const layout{phare_box_from<dim>(patch->getBox())};
            auto const ghostBox = grow(layout.AMRBox(), ghost_cells);
            for (auto const& box : ghostBox.remove(layout.AMRBox()))
                core::add_ghost_particles(pop.patchGhostParticles(), layout, box, ppc);
        }
        rm.setTime(ions, *patch, 1);
    }

    this->hierarchy.messenger->fillIonGhostParticles(ions, *lvl0, 0);

    auto n_ghost_cells_for_neighbours = [&](auto& patch) {
        auto domainSamBox    = patch->getBox();
        auto const domainBox = phare_box_from<dim>(domainSamBox);
        return core::sum_from(
            core::generate_from(
                [](auto const& el) { return phare_box_from<dim>(el); },
                SAMRAI::hier::HierarchyNeighbors{*this->hierarchy.basicHierarchy->hierarchy(),
                                                 patch->getPatchLevelNumber(),
                                                 patch->getPatchLevelNumber()}
                    .getSameLevelNeighbors(domainSamBox, patch->getPatchLevelNumber())),
            [&](auto& el) { return (*(grow(el, ghost_cells) * domainBox)).size(); });
    };

    for (auto& patch : *lvl0)
    {
        auto dataOnPatch       = rm.setOnPatch(*patch, ions);
        auto const domainBox   = phare_box_from<dim>(patch->getBox());
        auto const check       = [&](auto const& p) { EXPECT_TRUE(isIn(p, domainBox)); };
        auto const check_array = [&](auto const& array) {
            for (auto const& p : array)
                check(p);
        };
        auto const ncells = n_ghost_cells_for_neighbours(patch);
        for (auto& pop : ions)
        {
            EXPECT_EQ(pop.domainParticles().size(), ncells * ppc);

            if constexpr (core::is_tiled(ParticleArray_t::layout_mode))
                for (auto const& tile : pop.domainParticles()())
                    check_array(tile());
            else
                check_array(pop.domainParticles());
        }
    }
}

template<typename TestParam_>
struct ParticleScheduleL1HierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    auto constexpr static dim = TestParam::dim;
    using Hierarchy_t         = TestParam::Hierarchy_t;

    std::string configFile
        = "test_particles_schedules_inputs/" + std::to_string(dim) + "d_config.txt";
    Hierarchy_t hierarchy{configFile};
};

// clang-format off
using ParticlesDatasL1 = testing::Types<
    TestParam<SimOpts{}>
   ,TestParam<SimOpts{2}>
   // ,TestParam<SimOpts{3}>

PHARE_WITH_MKN_GPU(
  ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSTS}>
  ,TestParam<SimOpts{.dimension=2,.layout_mode=LayoutMode::AoSTS}>
  // ,TestParam<SimOpts{.dimension=3,.layout_mode=LayoutMode::AoSTS}>
)

>;
// clang-format on

TYPED_TEST_SUITE(ParticleScheduleL1HierarchyTest, ParticlesDatasL1);


TYPED_TEST(ParticleScheduleL1HierarchyTest, fillIonPopMomentGhostsContributesToBoundaryNodes)
{
    using GridLayout_t        = TypeParam::GridLayout_t;
    auto constexpr static dim = TypeParam::dim;

    auto& rm          = *this->hierarchy.resourcesManagerHybrid;
    auto& ions        = this->hierarchy.hybridModel->state.ions;
    auto& hybridModel = *this->hierarchy.hybridModel;
    auto hier         = this->hierarchy.basicHierarchy->hierarchy();

    auto lvl1 = hier->getPatchLevel(1);

    for (auto& patch : *lvl1)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        for (auto& pop : ions)
            EXPECT_GT(pop.levelGhostParticlesOld().size(), 0u);
    }

    this->hierarchy.messenger->firstStep(hybridModel, *lvl1, hier, 0., 0., 1.);

    for (auto& patch : *lvl1)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        resetMoments(ions);
    }

    this->hierarchy.messenger->fillIonPopMomentGhosts(ions, *lvl1, 0.);

    for (auto& patch : *lvl1)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        GridLayout_t const layout{phare_box_from<dim>(patch->getBox())};

        for (auto& pop : ions)
        {
            auto const& density            = reduce(pop.particleDensity());
            auto const local_domain        = layout.domainBoxFor(density);
            auto const local_domain_border = local_domain.remove(shrink(local_domain, 1));
            // ghost zones in local coords: avoids SAMRAI periodic-wrap in AMR indices
            auto const local_ghost_zones = layout.ghostBoxFor(density).remove(local_domain);
            std::size_t checked          = 0;
            for (auto const& border : local_domain_border)
            {
                // grow border outward (lower always >= nghosts >= 1, no uint32_t underflow)
                auto const grown_border = grow(border, 1);
                for (auto const& ghost : local_ghost_zones)
                {
                    if (grown_border * ghost)
                        for (auto const& bix : border)
                        {
                            EXPECT_GT(density(bix), 0.);
                            ++checked;
                        }
                }
            }
            EXPECT_GE(checked, local_domain.size() - shrink(local_domain, 1).size());
        }
    }
}

} // namespace PHARE::amr



int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
