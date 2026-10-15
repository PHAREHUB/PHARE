
#include "phare_core.hpp"


#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "amr/utilities/box/amr_box.hpp"

#include "simulator/simulator_def.hpp"

#include "tests/core/data/particles/test_particles.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"

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
    using GridLayout_t        = TestGridLayout<typename PhareTypes::Hybrid::GridLayout_t>;
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
    using ParticleArray_t             = TypeParam::TestParam::PhareTypes::Hybrid::ParticleArray_t;
    using GridLayout_t                = TypeParam::GridLayout_t;
    auto constexpr static dim         = TypeParam::dim;
    auto constexpr static ghost_cells = GridLayout_t::options.particle_ghost_width;

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
                core::add_ghost_particles(pop.patchGhostParticles(), box, ppc);
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

            if constexpr (ParticleArray_t::layout_mode == core::LayoutMode::AoSPCTS)
            {
                for (auto const& tile : pop.domainParticles()())
                    for (auto const& bix : tile().local_box())
                        check_array(tile()(bix));
            }
            else if constexpr (core::is_tiled(ParticleArray_t::layout_mode))
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

    // only level ghost particles are deposited: border nodes on a patch-patch boundary get
    // nothing, so restrict checks to nodes next to ghost zones no same-level neighbour covers
    auto const level_ghost_zones = [&](auto const& patch, auto const& amr_ghost_box,
                                       auto const& amr_domain) {
        auto zones = amr_ghost_box.remove(amr_domain);
        SAMRAI::hier::HierarchyNeighbors const neighbours{*hier, 1, 1};
        for (auto const& sam_neighbour : neighbours.getSameLevelNeighbors(patch.getBox(), 1))
        {
            auto const neighbour = grow(phare_box_from<dim>(sam_neighbour), 1); // cells -> nodes
            std::vector<core::Box<int, dim>> remaining;
            for (auto const& zone : zones)
                for (auto const& left : zone.remove(neighbour))
                    remaining.push_back(left);
            zones = std::move(remaining);
        }
        return zones;
    };

    for (auto& patch : *lvl1)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);
        GridLayout_t const layout{phare_box_from<dim>(patch->getBox())};

        for (auto& pop : ions)
        {
            auto const& density    = reduce(pop.particleDensity());
            auto const amr_domain  = layout.AMRBoxFor(density);
            auto const ghost_zones = level_ghost_zones(*patch, layout.AMRGhostBoxFor(density),
                                                       amr_domain);
            std::size_t checked = 0;
            for (auto const& border : amr_domain.remove(shrink(amr_domain, 1)))
                for (auto const& bix : border)
                {
                    auto const around = grow(core::Box<int, dim>{bix, bix}, 1);
                    if (std::none_of(ghost_zones.begin(), ghost_zones.end(),
                                     [&](auto const& zone) { return bool{around * zone}; }))
                        continue;
                    EXPECT_GT(density(layout.AMRToLocal(bix)), 0.) << "node " << bix;
                    ++checked;
                }
            if (ghost_zones.size())
            {
                EXPECT_GT(checked, 0u);
            }
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
