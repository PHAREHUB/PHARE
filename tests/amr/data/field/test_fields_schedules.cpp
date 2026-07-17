

#include "phare_core.hpp"

#include "core/utilities/types.hpp"
#include "core/utilities/box/box.hpp"
#include "core/data/grid/grid_tiles.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/tensorfield/tensorfield.hpp"

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"


#include <SAMRAI/pdat/CellGeometry.h>
#include <SAMRAI/hier/HierarchyNeighbors.h>


#include "gtest/gtest.h"

namespace PHARE::amr
{

static constexpr std::size_t ppc = 100;


template<SimOpts opts>
struct TestParam
{
    auto constexpr static dim    = opts.dimension;
    auto constexpr static interp = opts.interp_order;

    using PhareTypes       = core::PHARE_Types<opts>;
    using Field_t          = PhareTypes::Field_t;
    using GridLayout_t     = PhareTypes::GridLayout_t;
    using TestGridLayout_t = TestGridLayout<GridLayout_t>;
    using Box_t            = core::Box<int, dim>;
    using ParticleArray_t  = PhareTypes::ParticleArray_t;

    using Hierarchy_t = AfullHybridBasicHierarchy<opts>;
};

template<typename TestParam_>
struct FieldScheduleHierarchyTest : public ::testing::Test
{
    using TestParam           = TestParam_;
    using Hierarchy_t         = typename TestParam::Hierarchy_t;
    using ResourceManager_t   = typename Hierarchy_t::ResourcesManagerT;
    auto constexpr static dim = TestParam::dim;

    std::string const configFile
        = "test_fields_schedules_inputs/" + std::to_string(dim) + "d_config.txt";
    Hierarchy_t hierarchy{configFile};
};


// clang-format off
using FieldDatas = testing::Types<
    TestParam<SimOpts{}>
   // ,TestParam<SimOpts{2}>
PHARE_WITH_MKN_GPU(
   ,TestParam<SimOpts{.dimension=1,.layout_mode=LayoutMode::AoSTS}>
   // ,TestParam<SimOpts{.dimension=2,.layout_mode=LayoutMode::AoSTS}>
   // ,TestParam<SimOpts{.dimension=2, .layout_mode=LayoutMode::AoSTS}>
)

>;
// clang-format on

TYPED_TEST_SUITE(FieldScheduleHierarchyTest, FieldDatas, );


TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_schedules)
{
    auto constexpr static dim = TypeParam::dim;
    using TestParam           = TestFixture::TestParam;
    using ParticleArray_t     = TestParam::ParticleArray_t;
    using GridLayout_t        = TestParam::GridLayout_t;
    using Interpolating_t
        = core::Interpolating<ParticleArray_t, TestParam::interp, /*atomic_interp*/ false>;


    auto constexpr static interp      = GridLayout_t::interp_order;
    auto constexpr static ghost_cells = GridLayout_t::nbrGhosts();

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    Interpolating_t interpolate;

    for (auto& patch : *lvl0)
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, ions);
        resetMoments(ions);

        core::depositParticles(ions, layout, interpolate, core::DomainDeposit{});

        for (auto& pop : ions)
            if (pop.name() == "protons")
            {
                PHARE_LOG_LINE_SS(core::sum_field(reduce(pop.particleDensity())));
            }
    }

    this->hierarchy.messenger->fillDensityBorders(ions, *lvl0, 0);
    this->hierarchy.messenger->fillFluxBorders(ions, *lvl0, 0);


    for (auto& patch : *lvl0)
    {
        auto dataOnPatch = rm.setOnPatch(*patch, ions);

        if (mpi::rank() == 0)
            for (auto& pop : ions)
                if (pop.name() == "protons")
                {
                    auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
                    auto const domGhostBox
                        = layout.AMRGhostBoxFor(pop.flux()[0].physicalQuantity());

                    PHARE_LOG_LINE_SS(core::sum_field(reduce_single(pop.particleDensity())));
                }
    }
}




TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_field_refine_schedules)
{
    auto constexpr static dim = TypeParam::dim;
    using GridLayout_t        = TestFixture::TestParam::GridLayout_t;

    auto constexpr static interp      = GridLayout_t::interp_order;
    auto constexpr static ghost_cells = GridLayout_t::nbrGhosts();

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto lvl1  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(1);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& em   = this->hierarchy.hybridModel->state.electromag;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    auto has_nan = [](auto const& field) { return std::isnan(core::sum_field(field)); };

    for (auto& patch : rm.enumerate(*lvl0, *this->hierarchy.hybridModel))
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        fill_domain_box(em.E[0], layout, 1);
        fill_domain_box(em.E[1], layout, 1);
        fill_domain_box(em.E[2], layout, 1);
    }

    for (auto& patch : rm.enumerate(*lvl0, *this->hierarchy.hybridModel))
    {
        if (mpi::rank() == 0)
        {
            auto const layout = layoutFromPatch<GridLayout_t>(*patch);
            PHARE_LOG_LINE_SS(layout.AMRBox());
            PHARE_LOG_LINE_SS(em.E);
        }
    }

    this->hierarchy.messenger->fillElectricGhosts(em.E, *lvl0, 0);
    this->hierarchy.messenger->fillMagneticGhosts(em.B, *lvl0, 0);

    for (auto& patch : rm.enumerate(*lvl0, *this->hierarchy.hybridModel))
    {
        if (mpi::rank() == 0)
        {
            auto const layout = layoutFromPatch<GridLayout_t>(*patch);
            PHARE_LOG_LINE_SS(layout.AMRBox());
            PHARE_LOG_LINE_SS(em.E);
        }
    }

    // for (auto& patch : rm.enumerate(*lvl0, *this->hierarchy.hybridModel))
    //     for (auto& field : em.E)
    //         EXPECT_FALSE(has_nan(field)) << "NaN in lvl0 E ghost cells after fillElectricGhosts";


    for (auto& patch : rm.enumerate(*lvl1, *this->hierarchy.hybridModel))
    {
        // em.E[0].fill(0);
        // em.E[1].fill(0);
        // em.E[2].fill(0);
        if (mpi::rank() == 0)
        {
            // PHARE_LOG_LINE_SS(em.B);
            // PHARE_LOG_LINE_SS(em.E);
        }
        // check_tensor_field(electromag.B);
    }

    this->hierarchy.messenger->fillElectricGhosts(em.E, *lvl1, 0);
    this->hierarchy.messenger->fillMagneticGhosts(em.B, *lvl1, 0);

    for (auto& patch : rm.enumerate(*lvl1, *this->hierarchy.hybridModel))
    {
        auto const layout = layoutFromPatch<GridLayout_t>(*patch);
        if (mpi::rank() == 0)
        {
            PHARE_LOG_LINE_SS(layout.AMRBox());
            // PHARE_LOG_LINE_SS(em.B);
            PHARE_LOG_LINE_SS(em.E);
        }

        // check_tensor_field(electromag.B);
    }
}



} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
