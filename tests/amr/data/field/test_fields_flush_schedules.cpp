
#include "phare_core.hpp"

#include <core/utilities/types.hpp>
#include <core/utilities/box/box.hpp>
#include "core/data/grid/grid_tiles.hpp"
#include <core/data/ndarray/ndarray_vector.hpp>

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
        = "test_fields_schedules_inputs/" + std::to_string(dim) + "d_flush.txt";
    Hierarchy_t hierarchy{configFile};
};


// clang-format off
using FieldDatas = testing::Types<
    // TestParam<SimOpts{}>
   /*,*/TestParam<SimOpts{2}>
PHARE_WITH_MKN_GPU(
   // ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSTS}>
   // ,TestParam<SimOpts{.dimension=2, .layout_mode=LayoutMode::AoSTS}>
)

>;
// clang-format on

TYPED_TEST_SUITE(FieldScheduleHierarchyTest, FieldDatas, );



TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_field_refine_schedules)
{
    auto constexpr static dim    = TypeParam::dim;
    using GridLayout_t           = TestFixture::TestParam::GridLayout_t;
    using Field_t                = TestFixture::TestParam::Field_t;
    using FieldData_t            = TestFixture::ResourceManager_t::UserField_t::patch_data_type;
    auto constexpr static interp = GridLayout_t::interp_order;
    auto constexpr static ghost_cells = GridLayout_t::nbrGhosts();

    auto lvl0        = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto lvl1        = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(1);
    auto lvl2        = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(2);
    auto& rm         = *this->hierarchy.resourcesManagerHybrid;
    auto& electromag = this->hierarchy.hybridModel->state.electromag;
    auto& ions       = this->hierarchy.hybridModel->state.ions;


    for (auto& patch : *lvl0)
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, electromag);

        electromag.E[0].fill(100);
        electromag.E[1].fill(100);
    }


    for (auto& patch : *lvl1)
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, electromag);

        electromag.E[0].fill(200);
        electromag.E[1].fill(200);
    }

    for (auto& patch : *lvl2)
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        auto dataOnPatch  = rm.setOnPatch(*patch, electromag);

        electromag.E[0].fill(300);
        electromag.E[1].fill(300);
    }

    this->hierarchy.messenger->fillElectricGhosts(electromag.E, 2, 0);

    for (auto& patch : *lvl2)
    {
        PHARE_LOG_LINE_SS(patch->getGlobalId());
        auto dataOnPatch = rm.setOnPatch(*patch, *this->hierarchy.hybridModel);

        if (mpi::rank() == 0)
        {
            PHARE_LOG_LINE_SS("\n" << electromag.E[0]);
            PHARE_LOG_LINE_SS("\n" << electromag.E[1]);
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
