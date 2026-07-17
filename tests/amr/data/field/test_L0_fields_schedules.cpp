
#include "phare_core.hpp"
// #include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/field/field_overlaps.hpp"
#include "core/data/particles/particle_array_def.hpp"

#include <core/utilities/types.hpp>
#include <core/utilities/box/box.hpp>
// #include "core/data/grid/grid_tiles.hpp"
#include <core/data/ndarray/ndarray_vector.hpp>

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"
#include "tests/core/data/gridlayout/test_gridlayout.hpp"


// #include <SAMRAI/pdat/CellGeometry.h>
// #include <SAMRAI/hier/HierarchyNeighbors.h>


#include "gtest/gtest.h"

namespace PHARE::amr
{

static constexpr std::size_t ppc = 100;


template<SimOpts _opts>
struct TestParam
{
    auto constexpr static opts   = _opts;
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
        = "test_fields_schedules_inputs/" + std::to_string(dim) + "d_L0.txt";
    Hierarchy_t hierarchy{configFile, ppc};
};


// clang-format off
using FieldDatas = testing::Types<
    TestParam<SimOpts{1}>,
    TestParam<SimOpts{2}>,
    TestParam<SimOpts{3}>
// PHARE_WITH_MKN_GPU(
//    ,TestParam<SimOpts{.dimension=1, .layout_mode=LayoutMode::AoSTS}>
//    ,TestParam<SimOpts{.dimension=2, .layout_mode=LayoutMode::AoSTS}>
//    ,TestParam<SimOpts{.dimension=3, .layout_mode=LayoutMode::AoSTS}>
// )

>;
// clang-format on

TYPED_TEST_SUITE(FieldScheduleHierarchyTest, FieldDatas, );

template<typename GridLayout_t>
void setup_field(GridLayout_t const& layout, auto& field)
{
    auto constexpr static dim        = std::decay_t<decltype(field)>::dimension;
    auto constexpr static field_opts = FieldOpts<HybridQuantity::Scalar, double>{dim};
    using FieldOverlaps              = FieldTileOverlaps<GridLayout_t, field_opts>;
    FieldOverlaps::getOrCreateQuantity(layout, field);
}

void setup_vecfield(auto const& layout, auto& vecfield)
{
    for (auto& field : vecfield)
        setup_field(layout, field);
}

TYPED_TEST(FieldScheduleHierarchyTest, testing_hyhy_schedules)
{
    PHARE_FN_TIMER("FieldScheduleHierarchyTest::testing_hyhy_schedules");
    auto constexpr static dim = TypeParam::dim;
    using TestParam           = TestFixture::TestParam;
    using ParticleArray_t     = TestParam::ParticleArray_t;
    using GridLayout_t        = TestParam::GridLayout_t;

    using FieldData_t = TestFixture::ResourceManager_t::UserField_t::patch_data_type;
    using Interpolating_t
        = core::Interpolating<ParticleArray_t, TestParam::interp, /*atomic_interp*/ false>;

    auto constexpr function_id = join_string_views_v<>;

    auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
    auto& rm   = *this->hierarchy.resourcesManagerHybrid;
    auto& ions = this->hierarchy.hybridModel->state.ions;

    Interpolating_t interpolate;

    for (auto patch : rm.enumerate(*lvl0, ions))
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        resetMoments(ions);
        core::depositParticles(ions, layout, interpolate, core::DomainDeposit{});
        if constexpr (TestParam::opts.layout_mode == LayoutMode::AoSTS)
        {
            for (auto& pop : ions)
            {
                setup_field(layout, pop.particleDensity());
                setup_vecfield(layout, pop.flux());
                PHARE_LOG_LINE_SS(sum_field(reduce(pop.particleDensity())));
            }
        }
    }

    PHARE_LOG_SCOPE(1, ParticleArray_t::type_id);
    PHARE_FN_TIMER("FieldScheduleHierarchyTest::testing_hyhy_schedules::schedules");
    this->hierarchy.messenger->fillDensityBorders(ions, *lvl0, 0);
    this->hierarchy.messenger->fillFluxBorders(ions, *lvl0, 0);

    for (auto patch : rm.enumerate(*lvl0, ions))
    {
        auto const layout = PHARE::amr::layoutFromPatch<GridLayout_t>(*patch);
        PHARE_LOG_LINE_SS(layout.AMRBox());

        for (auto const& pop : ions)
        {
            PHARE_LOG_LINE_SS(sum_field(reduce_single(pop.particleDensity())));
        }
    }
}




} // namespace PHARE::amr

int main(int argc, char** argv)
{
    PHARE_WITH_PHLOP({
        phlop::threaded::ScopeTimerMan::INSTANCE()
            .file_name(".phare/timings/tests/amr/data/field/test_L0_fields_schedules.cpp.txt")
            .init();
    })
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    auto const ret = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
    return ret;
}
