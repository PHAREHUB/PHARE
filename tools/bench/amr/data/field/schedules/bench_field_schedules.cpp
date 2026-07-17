
#include "phare_core.hpp"

#include <core/utilities/types.hpp>
#include <core/utilities/box/box.hpp>

#include "tests/amr/amr.hpp"
#include "tests/amr/test_hierarchy_fixtures.hpp"

#include "gtest/gtest.h"

namespace PHARE::amr
{

static constexpr std::size_t ppc = 0; // not a particle test


template<SimOpts opts>
struct TestParam
{
    auto constexpr static dim    = opts.dimension;
    auto constexpr static interp = opts.interp_order;

    using PhareTypes  = core::PHARE_Types<opts>;
    using Hierarchy_t = AfullHybridBasicHierarchy<opts>;
};

template<typename TestParam_>
struct FieldScheduleHierarchy : public ::testing::Test
{
    using TestParam           = TestParam_;
    using Hierarchy_t         = typename TestParam::Hierarchy_t;
    auto constexpr static dim = TestParam::dim;

    std::string const configFile = "bench_field_schedules/config_" + std::to_string(dim) + "d.txt";
    Hierarchy_t hierarchy{configFile, ppc};
};


// clang-format off
using FieldDatas = testing::Types<
   //  TestParam<SimOpts{}>
   // ,TestParam<SimOpts{2}>
   // ,TestParam<SimOpts{3}>
// PHARE_WITH_MKN_GPU(
   // ,TestParam<SimOpts{.layout_mode=LayoutMode::AoSTS}>
   // ,TestParam<SimOpts{.dimension=2, .layout_mode=LayoutMode::AoSTS}>
   /*,*/TestParam<SimOpts{.dimension=3, .layout_mode=LayoutMode::AoSTS}>
// )

>;
// clang-format on

TYPED_TEST_SUITE(FieldScheduleHierarchy, FieldDatas, );

TYPED_TEST(FieldScheduleHierarchy, fillElectricGhosts)
{
    PHARE_LOG_SCOPE(1, "FieldScheduleHierarchy::fillElectricGhosts");
    auto& electromag = this->hierarchy.hybridModel->state.electromag;

    this->hierarchy.messenger->fillElectricGhosts(electromag.E, /*lvl=*/0, /*time=*/0);
}

// TYPED_TEST(FieldScheduleHierarchy, fillDensityBorders)
// {
//     auto& electromag = this->hierarchy.hybridModel->state.electromag;

//     auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
//     auto& ions = this->hierarchy.hybridModel->state.ions;

//     this->hierarchy.messenger->fillDensityBorders(ions, *lvl0, 0);
// }

// TYPED_TEST(FieldScheduleHierarchy, fillFluxBorders)
// {
//     auto& electromag = this->hierarchy.hybridModel->state.electromag;

//     auto lvl0  = this->hierarchy.basicHierarchy->hierarchy()->getPatchLevel(0);
//     auto& ions = this->hierarchy.hybridModel->state.ions;

//     this->hierarchy.messenger->fillFluxBorders(ions, *lvl0, 0);
// }


} // namespace PHARE::amr


int main(int argc, char** argv)
{
    PHARE_WITH_PHLOP( //
        phlop::threaded::ScopeTimerMan::INSTANCE()
            .file_name(
                ".phare/timings/bench/amr/data/field/schedules/bench_field_schedules.cpp.txt")
            .init(); //
    )
    PHARE::test::amr::SamraiLifeCycle samsam{argc, argv};
    ::testing::InitGoogleTest(&argc, argv);
    // PHARE_LOG_SCOPE(1, "WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    auto const ret = RUN_ALL_TESTS();
    PHARE_WITH_PHLOP(phlop::threaded::ScopeTimerMan::reset());
    return ret;
}
