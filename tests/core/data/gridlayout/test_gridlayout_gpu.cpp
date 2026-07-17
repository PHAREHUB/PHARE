#include "phare_core.hpp"
#include "test_gridlayout.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

template<typename T>
using ManagedVector = std::vector<T, mkn::gpu::ManagedAllocator<T>>;

TEST(GridLayoutGPUTest, copyWorksOnGPU)
{
    std::size_t constexpr static N       = 1024;
    std::size_t constexpr static dim     = 3;
    std::size_t constexpr static interp  = 1;
    std::uint32_t constexpr static cells = 30;

    using PHARE_Types  = PHARE::core::PHARE_Types<SimOpts{dim, interp}>;
    using GridLayout_t = typename PHARE_Types::GridLayout_t;

    TestGridLayout<GridLayout_t> layout{cells};

    ManagedVector<std::array<std::uint32_t, 3>> icells(N, ConstArray<std::uint32_t, 3>(1));

    mkn::gpu::GDLauncher{N}(
        [=] __device__(auto ics, auto& l) {
            ics[mkn::gpu::idx()] = l.AMRToLocal(Point{ConstArray<int, 3>(2)}).toArray();
        },
        icells, layout);

    auto expect = ConstArray<std::uint32_t, 3>(4); // 2 + 2 ghostnodes
    for (std::size_t i = 0; i < N; ++i)
        EXPECT_EQ(icells[i], expect);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
