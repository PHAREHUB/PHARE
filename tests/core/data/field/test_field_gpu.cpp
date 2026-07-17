#include "phare_core.hpp"
// #include "test_gridlayout.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <vector>
#include <cstddef>

using namespace PHARE::core;

template<std::size_t dim>
using NdArray_t = typename PHARE_Types<PHARE::SimOpts{dim}>::Array_t;

template<typename T>
using ManagedVector = std::vector<T, mkn::gpu::ManagedAllocator<T>>;

TEST(FieldGPUTest, copyWorksOnGPU)
{
    static_assert(std::is_same_v<PHARE::Allocator::allocator_type<double>,
                                 mkn::gpu::ManagedAllocator<double>>);

    std::size_t constexpr static N       = 1024;
    std::size_t constexpr static dim     = 3;
    std::size_t constexpr static interp  = 1;
    std::uint32_t constexpr static cells = 30;
    using Field_t                        = Field<NdArray_t<dim>, HybridQuantity::Scalar>;

    Field_t rho{"test", HybridQuantity::Scalar::rho, ConstArray<std::uint32_t, dim>(cells)};
    auto v_rho = rho.view();
    assert(mkn::gpu::Pointer{rho.data()}.is_managed_ptr());
    assert(mkn::gpu::Pointer{v_rho.data()}.is_managed_ptr());
    KLOG(INF) << rho.size();
    for (std::size_t i = 0; i < N; ++i)
        rho.data()[i] = i;

    {
        mkn::gpu::Pointer a{rho.data()};
        KLOG(INF) << a.attributes.memoryType;
        KLOG(INF) << a.attributes.isManaged;
        KLOG(INF) << a.attributes.devicePointer;
        KLOG(INF) << a.attributes.hostPointer;
        KLOG(INF) << a.is_host_ptr();
        KLOG(INF) << a.is_managed_ptr();
        KLOG(INF) << a.is_device_ptr();
    }

    ManagedVector<std::array<std::uint32_t, 3>> icells(N, ConstArray<std::uint32_t, dim>(1));

    auto rd = v_rho.data();
    // mkn::gpu::GDLauncher{N}([=] __device__(auto& r) { r.data()[mkn::gpu::idx()] += 1; }, v_rho);
    mkn::gpu::GDLauncher{N}([=] __device__() { rd[mkn::gpu::idx()] += 1; });

    auto expect = ConstArray<std::uint32_t, dim>(4); // 2 + 2 ghostnodes
    for (std::size_t i = 0; i < N; ++i)
        EXPECT_EQ(rho.data()[i], i + 1);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
