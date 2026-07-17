// expects to be included
#include <thrust/partition.h>
#include <thrust/execution_policy.h>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>


// TYPED_TEST(VectorTest, thrust_gpu_vector_can_push_back_atomically)
// {
//     std::size_t static constexpr N = 1024;

//     thrust::universal_vector<int> vec(N);

//     mkn::gpu::GDLauncher{N}([=] _PHARE_ALL_FN_() mutable {
//         auto idx = mkn::gpu::idx();
//         // fv[mkn::gpu::idx()] = cf(ps.iCell(idx));
//         if (N % 100)
//             vec.push_back(idx);
//     });

//     EXPECT_EQ(vec.size(), N + 102);
// }


TYPED_TEST(VectorTest, is_unified_usable_on_host)
{
    static_assert(TypeParam::allocator_mode == PHARE::AllocatorMode::GPU_UNIFIED);
    auto vec = TypeParam::make(10);
    vec[0]   = 2;
    EXPECT_EQ(PHARE::core::sum(vec), 2);
}


TYPED_TEST(VectorTest, is_host_to_device_copyable)
{
    using T = typename TypeParam::value_type;

    std::vector<T> vec0(10, 1);
    auto vec1 = TypeParam::make(10);
    EXPECT_EQ(PHARE::core::sum(vec1), 0);

    TypeParam::copy(vec1, vec0);
    EXPECT_EQ(PHARE::core::sum(vec1), 10);
}

TYPED_TEST(VectorTest, is_device_to_host_copyable)
{
    using T = typename TypeParam::value_type;

    auto vec0 = TypeParam::make(10);
    TypeParam::fill(vec0, 1);
    EXPECT_EQ(PHARE::core::sum(vec0), 10);

    std::vector<T> vec1(10, 0);
    EXPECT_EQ(PHARE::core::sum(vec1), 0);
    PHARE::Vector<T>::copy(vec1, vec0);
    EXPECT_EQ(PHARE::core::sum(vec1), 10);
}

struct is_even
{
    __device__ bool operator()(const int& x) { return (x % 2) == 0; }
};

TYPED_TEST(VectorTest, thrust_device_partition)
{
    auto constexpr N = 1024;
    auto A           = TypeParam::make(N);

    assert(mkn::gpu::Pointer{A.data()}.is_managed_ptr());
    for (std::size_t i = 0; i < N; ++i)
        A[i] = i;

    // there might be other versions that work
    thrust::stable_partition(A.data(), A.data() + N, is_even());

    for (std::size_t i = 0; i < N / 2; ++i)
        EXPECT_EQ(static_cast<int>(A[i]) % 2, 0);
    for (std::size_t i = N / 2; i < N; ++i)
        EXPECT_EQ(static_cast<int>(A[i]) % 2, 1);
}
