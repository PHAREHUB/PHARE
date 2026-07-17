// This test will segfault with a regular MPI not setup for GPU/ROCM/CUDA explictely
// check output of `ompi_info | grep Ext` for supported modules

#include "core/utilities/span.hpp"
#include "core/utilities/types.hpp"
#include "core/utilities/mpi_utils.hpp"

#include "gtest/gtest.h"

struct MpiCommsTest : public ::testing::Test
{
};

TEST(MpiCommsTest, collect_gpu_mem)
{
    // using T             = int;
    // using ManagedVector = std::vector<T, mkn::gpu::ManagedAllocator<T>>;

    std::size_t static constexpr SIZE = 10;
    auto mpi_size                     = PHARE::core::mpi::size();
    auto mpi_rank                     = PHARE::core::mpi::rank();
    std::size_t full_size             = SIZE * mpi_size;
    mkn::gpu::DeviceMem<int> devA{std::vector<int>(SIZE, mpi_rank + 1)}, datas{full_size};
    PHARE::core::mpi::collect(devA, datas);
    int expected = 0;
    for (std::int32_t i = 0; i < mpi_size; ++i)
        expected += (i + 1) * SIZE;
    EXPECT_EQ(expected, PHARE::core::sum(datas()));
}

TEST(MpiCommsTest, is_pure_gpu_mem_mpi_sendable)
{
    std::size_t static constexpr SIZE = 10;
    auto mpi_size                     = PHARE::core::mpi::size();
    auto mpi_rank                     = PHARE::core::mpi::rank();
    std::size_t full_size             = SIZE * mpi_size;
    mkn::gpu::DeviceMem<int> devA{std::vector<int>(SIZE, mpi_rank + 1)}, datas{full_size};
    PHARE::core::mpi::_collect(devA.data(), datas, devA.size(), devA.size());
    int expected = 0;
    for (std::int32_t i = 0; i < mpi_size; ++i)
        expected += (i + 1) * SIZE;
    EXPECT_EQ(expected, PHARE::core::sum(datas()));
}

TEST(MpiCommsTest, is_pure_gpu_mem_mpi_sendable_retrievable_on_cpu)
{
    std::size_t static constexpr SIZE = 10;
    auto mpi_size                     = PHARE::core::mpi::size();
    auto mpi_rank                     = PHARE::core::mpi::rank();
    std::size_t full_size             = SIZE * mpi_size;
    mkn::gpu::DeviceMem<int> devA{std::vector<int>(SIZE, mpi_rank + 1)};
    std::vector<int> datas(full_size);
    PHARE::core::mpi::_collect(devA.data(), datas, devA.size(), devA.size());
    int expected = 0;
    for (std::int32_t i = 0; i < mpi_size; ++i)
        expected += (i + 1) * SIZE;
    EXPECT_EQ(expected, PHARE::core::sum(datas));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    auto res = MPI_Init(&argc, &argv);
    if (res != MPI_SUCCESS)
        throw std::runtime_error("MPI init failed");
    res = RUN_ALL_TESTS();
    MPI_Finalize();
    return res;
}
