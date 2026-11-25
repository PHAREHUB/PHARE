
#include "mpi/mpi_utils.hpp"

#include "gtest/gtest.h"

#include <array>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>

// every rank sends values made from its rank, so the checks cover both the content and the
// order of what is gathered; run with any number of ranks (ctest runs it with mpirun -n 2
// when configured with -DtestMPI=ON)

using namespace PHARE;

namespace
{
template<typename T>
T value_for(int const rank)
{
    if constexpr (std::is_same_v<T, char>)
        return static_cast<char>('a' + rank);
    else
        return static_cast<T>(rank + 1);
}

std::vector<double> vector_for(int const rank) // rank 0 sends nothing
{
    std::vector<double> v(rank);
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = rank * 10 + i;
    return v;
}
} // namespace


template<typename T>
struct MPICollectScalar : public ::testing::Test
{
};

using Scalars = ::testing::Types<double, float, int, std::uint32_t, std::int64_t, std::uint8_t,
                                 std::size_t, char>;
TYPED_TEST_SUITE(MPICollectScalar, Scalars);

TYPED_TEST(MPICollectScalar, gathers_one_value_per_rank_in_rank_order)
{
    auto const values = mpi::collect(value_for<TypeParam>(mpi::rank()));

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], value_for<TypeParam>(i));
}


TEST(MPICollect, strings_of_different_lengths_including_empty)
{
    auto const str_for = [](int const rank) { return std::string(rank, 'a' + rank); };

    auto const values = mpi::collect(str_for(mpi::rank()));

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], str_for(i));
}

TEST(MPICollect, vectors_of_different_sizes_including_empty)
{
    auto const values = mpi::collect(vector_for(mpi::rank()));

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], vector_for(i));
}

TEST(MPICollect, array_of_one)
{
    std::array<double, 1> const arr{mpi::rank() + .5};

    auto const values = mpi::collect(arr);

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], (std::array<double, 1>{i + .5}));
}

namespace
{
std::array<std::size_t, 3> array_for(int const rank)
{
    auto const r = static_cast<std::size_t>(rank);
    return {r * 3 + 1, r * 3 + 2, r * 3 + 1000}; // last value > 255: needs more than one byte
}
} // namespace

TEST(MPICollect, array_of_three)
{
    auto const values = mpi::collect(array_for(mpi::rank()));

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], array_for(i));
}

TEST(MPICollect, collectArrays_copies_every_element) // #1571
{
    auto const values = mpi::collectArrays(array_for(mpi::rank()), mpi::size());

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], array_for(i));
}

TEST(MPICollect, collectVector)
{
    auto const values = mpi::collectVector(vector_for(mpi::rank()));

    ASSERT_EQ(values.size(), static_cast<std::size_t>(mpi::size()));
    for (int i = 0; i < mpi::size(); ++i)
        EXPECT_EQ(values[i], vector_for(i));
}

TEST(MPICollect, collect_raw_and_collectSpanSet)
{
    auto const check = [](auto const& spans) {
        ASSERT_EQ(spans.sizes.size(), static_cast<std::size_t>(mpi::size()));
        for (int i = 0; i < mpi::size(); ++i)
        {
            auto const expected = vector_for(i);
            auto const span     = spans[i];
            ASSERT_EQ(static_cast<std::size_t>(span.size()), expected.size());
            EXPECT_TRUE(std::equal(span.begin(), span.end(), expected.begin()));
        }
    };

    auto const local = vector_for(mpi::rank());
    check(mpi::collect_raw(local, mpi::size()));
    check(mpi::collectSpanSet<double>(local));
}


TEST(MPIReduce, max_min_any)
{
    auto const last = mpi::size() - 1;

    EXPECT_EQ(mpi::max(static_cast<std::size_t>(mpi::rank())), static_cast<std::size_t>(last));
    EXPECT_EQ(mpi::min(mpi::rank() + 1.), 1.);
    EXPECT_TRUE(mpi::any(mpi::rank() == last));
    EXPECT_FALSE(mpi::any(false));
}


TEST(MPIAllGetFrom, scalar_string_and_vector_from_rank_0)
{
    EXPECT_EQ(mpi::all_get_from_rank_0([]() { return value_for<int>(mpi::rank()); }),
              value_for<int>(0));
    EXPECT_EQ(mpi::all_get_from_rank_0([]() { return std::string(5, 'a' + mpi::rank()); }),
              std::string(5, 'a'));
    EXPECT_EQ(mpi::all_get_from_rank_0([]() { return vector_for(mpi::rank() + 1); }),
              vector_for(1));
}

TEST(MPIAllGetFrom, scalar_string_and_vector_from_last_rank)
{
    auto const last = mpi::size() - 1;

    EXPECT_EQ(mpi::all_get_from(last, []() { return value_for<int>(mpi::rank()); }),
              value_for<int>(last));
    EXPECT_EQ(mpi::all_get_from(last, []() { return std::string(5, 'a' + mpi::rank()); }),
              std::string(5, 'a' + last));
    EXPECT_EQ(mpi::all_get_from(last, []() { return vector_for(mpi::rank() + 1); }),
              vector_for(last + 1));
}


TEST(MPIEnv, size_rank_barrier)
{
    EXPECT_GE(mpi::size(), 1);
    EXPECT_GE(mpi::rank(), 0);
    EXPECT_LT(mpi::rank(), mpi::size());
    mpi::barrier();
}

TEST(MPIEnv, is_root_rank_for_host_is_the_lowest_rank_per_host)
{
    auto const hosts = mpi::collect(mpi::hostname());
    auto const roots = mpi::collect(static_cast<int>(mpi::is_root_rank_for_host()));

    for (int i = 0; i < mpi::size(); ++i)
    {
        auto const first = std::find(hosts.begin(), hosts.end(), hosts[i]) - hosts.begin();
        EXPECT_EQ(roots[i], static_cast<int>(first == i));
    }
}

TEST(MPIEnv, timestamps_agree_across_ranks)
{
    auto const now   = mpi::collect(mpi::unix_timestamp_now());
    auto const dates = mpi::collect(mpi::date_time());

    for (int i = 1; i < mpi::size(); ++i)
    {
        EXPECT_EQ(now[i], now[0]);
        EXPECT_EQ(dates[i], dates[0]);
    }
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    auto const ret = RUN_ALL_TESTS();
    MPI_Finalize();
    return ret;
}
