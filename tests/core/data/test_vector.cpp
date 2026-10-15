#include "core/vector.hpp"
#include "core/utilities/types.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


auto constexpr static allocator_mode = PHARE::default_allocator_mode();

template<typename T>
using Vector  = PHARE::Vector<T, allocator_mode>;
using Vectors = ::testing::Types<Vector<int>, Vector<double>>;

template<typename Vector>
struct VectorTest : public ::testing::Test
{
};

TYPED_TEST_SUITE(VectorTest, Vectors, );


template<std::size_t dim>
using DimConst  = PHARE::core::DimConst<dim>;
using DimConsts = ::testing::Types<DimConst<1>, DimConst<2>>;

template<typename Vector>
struct DimConstTest : public ::testing::Test
{
};

TYPED_TEST_SUITE(DimConstTest, DimConsts, );

//
#if PHARE_HAVE_GPU && PHARE_HAVE_THRUST
#include "test_gpu_vector.hpp"
#endif

TYPED_TEST(VectorTest, is_constructable)
{
    auto vec = TypeParam::make(10);
    EXPECT_EQ(vec.size(), 10ull);
    EXPECT_EQ(PHARE::core::sum(vec), 0);
}

TYPED_TEST(VectorTest, is_fillable)
{
    auto vec = TypeParam::make(10);
    EXPECT_EQ(PHARE::core::sum(vec), 0);

    TypeParam::fill(vec, 1);
    EXPECT_EQ(PHARE::core::sum(vec), 10);
}

TYPED_TEST(VectorTest, is_copyable)
{
    auto vec0 = TypeParam::make(10);
    TypeParam::fill(vec0, 1);
    auto vec1 = vec0;
    EXPECT_EQ(vec0, vec1);
}

// TYPED_TEST(VectorTest, is_movable)
// {
//     auto vec0 = TypeParam::make(10);
//     TypeParam::fill(vec0, 1);
//     auto vec1 = std::move(vec0);

//     EXPECT_EQ(vec0.size(), 0);
//     EXPECT_EQ(vec1.size(), 10);
//     EXPECT_EQ(PHARE::core::sum(vec1), 10);
// }


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
