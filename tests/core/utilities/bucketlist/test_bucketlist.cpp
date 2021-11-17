#include <string>
#include <vector>

#include "core/utilities/bucketlist.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;

struct Obj
{
};


TEST(BucketList, canBeGivenAParticle)
{
    BucketList<100, Obj> bl;
    bl.add(Obj{});
}



TEST(BucketList, registerMoreThanBucketSize)
{
    BucketList<2, Obj> bl;
    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});
}



TEST(BucketList, sizeEqualsNumberOfRegisteredItems)
{
    BucketList<3, Obj> bl;
    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});

    EXPECT_EQ(4, bl.size());

    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});

    EXPECT_EQ(8, bl.size());
}


TEST(BucketList, capcityEqualsTotalNbrOfMemorySlots)
{
    BucketList<3, Obj> bl;
    EXPECT_EQ(3, bl.capacity());
    BucketList<100, Obj> bl100;
    EXPECT_EQ(100, bl100.capacity());

    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});
    EXPECT_EQ(3, bl.capacity());
    EXPECT_EQ(bl.size(), bl.capacity());


    bl.add(Obj{});
    EXPECT_EQ(6, bl.capacity());
}



TEST(BucketList, emptySetSizeZeroLeavingCapacityUnchanged)
{
    BucketList<3, Obj> bl;
    EXPECT_EQ(3, bl.capacity());
    EXPECT_EQ(0, bl.size());

    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});
    bl.add(Obj{});

    EXPECT_EQ(4, bl.size());
    EXPECT_FALSE(bl.is_empty());

    bl.empty();
    EXPECT_TRUE(bl.is_empty());

    EXPECT_EQ(0, bl.size());
    EXPECT_EQ(6, bl.capacity());
}




TEST(BucketList, beginReturnsIteratorOnFirstElement)
{
    BucketList<3, int> bl;
    int a = 3, b = 4, c = 5, d = 6;
    bl.add(a);
    bl.add(b);
    bl.add(c);
    bl.add(d);
    int first = **(std::begin(bl));
    EXPECT_EQ(3, first);
}


TEST(BucketList, endReturnsIteratorOnLastPlusOneElement)
{
    std::vector<int> values({1, 2, 3, 4, 5});
    BucketList<3, int> bl;
    for (auto const& value : values)
        bl.add(value);

    int actual;
    for (auto const& v : bl)
    {
        actual = *v;
    }
    EXPECT_EQ(actual, 5);
}




TEST(BucketList, trimRemovesTheNLastEmptyBuckets)
{
    std::vector<int> values({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
    BucketList<3, int> bl;
    for (auto const& value : values)
        bl.add(value);

    EXPECT_EQ(4 * 3, bl.capacity());
    EXPECT_EQ(bl.size(), bl.capacity());

    bl.empty();
    std::vector<int> others({1, 2, 3, 4});
    for (auto const& value : others)
        bl.add(value);

    EXPECT_EQ(4, bl.size());
    EXPECT_EQ(12, bl.capacity());

    bl.trim(1);
    EXPECT_EQ(9, bl.capacity());


    bl.empty();
    for (auto const& value : values)
        bl.add(value);
    EXPECT_EQ(12, bl.capacity());
    bl.trim(0);
    EXPECT_EQ(12, bl.capacity());
}



int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}

