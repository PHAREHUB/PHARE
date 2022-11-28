#include <string>
#include <vector>
#include <random>

#include "core/utilities/indexer.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using namespace PHARE::core;



TEST(Indexer, canAddAnIndex)
{
    Indexer indexer;
    indexer.add(0);
}



TEST(Indexer, sizeEqualsNumberOfRegisteredItems)
{
    Indexer indexer;
    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    indexer.add(3);

    EXPECT_EQ(4, indexer.size());

    indexer.add(4);
    indexer.add(5);
    indexer.add(6);
    indexer.add(7);

    EXPECT_EQ(8, indexer.size());
}

TEST(Indexer, canBeSorted)
{
    Indexer indexer;
    indexer.add(3);
    indexer.add(1);
    indexer.add(2);
    indexer.add(0);
    indexer.add(7);
    indexer.add(5);
    indexer.add(4);
    indexer.add(6);

    indexer.sort();
    auto idx                          = std::begin(indexer);
    std::array<int, 8> expectedSorted = {0, 1, 2, 3, 4, 5, 6, 7};
    for (std::size_t i = 0; i < expectedSorted.size(); ++i)
    {
        EXPECT_EQ(expectedSorted[i], *idx);
        ++idx;
    }
}


TEST(Indexer, canBeSortedLarge)
{
    Indexer indexer;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, 1000);
    std::unordered_map<int, int> counter;
    std::vector<int> indexes;
    for (std::size_t cpt = 0; cpt < 1000; ++cpt)
    {
        auto r = dis(gen);
        if (counter.count(r))
            counter[r]++;
        else
            counter[r] = 0;
    }
    for (auto const& [r, c] : counter)
    {
        if (counter[r] == 1)
        {
            indexes.push_back(r);
            indexer.add(r);
        }
    }
    indexer.sort();
    std::sort(std::begin(indexes), std::end(indexes));
    auto idx = std::begin(indexer);
    for (std::size_t i = 0; i < indexes.size(); ++i)
    {
        EXPECT_EQ(indexes[i], *idx);
        ++idx;
    }
}


TEST(Indexer, capacityIsLargerOrEqualThanSize)
{
    Indexer indexer;
    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    EXPECT_LE(indexer.size(), indexer.capacity());
}



TEST(Indexer, emptySetsSizeZeroLeavingCapacityUnchanged)
{
    Indexer indexer;
    EXPECT_EQ(0, indexer.capacity());
    EXPECT_EQ(0, indexer.size());

    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    indexer.add(3);

    EXPECT_EQ(4, indexer.size());
    EXPECT_FALSE(indexer.is_empty());
    auto capa = indexer.capacity();

    indexer.empty();
    EXPECT_TRUE(indexer.is_empty());

    EXPECT_EQ(0, indexer.size());
    EXPECT_EQ(capa, indexer.capacity());
}




TEST(Indexer, beginReturnsIteratorOnFirstElement)
{
    Indexer indexer;
    std::array<int, 4> values = {3, 4, 5, 6};
    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    indexer.add(3);
    int first = values[*std::begin(indexer)];
    EXPECT_EQ(3, first);
}

TEST(Indexer, iteratorPlusEqual)
{
    Indexer indexer;
    std::array<std::size_t, 4> indexes{4, 5, 6, 7};
    indexer.add(4);
    indexer.add(5);
    indexer.add(6);
    indexer.add(7);
    auto it = std::begin(indexer);
    EXPECT_EQ(4, *it);
    it = it + 1;
}



TEST(Indexer, iteratorPlusInt)
{
    Indexer indexer;
    std::array<std::size_t, 4> indexes{4, 5, 6, 7};
    indexer.add(4);
    indexer.add(5);
    indexer.add(6);
    indexer.add(7);
    auto it = std::begin(indexer);
    EXPECT_EQ(4, *it);
    it = it + 1;
    EXPECT_EQ(5, *it);
    it = std::begin(indexer);
    EXPECT_EQ(4, *it);
    for (std::size_t i = 0; i < 4; ++i)
    {
        EXPECT_EQ(*it, indexes[i]);
        it = it + 1;
    }

    it = std::begin(indexer);
    it = it + 3;
    EXPECT_EQ(7, *it);
    int a = -2;
    it    = it + a;
    EXPECT_EQ(5, *it);
}



TEST(Indexer, iteratorDecrement)
{
    Indexer indexer;
    std::array<std::size_t, 4> indexes{4, 5, 6, 7};
    indexer.add(4);
    indexer.add(5);
    indexer.add(6);
    indexer.add(7);
    auto it = std::end(indexer);
    --it;
    EXPECT_EQ(7, *it);
    --it;
    EXPECT_EQ(6, *it);
}


TEST(Indexer, iteratorMinus)
{
    Indexer indexer;
    std::array<std::size_t, 4> indexes{4, 5, 6, 7};
    indexer.add(4);
    indexer.add(5);
    indexer.add(6);
    indexer.add(7);
    auto it = std::end(indexer);
    it      = it - 1;
    EXPECT_EQ(7, *it);
    it = it - 1;
    EXPECT_EQ(6, *it);
    it = std::begin(indexer);
}

TEST(Indexer, iteratorLessThan)
{
    Indexer indexer;
    std::array<std::size_t, 4> indexes{4, 5, 6, 7};
    indexer.add(4);
    indexer.add(5);
    indexer.add(6);
    indexer.add(7);
    auto it  = std::begin(indexer);
    auto it2 = it + 1;
    EXPECT_LT(it, it2);
    auto it3 = it + 3;
    EXPECT_LT(it, it3);
}



TEST(Indexer, iteratorDifference)
{
    Indexer indexer;
    std::array<int, 4> values = {3, 4, 5, 6};
    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    indexer.add(3);
    indexer.add(4);
    auto it1 = std::begin(indexer);
    auto it2 = std::end(indexer);
    auto v   = it2 - it1;
    EXPECT_EQ(std::distance(std::begin(indexer), std::end(indexer)), indexer.size());
    EXPECT_EQ(it1 - it1 + 1, 1);
    auto it3 = it1 + 2;
    EXPECT_EQ(it3 - it1, 2);
}


TEST(Indexer, iteratorEquality)
{
    Indexer indexer;
    std::array<int, 4> values = {3, 4, 5, 6};
    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    indexer.add(3);
    indexer.add(4);
    auto it1 = std::begin(indexer);
    auto it2 = it1 + 2;
    EXPECT_EQ(it2, it1 + 2);
}


TEST(Indexer, endReturnsIteratorOnLastPlusOneElement)
{
    std::vector<int> values({1, 2, 3, 4, 5});
    Indexer indexer;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        indexer.add(itemIndex);

    int actual = 0;
    for (auto const& itemIndex : indexer)
    {
        actual = values[itemIndex];
    }
    EXPECT_EQ(actual, 5);
}

TEST(Indexer, singleElementIndexerHasCorrectEnd)
{
    std::vector<int> values({18});
    Indexer indexer;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        indexer.add(itemIndex);

    int actual = 0;
    for (auto const& itemIndex : indexer)
    {
        actual = values[itemIndex];
    }
    EXPECT_EQ(actual, 18);
}


TEST(Indexer, stdFindUsageOnIndexer)
{
    std::vector<int> values({18, 22, 43, 24});
    Indexer indexer;
    for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
        indexer.add(itemIndex);

    for (std::size_t idx = 0; idx < values.size(); ++idx)
    {
        auto it = std::find(std::begin(indexer), std::end(indexer), idx);
        auto i  = *it;
        EXPECT_EQ(values[i], values[idx]);
    }

    indexer.empty();
    indexer.add(2);
    auto it = std::find(std::begin(indexer), std::end(indexer), 2);
    EXPECT_EQ(*it, 2);
    auto it2 = std::find(std::begin(indexer), std::end(indexer), 3);
    EXPECT_FALSE(it2 != std::end(indexer));
}


TEST(Indexer, stdFindOnSingleItemBucket)
{
    Indexer indexer;
    indexer.add(2);
    auto it = std::find(std::begin(indexer), std::end(indexer), 2);
    EXPECT_EQ(*it, 2);
}

TEST(Indexer, loopOverEmptyBucketMakesNoIteration)
{
    std::vector<int> values({1, 2, 3, 4, 5});
    Indexer indexer;
    auto cpt = 0;
    for (auto const& itemIndex : indexer)
    {
        ++cpt;
    }
    EXPECT_EQ(0, cpt);
}




// TEST(Indexer, trimRemovesTheNLastEmptyBuckets)
//{
//     std::vector<int> values({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12});
//     Indexer bl;
//     for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
//         bl.add(itemIndex);
//
//     EXPECT_EQ(4 * 3, bl.capacity());
//     EXPECT_EQ(bl.size(), bl.capacity());
//
//     bl.empty();
//     std::vector<int> others({1, 2, 3, 4});
//     for (std::size_t itemIndex = 0; itemIndex < others.size(); ++itemIndex)
//         bl.add(itemIndex);
//
//     EXPECT_EQ(4, bl.size());
//     EXPECT_EQ(12, bl.capacity());
//
//     bl.trim(1);
//     EXPECT_EQ(9, bl.capacity());
//
//
//     bl.empty();
//     for (std::size_t itemIndex = 0; itemIndex < values.size(); ++itemIndex)
//         bl.add(itemIndex);
//     EXPECT_EQ(12, bl.capacity());
//     bl.trim(0);
//     EXPECT_EQ(12, bl.capacity());
// }



TEST(Indexer, removeAnElement)
{
    Indexer indexer;
    std::array<int, 4> values = {3, 4, 5, 6};

    indexer.add(0);
    indexer.add(1);
    indexer.add(2);
    indexer.add(3);
    indexer.remove(1);
    std::array<int, 3> expected{3, 5, 6};
    std::size_t i = 0;
    for (auto index : indexer)
    {
        EXPECT_EQ(values[index], expected[i++]);
    }
    EXPECT_EQ(indexer.size(), expected.size());
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
