#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <functional>
#include <numeric>

#include "data/field/field_data_linear_refine.h"
#include "data/refine/field_linear_refine.h"



using namespace PHARE;

using GridYee1DO1 = GridLayoutImplYee<1, 1>;
using Field1D     = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

TEST(FieldDataLinearRefine, CanBeCreated)
{
    FieldDataLinearRefine<GridYee1DO1, Field1D> linearRefine{};
}

TEST(FieldLinearRefine, CanBeCreated)
{
    constexpr std::size_t dimension{1};

    SAMRAI::tbox::Dimension dim{dimension};

    std::array<QtyCentering, dimension> centering = {{QtyCentering::primal}};

    SAMRAI::hier::Box destinationGhostBox{dim};

    SAMRAI::hier::Box sourceGhostBox{dim};

    SAMRAI::hier::IntVector ratio{dim, 2};


    FieldLinearRefine<1> fieldLinearRefine{centering, destinationGhostBox, sourceGhostBox, ratio};
}

TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForPrimalEvenRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 6};

    FieldLinearRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    Point<int, dimension> fineIndex{12};

    Point<int, dimension> expectedStartIndex{2};

    auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    do
    {
        ++fineIndex[dirX];

        startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

        EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    } while (fineIndex[dirX] < 17);

    ++fineIndex[dirX]; // 18

    ++expectedStartIndex[dirX]; // so now it is the next of 2

    startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
}

TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForPrimalOddRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 9};

    FieldLinearRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    Point<int, dimension> fineIndex{18};

    Point<int, dimension> expectedStartIndex{2};

    auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    do
    {
        ++fineIndex[dirX];

        startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

        EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    } while (fineIndex[dirX] < 26);

    ++fineIndex[dirX]; //

    ++expectedStartIndex[dirX]; // so now it is the next of 2

    startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
}

TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForDualEvenRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 4};

    FieldLinearRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    Point<int, dimension> fineIndex{8};

    Point<int, dimension> expectedStartIndex{1};

    auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    do
    {
        ++fineIndex[dirX];

        startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

        EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    } while (fineIndex[dirX] < 8);

    ++fineIndex[dirX]; //

    ++expectedStartIndex[dirX]; // so now it is the next of 2

    startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
}

TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForDualOddRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 5};

    FieldLinearRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    Point<int, dimension> fineIndex{10};

    Point<int, dimension> expectedStartIndex{1};

    auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    do
    {
        ++fineIndex[dirX];

        startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

        EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

    } while (fineIndex[dirX] < 12);

    ++fineIndex[dirX]; //

    ++expectedStartIndex[dirX]; // so now it is the next of 2

    startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

    EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
}


int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    SAMRAI::tbox::SAMRAI_MPI::init(&argc, &argv);
    SAMRAI::tbox::SAMRAIManager::initialize();
    SAMRAI::tbox::SAMRAIManager::startup();


    int testResult = RUN_ALL_TESTS();

    // Finalize
    SAMRAI::tbox::SAMRAIManager::shutdown();
    SAMRAI::tbox::SAMRAIManager::finalize();
    SAMRAI::tbox::SAMRAI_MPI::finalize();

    return testResult;
}
