#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "data/field/refine/field_data_linear_refine.h"
#include "data/field/refine/field_linear_refine.h"
#include "data/grid/gridlayout.h"

#include "test_basic_hierarchy.h"


#include <functional>
#include <iostream>
#include <numeric>

using namespace PHARE;

using GridYee1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using Field1D     = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

using testing::Eq;



TEST(primalRefinementWithEvenRatio, givesRatioPlusOneRefinementPoints)
{
    auto centering   = QtyCentering::primal;
    auto ratio       = 4;
    auto isEvenRatio = true;
    EXPECT_EQ(ratio, nbrRefinePoints(centering, isEvenRatio, ratio));
}



TEST(dualRefinementWithEvenRatio, givesRatioRefinementPoints)
{
    auto centering   = QtyCentering::dual;
    auto ratio       = 8;
    auto isEvenRatio = true;
    EXPECT_EQ(ratio, nbrRefinePoints(centering, isEvenRatio, ratio));
}



TEST(primalRefinementWithOddRatio, givesRatioPlusOneRefinementPoints)
{
    auto centering   = QtyCentering::primal;
    auto ratio       = 5;
    auto isEvenRatio = false;
    EXPECT_EQ(ratio, nbrRefinePoints(centering, isEvenRatio, ratio));
}



TEST(dualRefinementWithOddRatio, givesRatioPlusOneRefinementPoints)
{
    auto centering   = QtyCentering::primal;
    auto ratio       = 5;
    auto isEvenRatio = false;
    EXPECT_EQ(ratio, nbrRefinePoints(centering, isEvenRatio, ratio));
}



TEST(UniformIntervalPartition, givesCorrectPartitionsForPrimalEven)
{
    auto nbrPoints = nbrRefinePoints(QtyCentering::primal, true, 2);
    UniformIntervalPartitionWeight uipw{QtyCentering::primal, 2, nbrPoints};

    std::array<double, 2> expectedDistances{0, 0.5};
    auto const& actualDistances = uipw.getUniformDistances();

    for (auto i = 0u; i < 2; ++i)
    {
        EXPECT_DOUBLE_EQ(expectedDistances[i], actualDistances[i]);
    }
}



TEST(UniformIntervalPartition, givesCorrectPartitionsForPrimalOdd)
{
    UniformIntervalPartitionWeight uipw{QtyCentering::primal, 5,
                                        nbrRefinePoints(QtyCentering::primal, false, 5)};

    std::array<double, 5> expectedDistances{0., 1. / 5., 2. / 5., 3. / 5., 4. / 5.};
    auto const& actualDistances = uipw.getUniformDistances();

    for (auto i = 0u; i < 5; ++i)
    {
        EXPECT_DOUBLE_EQ(expectedDistances[i], actualDistances[i]);
    }
}


TEST(UniformIntervalPartition, givesCorrectPartitionsForDualEven)
{
    auto const ratio = 4u;
    auto nbrPoints   = nbrRefinePoints(QtyCentering::dual, ratio % 2 == 0, ratio);

    UniformIntervalPartitionWeight uipw{QtyCentering::dual, ratio, nbrPoints};

    auto smallCellSize = 1. / ratio;
    std::array<double, 4> expectedDistances{5. / 2. * smallCellSize, 7. / 2. * smallCellSize,
                                            1. / 2. * smallCellSize, 3. / 2. * smallCellSize};

    auto const& actualDistances = uipw.getUniformDistances();

    for (auto i = 0u; i < nbrPoints; ++i)
    {
        EXPECT_DOUBLE_EQ(expectedDistances[i], actualDistances[i]);
    }
}



TEST(UniformIntervalPartition, givesCorrectPartitionsForDualOdd)
{
    auto const ratio = 5u;
    auto nbrPoints   = nbrRefinePoints(QtyCentering::dual, ratio % 2 == 0, ratio);

    UniformIntervalPartitionWeight uipw{QtyCentering::dual, ratio, nbrPoints};
    auto smallCellSize = 1. / ratio;

    std::array<double, 5> expectedDistances{smallCellSize * 3, smallCellSize * 4, 0,
                                            smallCellSize * 1, smallCellSize * 2};

    auto const& actualDistances = uipw.getUniformDistances();

    for (auto i = 0u; i < nbrPoints; ++i)
    {
        EXPECT_DOUBLE_EQ(expectedDistances[i], actualDistances[i]);
    }
}




TEST(FieldRefineOperator, CanBeCreated)
{
    FieldRefineOperator<GridYee1DO1, Field1D> linearRefine{};
}




TEST(FieldRefine, CanBeCreated)
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
    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    std::array<Point<int, dimension>, 2> fineIndexes{Point<int, dimension>{12},
                                                     Point<int, dimension>{18}};
    std::array<Point<int, dimension>, 2> expectedStartIndexes{Point<int, dimension>{2},
                                                              Point<int, dimension>{3}};

    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex          = fineIndexes[i];
        auto expectedStartIndex = expectedStartIndexes[i];
        do
        {
            auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);
            std::cout << "fineIndex = " << fineIndex[dirX] << " startIndex = " << startIndex[dirX]
                      << " expected = " << expectedStartIndex[dirX] << "\n";

            EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
            ++fineIndex[dirX];

        } while (fineIndex[dirX] < fineIndexes[i][dirX] + ratio(dirX));
    }
}




TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForPrimalOddRatio)
{
    std::size_t constexpr dimension{1};
    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 9};
    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    std::array<Point<int, dimension>, 2> fineIndexes{Point<int, dimension>{18},
                                                     Point<int, dimension>{27}};
    std::array<Point<int, dimension>, 2> expectedStartIndexes{Point<int, dimension>{2},
                                                              Point<int, dimension>{3}};

    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex          = fineIndexes[i];
        auto expectedStartIndex = expectedStartIndexes[i];
        do
        {
            auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);
            std::cout << "fineIndex = " << fineIndex[dirX] << " startIndex = " << startIndex[dirX]
                      << " expected = " << expectedStartIndex[dirX] << "\n";

            EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
            ++fineIndex[dirX];

        } while (fineIndex[dirX] < fineIndexes[i][dirX] + ratio(dirX));
    }
}




TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForDualEvenRatio)
{
    std::size_t constexpr dimension{1};
    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 4};
    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    std::array<Point<int, dimension>, 2> fineIndexes{Point<int, dimension>{12},
                                                     Point<int, dimension>{16}};
    std::array<Point<int, dimension>, 2> expectedStartIndexes{Point<int, dimension>{2},
                                                              Point<int, dimension>{3}};

    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex          = fineIndexes[i];
        auto expectedStartIndex = expectedStartIndexes[i];
        do
        {
            auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

            if (fineIndex[dirX] < fineIndexes[i][dirX] + ratio(dirX) / 2)
            {
                std::cout << "fineIndex = " << fineIndex[dirX]
                          << " startIndex = " << startIndex[dirX]
                          << " expected = " << expectedStartIndex[dirX] << "\n";
                EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
            }
            else
            {
                std::cout << "fineIndex = " << fineIndex[dirX]
                          << " startIndex = " << startIndex[dirX]
                          << " expected = " << expectedStartIndex[dirX] + 1 << "\n";
                EXPECT_EQ(expectedStartIndex[dirX] + 1, startIndex[dirX]);
            }

            ++fineIndex[dirX];

        } while (fineIndex[dirX] < fineIndexes[i][dirX] + ratio(dirX));
    }
}




TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectStartIndexForDualOddRatio)
{
    std::size_t constexpr dimension{1};
    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 7};
    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    std::array<Point<int, dimension>, 2> fineIndexes{Point<int, dimension>{21},
                                                     Point<int, dimension>{28}};
    std::array<Point<int, dimension>, 2> expectedStartIndexes{Point<int, dimension>{3},
                                                              Point<int, dimension>{4}};

    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex          = fineIndexes[i];
        auto expectedStartIndex = expectedStartIndexes[i];
        do
        {
            auto startIndex = indexesAndWeights.computeStartIndexes(fineIndex);

            auto midCell = static_cast<int>(fineIndexes[i][dirX] + ratio(dirX) / 2);
            if (fineIndex[dirX] < midCell)
            {
                std::cout << "fineIndex = " << fineIndex[dirX]
                          << " startIndex = " << startIndex[dirX]
                          << " expected = " << expectedStartIndex[dirX] - 1 << "\n";
                EXPECT_EQ(expectedStartIndex[dirX] - 1, startIndex[dirX]);
            }
            else
            {
                std::cout << "fineIndex = " << fineIndex[dirX]
                          << " startIndex = " << startIndex[dirX]
                          << " expected = " << expectedStartIndex[dirX] << "\n";
                EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
            }

            ++fineIndex[dirX];

        } while (fineIndex[dirX] < fineIndexes[i][dirX] + ratio(dirX));
    }

    /*std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 5};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


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
    */
}




// now the weights
TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectWeightsForDualOddRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 5};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    auto weights = indexesAndWeights.getWeights();

    auto const& xWeights = weights[dirX];

    EXPECT_DOUBLE_EQ(xWeights[2][0], 1.);

    EXPECT_DOUBLE_EQ(xWeights[0][1], 1. - xWeights[4][1]);
    EXPECT_DOUBLE_EQ(xWeights[1][1], 1. - xWeights[3][1]);

    // EXPECT_DOUBLE_EQ(xWeights[5][1], 1. - xWeights[2][1]);
}




TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectWeightsForPrimalOddRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 5};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    auto weights = indexesAndWeights.getWeights();

    auto const& xWeights = weights[dirX];

    EXPECT_DOUBLE_EQ(xWeights.front()[1], 0.);
    EXPECT_DOUBLE_EQ(xWeights.back()[1], 1.);
    EXPECT_DOUBLE_EQ(xWeights[1][1], 1. - xWeights[4][1]);
    EXPECT_DOUBLE_EQ(xWeights[2][1], 1. - xWeights[3][1]);
}




TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectWeightsForDualEvenRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 6};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    auto weights = indexesAndWeights.getWeights();

    auto const& xWeights = weights[dirX];
    EXPECT_DOUBLE_EQ(xWeights[0][1], 1. - xWeights[5][1]);
    EXPECT_DOUBLE_EQ(xWeights[1][1], 1. - xWeights[4][1]);
    EXPECT_DOUBLE_EQ(xWeights[2][1], 1. - xWeights[3][1]);
}




TEST(AFieldLinearRefineIndexesAndWeights1D, giveACorrectWeightsForPrimalEvenRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 6};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};


    auto weights = indexesAndWeights.getWeights();

    auto const& xWeights = weights[dirX];

    EXPECT_DOUBLE_EQ(xWeights[0][1], 0.);
    EXPECT_DOUBLE_EQ(xWeights[6][1], 1.);
    EXPECT_DOUBLE_EQ(xWeights[1][1], 1. - xWeights[5][1]);
    EXPECT_DOUBLE_EQ(xWeights[2][1], 1. - xWeights[4][1]);
    EXPECT_DOUBLE_EQ(xWeights[3][1], 1. - xWeights[3][1]);
}




TEST(AFieldLinearIndexesAndWeights1D, giveACorrectiWeightForPrimalEvenRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 6};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    Point<int, dimension> fineIndex{12};

    auto iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];

    int lastIndex = 6;

    EXPECT_THAT(iWeight, Eq(0));

    for (int i = 1; i < lastIndex; ++i)
    {
        ++fineIndex[dirX];
        iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
        EXPECT_EQ(i, iWeight);
    }


    ++fineIndex[dirX];
    iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
    EXPECT_THAT(iWeight, Eq(0));
}




TEST(AFieldLinearIndexesAndWeights1D, giveACorrectiWeightForPrimalOddRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::primal}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 5};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    Point<int, dimension> fineIndex{10};

    auto iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];

    int lastIndex = 5;

    EXPECT_THAT(iWeight, Eq(0));

    for (int i = 1; i < lastIndex; ++i)
    {
        ++fineIndex[dirX];
        iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
        EXPECT_EQ(i, iWeight);
    }


    ++fineIndex[dirX];
    iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
    EXPECT_THAT(iWeight, Eq(0));
}




TEST(AFieldLinearIndexesAndWeights1D, giveACorrectiWeightForDualEvenRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 6};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    Point<int, dimension> fineIndex{12};

    auto iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];

    int lastIndex = 6;

    EXPECT_THAT(iWeight, Eq(0));

    for (int i = 1; i < lastIndex; ++i)
    {
        ++fineIndex[dirX];
        iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
        EXPECT_EQ(i, iWeight);
    }


    ++fineIndex[dirX];
    iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
    EXPECT_THAT(iWeight, Eq(0));
}




TEST(AFieldLinearIndexesAndWeights1D, giveACorrectiWeightForDualOddRatio)
{
    std::size_t constexpr dimension{1};

    std::array<QtyCentering, dimension> centering{{QtyCentering::dual}};

    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dimension}, 7};

    FieldRefineIndexesAndWeights<dimension> indexesAndWeights{centering, ratio};

    Point<int, dimension> fineIndex{14};

    auto iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];

    int lastIndex = 7;

    EXPECT_THAT(iWeight, Eq(0));

    for (int i = 1; i < lastIndex; ++i)
    {
        ++fineIndex[dirX];
        iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
        EXPECT_EQ(i, iWeight);
    }


    ++fineIndex[dirX];
    iWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];
    EXPECT_THAT(iWeight, Eq(0));
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
