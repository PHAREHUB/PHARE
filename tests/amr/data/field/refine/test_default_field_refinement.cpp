
#include "core/def/phare_mpi.hpp"

#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"

#include "amr/data/field/refine/field_linear_refine.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/field_refiner.hpp"

#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/tbox/SAMRAIManager.h>

#include "gtest/gtest.h"



using namespace PHARE::core;
using namespace PHARE::amr;


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


TEST(UniformIntervalPartition, givesCorrectPartitionsForPrimal)
{
    LinearWeighter linearWeighter{QtyCentering::primal, 2};
    std::array<double, 2> expectedDistances{0, 0.5};

    auto const& actualDistances = linearWeighter.getUniformDistances();

    for (auto i = 0u; i < 2; ++i)
    {
        EXPECT_DOUBLE_EQ(expectedDistances[i], actualDistances[i]);
    }
}


TEST(UniformIntervalPartition, givesCorrectPartitionsForDual)
{
    LinearWeighter linearWeighter{QtyCentering::dual, 2};
    std::array<double, 2> expectedDistances{0.75, 0.25};

    auto const& actualDistances = linearWeighter.getUniformDistances();

    for (auto i = 0u; i < 2; ++i)
    {
        EXPECT_DOUBLE_EQ(expectedDistances[i], actualDistances[i]);
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct aFieldRefineOperator : public ::testing::Test
{
};

using aFieldRefineOperatorInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<1>, InterpConst<2>>,
                     std::pair<DimConst<1>, InterpConst<3>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<2>, InterpConst<2>>, std::pair<DimConst<2>, InterpConst<3>>,
                     std::pair<DimConst<3>, InterpConst<1>>, std::pair<DimConst<3>, InterpConst<2>>,
                     std::pair<DimConst<3>, InterpConst<3>>>;

TYPED_TEST_SUITE(aFieldRefineOperator, aFieldRefineOperatorInfos);


TYPED_TEST(aFieldRefineOperator, canBeCreated)
{
    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
    using GridT   = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

    FieldRefineOperator<GridYee, GridT, DefaultFieldRefiner<dim>> linearRefine{};
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


template<typename dimType>
struct aFieldRefine : public testing::Test
{
};

using WithAllDim = testing::Types<DimConst<1>, DimConst<2>, DimConst<3>>;

TYPED_TEST_SUITE(aFieldRefine, WithAllDim);


TYPED_TEST(aFieldRefine, canBeCreated)
{
    static constexpr auto dim = TypeParam{}();

    SAMRAI::tbox::Dimension dimension{dim};
    std::array<QtyCentering, dim> centering = {{QtyCentering::primal}};
    SAMRAI::hier::Box destinationGhostBox{dimension};
    SAMRAI::hier::Box sourceGhostBox{dimension};
    SAMRAI::hier::IntVector ratio{dimension, 2};

    DefaultFieldRefiner<dim> fieldLinearRefine{centering, destinationGhostBox, sourceGhostBox,
                                               ratio};
}


// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


template<typename dimType>
struct aFieldLinearRefineIndexesAndWeights : public testing::Test
{
};

using WithAllDim = testing::Types<DimConst<1>, DimConst<2>, DimConst<3>>;

TYPED_TEST_SUITE(aFieldLinearRefineIndexesAndWeights, WithAllDim);



template<int dim, int numOfIndexes>
constexpr std::array<Point<int, dim>, numOfIndexes>
makeArrayOfPoints(std::array<int, numOfIndexes> indexVal)
{
    std::array<Point<int, dim>, numOfIndexes> arrayOfPoints{};

    for (auto i = 0u; i < numOfIndexes; ++i)
    {
        int index = indexVal[i];

        arrayOfPoints[i] = ConstArray<int, dim>(index);
    }
    return arrayOfPoints;
}



TYPED_TEST(aFieldLinearRefineIndexesAndWeights, giveACorrectStartIndexForPrimalQty)
{
    static constexpr auto dim = TypeParam{}();

    auto constexpr centering = ConstArray<QtyCentering, dim>(QtyCentering::primal);
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dim}, 2};
    FieldRefineIndexesAndWeights<dim> indexesAndWeights{centering, ratio};

    constexpr std::array<Point<int, dim>, 4> fineIndexes = makeArrayOfPoints<dim, 4>({-1, 0, 1, 2});
    constexpr std::array<Point<int, dim>, 4> expectedStartIndexes
        = makeArrayOfPoints<dim, 4>({-1, 0, 0, 1});


    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex          = fineIndexes[i];
        auto expectedStartIndex = expectedStartIndexes[i];

        if constexpr (dim == 1)
        {
            auto startIndex = indexesAndWeights.coarseStartIndex(fineIndex);

            EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
        }
        if constexpr (dim == 2)
        {
            auto startIndex = indexesAndWeights.coarseStartIndex(fineIndex);

            EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
            EXPECT_EQ(expectedStartIndex[dirY], startIndex[dirY]);
        }
        if constexpr (dim == 3)
        {
            auto startIndex = indexesAndWeights.coarseStartIndex(fineIndex);

            EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);
            EXPECT_EQ(expectedStartIndex[dirY], startIndex[dirY]);
            EXPECT_EQ(expectedStartIndex[dirZ], startIndex[dirZ]);
        }
    }
}


TYPED_TEST(aFieldLinearRefineIndexesAndWeights, giveACorrectStartIndexForDualQty)
{
    static constexpr auto dim = TypeParam{}();

    auto constexpr centering = ConstArray<QtyCentering, dim>(QtyCentering::dual);
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dim}, 2};
    FieldRefineIndexesAndWeights<dim> indexesAndWeights{centering, ratio};

    constexpr std::array<Point<int, dim>, 4> fineIndexes = makeArrayOfPoints<dim, 4>({-1, 0, 1, 2});
    constexpr std::array<Point<int, dim>, 4> expectedStartIndexes
        = makeArrayOfPoints<dim, 4>({-1, -1, 0, 0});


    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex          = fineIndexes[i];
        auto expectedStartIndex = expectedStartIndexes[i];

        auto startIndex = indexesAndWeights.coarseStartIndex(fineIndex);

        EXPECT_EQ(expectedStartIndex[dirX], startIndex[dirX]);

        if constexpr (dim > 1)
        {
            EXPECT_EQ(expectedStartIndex[dirY], startIndex[dirY]);
        }

        if constexpr (dim > 2)
        {
            EXPECT_EQ(expectedStartIndex[dirZ], startIndex[dirZ]);
        }
    }
}


TYPED_TEST(aFieldLinearRefineIndexesAndWeights, giveACorrectWeightsForPrimalQty)
{
    static constexpr auto dim = TypeParam{}();

    auto constexpr centering = ConstArray<QtyCentering, dim>(QtyCentering::primal);
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dim}, 2};
    FieldRefineIndexesAndWeights<dim> indexesAndWeights{centering, ratio};

    std::size_t constexpr primal = 0;
    std::size_t constexpr dual   = 1;


    auto xWeights = indexesAndWeights.weights(Direction::X);

    EXPECT_DOUBLE_EQ(xWeights[primal][1], 0.);
    EXPECT_DOUBLE_EQ(xWeights[primal][0], 1.);

    EXPECT_DOUBLE_EQ(xWeights[dual][1], 0.5);
    EXPECT_DOUBLE_EQ(xWeights[dual][0], 0.5);

    if constexpr (dim > 1)
    {
        auto yWeights = indexesAndWeights.weights(Direction::Y);

        EXPECT_DOUBLE_EQ(yWeights[primal][1], 0.);
        EXPECT_DOUBLE_EQ(yWeights[primal][0], 1.);

        EXPECT_DOUBLE_EQ(yWeights[dual][1], 0.5);
        EXPECT_DOUBLE_EQ(yWeights[dual][0], 0.5);
    }
    if constexpr (dim > 2)
    {
        auto zWeights = indexesAndWeights.weights(Direction::Z);

        EXPECT_DOUBLE_EQ(zWeights[primal][1], 0.);
        EXPECT_DOUBLE_EQ(zWeights[primal][0], 1.);

        EXPECT_DOUBLE_EQ(zWeights[dual][1], 0.5);
        EXPECT_DOUBLE_EQ(zWeights[dual][0], 0.5);
    }
}


TYPED_TEST(aFieldLinearRefineIndexesAndWeights, giveACorrectWeightsForDualQty)
{
    static constexpr auto dim = TypeParam{}();

    auto constexpr centering = ConstArray<QtyCentering, dim>(QtyCentering::dual);
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dim}, 2};
    FieldRefineIndexesAndWeights<dim> indexesAndWeights{centering, ratio};

    std::size_t constexpr primal = 0;
    std::size_t constexpr dual   = 1;


    auto xWeights = indexesAndWeights.weights(Direction::X);

    EXPECT_DOUBLE_EQ(xWeights[primal][1], 0.75);
    EXPECT_DOUBLE_EQ(xWeights[primal][0], 0.25);

    EXPECT_DOUBLE_EQ(xWeights[dual][1], 0.25);
    EXPECT_DOUBLE_EQ(xWeights[dual][0], 0.75);

    if constexpr (dim > 1)
    {
        auto yWeights = indexesAndWeights.weights(Direction::Y);

        EXPECT_DOUBLE_EQ(yWeights[primal][1], 0.75);
        EXPECT_DOUBLE_EQ(yWeights[primal][0], 0.25);

        EXPECT_DOUBLE_EQ(yWeights[dual][1], 0.25);
        EXPECT_DOUBLE_EQ(yWeights[dual][0], 0.75);
    }
    if constexpr (dim > 2)
    {
        auto zWeights = indexesAndWeights.weights(Direction::Z);

        EXPECT_DOUBLE_EQ(zWeights[primal][1], 0.75);
        EXPECT_DOUBLE_EQ(zWeights[primal][0], 0.25);

        EXPECT_DOUBLE_EQ(zWeights[dual][1], 0.25);
        EXPECT_DOUBLE_EQ(zWeights[dual][0], 0.75);
    }
}


TYPED_TEST(aFieldLinearRefineIndexesAndWeights, giveACorrectWeightIndexesForPrimalQty)
{
    static constexpr auto dim = TypeParam{}();

    auto constexpr centering = ConstArray<QtyCentering, dim>(QtyCentering::primal);
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dim}, 2};
    FieldRefineIndexesAndWeights<dim> indexesAndWeights{centering, ratio};

    constexpr std::array<Point<int, dim>, 4> fineIndexes = makeArrayOfPoints<dim, 4>({-1, 0, 1, 2});
    constexpr std::array<int, 4> expectedWeightIndexes{1, 0, 1, 0};


    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex           = fineIndexes[i];
        auto expectedWeightIndex = expectedWeightIndexes[i];

        auto xWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];

        EXPECT_EQ(expectedWeightIndex, xWeight);

        if constexpr (dim > 1)
        {
            auto yWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirY];

            EXPECT_EQ(expectedWeightIndex, yWeight);
        }

        if constexpr (dim > 2)
        {
            auto zWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirZ];

            EXPECT_EQ(expectedWeightIndex, zWeight);
        }
    }
}


TYPED_TEST(aFieldLinearRefineIndexesAndWeights, giveACorrectWeightIndexesForDualQty)
{
    static constexpr auto dim = TypeParam{}();

    auto constexpr centering = ConstArray<QtyCentering, dim>(QtyCentering::dual);
    SAMRAI::hier::IntVector ratio{SAMRAI::tbox::Dimension{dim}, 2};
    FieldRefineIndexesAndWeights<dim> indexesAndWeights{centering, ratio};

    constexpr std::array<Point<int, dim>, 4> fineIndexes = makeArrayOfPoints<dim, 4>({-1, 0, 1, 2});
    constexpr std::array<int, 4> expectedWeightIndexes{1, 0, 1, 0};


    for (auto i = 0u; i < fineIndexes.size(); ++i)
    {
        auto fineIndex           = fineIndexes[i];
        auto expectedWeightIndex = expectedWeightIndexes[i];

        auto xWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirX];

        EXPECT_EQ(expectedWeightIndex, xWeight);

        if constexpr (dim > 1)
        {
            auto yWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirY];

            EXPECT_EQ(expectedWeightIndex, yWeight);
        }

        if constexpr (dim > 2)
        {
            auto zWeight = indexesAndWeights.computeWeightIndex(fineIndex)[dirZ];

            EXPECT_EQ(expectedWeightIndex, zWeight);
        }
    }
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
