
#include "core/def/phare_mpi.hpp"
#include "phare_core.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/grid/gridlayout.hpp"
#include "core/data/grid/gridlayoutimplyee.hpp"

#include "amr/data/field/refine/field_linear_refine.hpp"
#include "amr/data/field/refine/field_refine_operator.hpp"
#include "amr/data/field/refine/field_refiner.hpp"
#include "amr/data/field/refine/composite_field_refiner.hpp"
#include "amr/data/field/refine/magnetic_composite_refiner.hpp"

#include <SAMRAI/tbox/SAMRAI_MPI.h>
#include <SAMRAI/tbox/SAMRAIManager.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include <array>
#include <cmath>
#include <limits>



using namespace PHARE::core;
using namespace PHARE::amr;

using testing::Eq;



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

    using GridYee = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    using GridT   = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

    FieldRefineOperator<GridYee, GridT, DefaultFieldRefiner<dim>> linearRefine{};
}


// instantiation gate: forces full compilation of CompositeFieldRefiner<...,order> (vtable →
// refineBox) and the additive KernelFieldRefineOperator across all dim/interp, order 2.
TYPED_TEST(aFieldRefineOperator, kernelRefineOperatorCanBeCreated)
{
    static constexpr auto dim    = typename TypeParam::first_type{}();
    static constexpr auto interp = typename TypeParam::second_type{}();

    using GridYee = PHARE::core::PHARE_Types<PHARE::SimOpts{dim, interp}>::Hybrid::GridLayout_t;
    using GridT   = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

    auto linearKernel = makeRefineKernel<GridYee, GridT>(2);
    EXPECT_NE(linearKernel, nullptr);

    KernelFieldRefineOperator<GridYee, GridT> kernelRefine{linearKernel};

    auto magLinearKernel = makeMagneticRefineKernel<GridYee, GridT>(2);
    EXPECT_NE(magLinearKernel, nullptr);

    EXPECT_ANY_THROW((makeRefineKernel<GridYee, GridT>(3)));
    EXPECT_ANY_THROW((makeRefineKernel<GridYee, GridT>(4)));
    EXPECT_ANY_THROW((makeMagneticRefineKernel<GridYee, GridT>(0)));
    EXPECT_ANY_THROW((makeMagneticRefineKernel<GridYee, GridT>(4)));
}




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




// ----- value-level refineBox tests for the composite (limited) kernel ----------------------------
//
// Boxes are placed with lower=0 so AMR == local indexing; ratio 2. The fine destination is filled
// with NaN so the kernel's NaN-guard writes every targeted index. The numeric core is separable, so
// 1-D exercises both primitives (dual ±¼ ladder, primal half-point); one 2-D magnetic
// case covers the sharedFacesOnly skip + the runtime tensor product of limited rows.

namespace
{
    using GridYee1D = PHARE::core::PHARE_Types<PHARE::SimOpts{1, 1}>::Hybrid::GridLayout_t;
    using Grid1D    = Grid<NdArrayVector<1>, HybridQuantity::Scalar>;
    using GridYee2D = PHARE::core::PHARE_Types<PHARE::SimOpts{2, 1}>::Hybrid::GridLayout_t;
    using Grid2D    = Grid<NdArrayVector<2>, HybridQuantity::Scalar>;

    template<std::size_t dim>
    SAMRAI::hier::Box boxOf(std::array<int, dim> lo, std::array<int, dim> up)
    {
        SAMRAI::tbox::Dimension d{dim};
        SAMRAI::hier::Index loi{d, 0}, upi{d, 0};
        for (std::size_t k = 0; k < dim; ++k)
        {
            loi(k) = lo[k];
            upi(k) = up[k];
        }
        return SAMRAI::hier::Box{loi, upi, SAMRAI::hier::BlockId{0}};
    }

    SAMRAI::hier::IntVector ratio2(std::size_t dim)
    {
        return SAMRAI::hier::IntVector{SAMRAI::tbox::Dimension{static_cast<unsigned short>(dim)}, 2};
    }

    constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
} // namespace


// the two fine children of a coarse dual cell always mean back to its average (conservation),
// at every order.
TEST(compositeRefiner1D, dualChildrenConserveCoarseAverage)
{
    std::array<double, 8> coarse = {1.0, 1.3, 0.4, 2.0, -0.5, 0.7, 0.9, 0.2};
    std::array<QtyCentering, 1> centering{QtyCentering::dual};

    auto run = [&](auto refiner) {
        Grid1D src{"c", HybridQuantity::Scalar::rho, 8u};
        Grid1D dst{"f", HybridQuantity::Scalar::rho, 16u};
        for (std::size_t i = 0; i < 8; ++i)
            src(i) = coarse[i];
        for (std::size_t i = 0; i < 16; ++i)
            dst(i) = NaN;

        refiner.refineBox(src, dst, boxOf<1>({4}, {11}), centering, boxOf<1>({0}, {15}),
                          boxOf<1>({0}, {7}), ratio2(1));

        for (int I = 2; I <= 5; ++I)
            EXPECT_NEAR(0.5 * (dst(2 * I) + dst(2 * I + 1)), coarse[I], 1e-12);
    };

    run(CompositeFieldRefiner<GridYee1D, Grid1D, 2>{});
}


// row sum ≡ 1 (partition of unity, S7): on a constant coarse field every fine child equals the
// constant.
TEST(compositeRefiner1D, dualRowSumsToOne)
{
    constexpr double konst = 3.7;
    std::array<QtyCentering, 1> centering{QtyCentering::dual};

    auto run = [&](auto refiner) {
        Grid1D src{"c", HybridQuantity::Scalar::rho, 8u};
        Grid1D dst{"f", HybridQuantity::Scalar::rho, 16u};
        for (std::size_t i = 0; i < 8; ++i)
            src(i) = konst;
        for (std::size_t i = 0; i < 16; ++i)
            dst(i) = NaN;
        refiner.refineBox(src, dst, boxOf<1>({4}, {11}), centering, boxOf<1>({0}, {15}),
                          boxOf<1>({0}, {7}), ratio2(1));
        for (int f = 4; f <= 11; ++f)
            EXPECT_NEAR(dst(f), konst, 1e-14);
    };

    run(CompositeFieldRefiner<GridYee1D, Grid1D, 2>{});
}


// B (Bx: primal-x normal, dual-y tangential), sharedFacesOnly: the two tangential (y) children of a
// shared (even-x) face are antisymmetric about the coarse value ⇒ sum = 2·C ⇒ ∇·B-neutral.
// Interior (odd-x) faces are left untouched (NaN) by the Tóth-Roe owner.
TEST(magneticCompositeRefiner2D, sharedFaceTangentialChildrenAreDivBNeutral)
{
    auto src_at = [](int ix, int iy) { return 1.0 * ix + 0.5 * iy + 0.1 * ix * iy; };
    std::array<QtyCentering, 2> centering{QtyCentering::primal, QtyCentering::dual}; // Bx

    auto run = [&](auto refiner) {
        Grid2D src{"c", HybridQuantity::Scalar::Bx, 6u, 6u};
        Grid2D dst{"f", HybridQuantity::Scalar::Bx, 12u, 12u};
        for (std::size_t ix = 0; ix < 6; ++ix)
            for (std::size_t iy = 0; iy < 6; ++iy)
                src(ix, iy) = src_at(ix, iy);
        for (std::size_t ix = 0; ix < 12; ++ix)
            for (std::size_t iy = 0; iy < 12; ++iy)
                dst(ix, iy) = NaN;

        refiner.refineBox(src, dst, boxOf<2>({4, 4}, {5, 7}), centering, boxOf<2>({0, 0}, {11, 11}),
                          boxOf<2>({0, 0}, {5, 5}), ratio2(2));

        // fine x=4 is a shared (even) face, anchor Ix=2; x=5 is interior (odd) → skipped (NaN)
        for (int Iy = 2; Iy <= 3; ++Iy)
        {
            double const c0 = dst(4, 2 * Iy);
            double const c1 = dst(4, 2 * Iy + 1);
            EXPECT_NEAR(c0 + c1, 2.0 * src_at(2, Iy), 1e-12);
        }
        EXPECT_TRUE(std::isnan(dst(5, 4))); // interior-normal face untouched
    };

    run(MagneticCompositeRefiner<GridYee2D, Grid2D, 2>{});
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
