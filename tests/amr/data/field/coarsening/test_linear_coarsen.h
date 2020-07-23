#ifndef PHARE_TEST_LINEAR_COARSEN_H
#define PHARE_TEST_LINEAR_COARSEN_H


#include "amr/data/field/coarsening/field_coarsen_index_weight.h"
#include "core/data/grid/gridlayout.h"
#include "core/data/grid/gridlayout_impl.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"



using testing::DoubleEq;
using testing::DoubleNear;
using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;

using GridYee1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using Field1D     = Field<NdArrayVector<1>, HybridQuantity::Scalar>;


/** @brief This structure will contain data for one function to be coarse.
 *  We get the fine and coarse data an the expected after coarse field
 *  for two centering : primal and dual
 */
struct FieldCoarsenTestData
{
    static constexpr int dimension{1};

    std::array<double, dimension> meshSizeCoarse{0.2};
    std::array<double, dimension> meshSizeFine{0.1};

    std::array<int, 2> coarseIndexesX{{5, 45}};
    std::array<int, 2> fineIndexesX{18, 28};


    Point<double, dimension> coarseOrigin{meshSizeCoarse[dirX] * coarseIndexesX[0]};
    Point<double, dimension> fineOrigin{meshSizeFine[dirX] * fineIndexesX[0]};

    int ghostWidth{5};

    std::array<std::uint32_t, dimension> nbrCellCoarse{40};
    int exCoarseNbrCell{40 + ghostWidth * 2};
    int eyCoarseNbrCell{41 + ghostWidth * 2};

    std::array<std::uint32_t, dimension> nbrCellFine{20};
    int exFineNbrCell{20 + ghostWidth * 2};
    int eyFineNbrCell{21 + ghostWidth * 2};


    std::shared_ptr<GridYee1DO1> coarseLayout;
    std::shared_ptr<GridYee1DO1> fineLayout;

    std::shared_ptr<Field1D> exFineValue;
    std::shared_ptr<Field1D> eyFineValue;

    std::shared_ptr<Field1D> exCoarseValue;
    std::shared_ptr<Field1D> eyCoarseValue;

    std::shared_ptr<Field1D> expectedExCoarseValue;
    std::shared_ptr<Field1D> expectedEyCoarseValue;

    HybridQuantity::Scalar exQuantity{HybridQuantity::Scalar::Ex};
    HybridQuantity::Scalar eyQuantity{HybridQuantity::Scalar::Ey};

    void init()
    {
        coarseLayout = std::make_shared<GridYee1DO1>(meshSizeCoarse, nbrCellCoarse, coarseOrigin);


        exCoarseValue
            = std::make_shared<Field1D>("Ex", exQuantity, coarseLayout->allocSize(exQuantity));

        expectedExCoarseValue
            = std::make_shared<Field1D>("Ex", exQuantity, coarseLayout->allocSize(exQuantity));

        eyCoarseValue
            = std::make_shared<Field1D>("Ey", eyQuantity, coarseLayout->allocSize(eyQuantity));

        expectedEyCoarseValue
            = std::make_shared<Field1D>("Ey", eyQuantity, coarseLayout->allocSize(eyQuantity));




        fineLayout = std::make_shared<GridYee1DO1>(meshSizeFine, nbrCellFine, fineOrigin);

        exFineValue
            = std::make_shared<Field1D>("Ex", exQuantity, fineLayout->allocSize(exQuantity));

        eyFineValue
            = std::make_shared<Field1D>("Ey", eyQuantity, fineLayout->allocSize(eyQuantity));
    }
};




struct AFieldCoarsenOperator : public testing::TestWithParam<FieldCoarsenTestData>
{
    void SetUp() override { param = GetParam(); }

    FieldCoarsenTestData param;

    double absError = 1.e-8;
};




TEST_P(AFieldCoarsenOperator, doTheExpectedCoarseningForEx)
{
    auto& layout = param.coarseLayout;

    auto iStartX = layout->ghostStartIndex(param.exQuantity, Direction::X);
    auto iEndX   = layout->ghostEndIndex(param.exQuantity, Direction::X);

    auto const& exCoarseValue         = *param.exCoarseValue;
    auto const& expectedExCoarseValue = *param.expectedExCoarseValue;


    for (auto ix = iStartX; ix <= iEndX; ++ix)
    {
        EXPECT_THAT(exCoarseValue(ix), DoubleNear(expectedExCoarseValue(ix), absError));
    }
}




TEST_P(AFieldCoarsenOperator, doTheExpectedCoarseningForEy)
{
    auto& layout = param.coarseLayout;

    auto iStartX = layout->ghostStartIndex(param.eyQuantity, Direction::X);
    auto iEndX   = layout->ghostEndIndex(param.eyQuantity, Direction::X);

    auto const& eyCoarseValue         = *param.eyCoarseValue;
    auto const& expectedEyCoarseValue = *param.expectedEyCoarseValue;


    for (auto ix = iStartX; ix <= iEndX; ++ix)
    {
        EXPECT_THAT(eyCoarseValue(ix), DoubleNear(expectedEyCoarseValue(ix), absError));
    }
}




std::vector<FieldCoarsenTestData> createParam()
{
    std::vector<FieldCoarsenTestData> coarsenData;

    std::ifstream exFine{"dual_fine_original1d.txt"};
    std::ifstream eyFine{"primal_fine_original1d.txt"};

    std::ifstream exCoarse{"dual_coarse_original1d.txt"};
    std::ifstream eyCoarse{"primal_coarse_original1d.txt"};

    std::ifstream exCoarseAfterCoarsening{"dual_coarse_linear_coarsed_1d.txt"};
    std::ifstream eyCoarseAfterCoarsening{"primal_coarse_linear_coarsed_1d.txt"};


    auto canContinue = [&exFine, &eyFine, &exCoarse, &eyCoarse, &exCoarseAfterCoarsening,
                        &eyCoarseAfterCoarsening]() {
        //
        bool eofFlags = exFine.eof() || eyFine.eof() || exCoarse.eof() || eyCoarse.eof()
                        || exCoarseAfterCoarsening.eof() || eyCoarseAfterCoarsening.eof();

        bool badFlags = exFine.bad() || eyFine.bad() || exCoarse.bad() || eyCoarse.bad()
                        || exCoarseAfterCoarsening.bad() || eyCoarseAfterCoarsening.bad();

        return !eofFlags && !badFlags;
    };


    while (canContinue())
    {
        coarsenData.emplace_back();

        auto& param = coarsenData.back();

        param.init();

        auto& exFineValue = *param.exFineValue;
        auto& eyFineValue = *param.eyFineValue;

        auto& exCoarseValue = *param.exCoarseValue;
        auto& eyCoarseValue = *param.eyCoarseValue;

        auto& expectedExCoarseValue = *param.expectedExCoarseValue;
        auto& expectedEyCoarseValue = *param.expectedEyCoarseValue;


        for (std::uint32_t ix = 0; ix < static_cast<std::uint32_t>(param.exFineNbrCell); ++ix)
        {
            double value;
            exFine >> value;
            exFineValue(ix) = value;
        }


        // did we encounter error during reading (including reaching end of file)
        if (!canContinue())
        {
            coarsenData.pop_back();
            break;
        }

        for (std::uint32_t ix = 0; ix < static_cast<std::uint32_t>(param.eyFineNbrCell); ++ix)
        {
            double value;
            eyFine >> value;
            eyFineValue(ix) = value;
        }

        for (std::uint32_t ix = 0; ix < static_cast<std::uint32_t>(param.exCoarseNbrCell); ++ix)
        {
            double value;

            exCoarse >> value;
            exCoarseValue(ix) = value;

            exCoarseAfterCoarsening >> value;
            expectedExCoarseValue(ix) = value;
        }
        for (std::uint32_t ix = 0; ix < static_cast<std::uint32_t>(param.eyCoarseNbrCell); ++ix)
        {
            double value;

            eyCoarse >> value;
            eyCoarseValue(ix) = value;

            eyCoarseAfterCoarsening >> value;
            expectedEyCoarseValue(ix) = value;
        }



        // finally apply the coarsening

        constexpr int dimension{1};
        SAMRAI::tbox::Dimension dim{dimension};


        // Here we only need the lower index
        // so we don't have to set the upper
        // of this boxes
        SAMRAI::hier::Box fineBoxEx{dim};
        SAMRAI::hier::Box coarseBoxEx{dim};

        SAMRAI::hier::Box fineBoxEy{dim};
        SAMRAI::hier::Box coarseBoxEy{dim};


        auto centeringEx = param.coarseLayout->centering(param.exQuantity);

        auto centeringEy = param.coarseLayout->centering(param.eyQuantity);


        fineBoxEx.setLower(dirX,
                           param.fineIndexesX[0]
                               - static_cast<int>(param.fineLayout->nbrGhosts(centeringEx[dirX])));
        coarseBoxEx.setLower(
            dirX, param.coarseIndexesX[0]
                      - static_cast<int>(param.coarseLayout->nbrGhosts(centeringEx[dirX])));



        fineBoxEy.setLower(dirX,
                           param.fineIndexesX[0]
                               - static_cast<int>(param.fineLayout->nbrGhosts(centeringEy[dirX])));
        coarseBoxEy.setLower(
            dirX, param.coarseIndexesX[0]
                      - static_cast<int>(param.coarseLayout->nbrGhosts(centeringEy[dirX])));


        SAMRAI::hier::IntVector ratio{dim, 2};


        FieldCoarsener<dimension> coarseItEx{centeringEx, fineBoxEx, coarseBoxEx, ratio};

        FieldCoarsener<dimension> coarseItEy{centeringEy, fineBoxEy, coarseBoxEy, ratio};



        std::vector<int> coarseIndexesEx;
        std::vector<int> coarseIndexesEy;

        // see python value for coarseIndex
        for (int index = 10; index < 15; ++index)
        {
            coarseIndexesEx.push_back(index);
        }
        for (int index = 10; index < 16; ++index)
        {
            coarseIndexesEy.push_back(index);
        }


        for (auto coarseIndex : coarseIndexesEx)
        {
            coarseItEx(*param.exFineValue, *param.exCoarseValue,
                       Point<int, dimension>{coarseIndex});
        }
        for (auto coarseIndex : coarseIndexesEy)
        {
            coarseItEy(*param.eyFineValue, *param.eyCoarseValue,
                       Point<int, dimension>{coarseIndex});
        }
    }

    return coarsenData;
}



INSTANTIATE_TEST_SUITE_P(TestWithMultipleFunctionThat, AFieldCoarsenOperator,
                         testing::ValuesIn(createParam()));



#endif
