#include "test_basic_hierarchy.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "data/field/field.h"
#include "data/field/field_geometry.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutdefs.h"
#include "tools/amr_utils.h"

using GridYee1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using Field1D     = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;


template<typename GridLayoutT, typename FieldT>
class ALinearFieldCoarsen : public testing::TestWithParam<int>
{
public:
    void SetUp() override
    {
        // here we create a BasicHierarchy with a ratio GetParam()
        basicHierarchy_ = std::make_shared<BasicHierarchy<GridLayoutT, FieldT>>(this->GetParam());
    }

    void TearDown() override { basicHierarchy_->TearDown(); }

protected:
    //
    std::shared_ptr<BasicHierarchy<GridLayoutT, FieldT>> basicHierarchy_;
};

using ALinearFieldCoarsen1DO1 = ALinearFieldCoarsen<GridYee1DO1, Field1D>;




TEST_P(ALinearFieldCoarsen1DO1, conserveLinearFunction)
{
    auto& basicHierarchy = *basicHierarchy_;

    auto& hierarchy = basicHierarchy.getHierarchy();

    // init

    int const maxLevel = hierarchy.getNumberOfLevels();


    // fill root level(if withCoarseLevel) and upper level with a given function
    auto fill = [maxLevel, &hierarchy, &basicHierarchy](auto fillFunction,
                                                        bool withCoarseLevel = true) {
        int startLevel = 0;
        if (!withCoarseLevel)
        {
            startLevel = 1;
        }

        for (int iLevel = startLevel; iLevel < maxLevel; ++iLevel)
        {
            auto level = hierarchy.getPatchLevel(iLevel);

            for (auto& patch : *level)
            {
                for (auto const& variablesId : basicHierarchy.getVariables())
                {
                    auto const& dataId = variablesId.second;
                    auto fieldData     = std::dynamic_pointer_cast<FieldData<GridYee1DO1, Field1D>>(
                        patch->getPatchData(dataId));

                    auto& layout = fieldData->gridLayout;
                    auto& field  = fieldData->field;

                    // here we are 1D

                    uint32 iStartX = layout.ghostStartIndex(field, Direction::X);
                    uint32 iEndX   = layout.ghostEndIndex(field, Direction::X);

                    for (uint32 ix = iStartX; ix <= iEndX; ++ix)
                    {
                        auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);
                        field(ix)     = fillFunction(position);
                    }
                }
            }
        }
    };

    auto nullFunctionFill = [](Point<double, 1>) { return 0; };

    fill(nullFunctionFill);



    //  Initialize with a affine function : ax + b
    // where a = 0.5, b = 2.0

    auto affineFill = [](Point<double, 1> position) {
        double a = 0.5;
        double b = 2.0;
        return a * position[dirX] + b;
    };

    bool const withCoarseLevel = true;

    // we fill the fine level only with the affine function
    fill(affineFill, !withCoarseLevel);



    // apply the coarse operation


    basicHierarchy.coarsify();



    // test coarse operation
    auto level = hierarchy.getPatchLevel(0);

    for (auto& patch : *level)
    {
        for (auto const& variablesId : basicHierarchy.getVariables())
        {
            auto const& dataId = variablesId.second;
            auto fieldData     = std::dynamic_pointer_cast<FieldData<GridYee1DO1, Field1D>>(
                patch->getPatchData(dataId));

            auto& layout = fieldData->gridLayout;
            auto& field  = fieldData->field;

            auto qty = field.physicalQuantity();

            // here we are 1D

            uint32 iStartX = layout.physicalStartIndex(field, Direction::X);
            uint32 iEndX   = layout.physicalEndIndex(field, Direction::X);


            auto const& currentBox = patch->getBox();

            auto coarseBox1 = currentBox;
            auto coarseBox2 = currentBox;


            // Here the box are the one of the refine boxes tags
            coarseBox1.setLower(dirX, 4);
            coarseBox1.setUpper(dirX, 15);


            coarseBox2.setLower(dirX, 30);
            coarseBox2.setUpper(dirX, 50);

            bool const withGhost{true};

            auto box1Layout
                = FieldGeometry<GridYee1DO1, decltype(qty)>::layoutFromBox(coarseBox1, layout);
            auto box2Layout
                = FieldGeometry<GridYee1DO1, decltype(qty)>::layoutFromBox(coarseBox2, layout);


            // we have to consider them in "FieldGeometry space "

            coarseBox1 = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(
                coarseBox1, qty, box1Layout, !withGhost);

            coarseBox2 = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(
                coarseBox2, qty, box2Layout, !withGhost);




            auto fieldBoxWithoutGhost = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(
                fieldData->getBox(), qty, layout, !withGhost);

            auto fieldBox = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(
                fieldData->getBox(), qty, layout, withGhost);


            auto const coarseBox1Restrict = fieldBox * coarseBox1;
            auto const coarseBox2Restrict = fieldBox * coarseBox2;

            // We coarse value on the interior of the patch

            bool const coarseBox1OverlapWithFine{coarseBox1.intersects(fieldBoxWithoutGhost)};
            bool const coarseBox2OverlapWithFine{coarseBox2.intersects(fieldBoxWithoutGhost)};

            bool const isOverlapWithFine = coarseBox1OverlapWithFine || coarseBox2OverlapWithFine;


            auto localCoarseFieldBox1 = AMRToLocal(coarseBox1Restrict, fieldBox);

            auto localCoarseFieldBox2 = AMRToLocal(coarseBox2Restrict, fieldBox);


            SAMRAI::tbox::Dimension dim{1};


            for (uint32 ix = iStartX; ix <= iEndX; ++ix)
            {
                auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);


                auto isInCoarseBox1 = [&localCoarseFieldBox1, ix, &dim]() {
                    return localCoarseFieldBox1.contains(
                        SAMRAI::hier::Index{dim, static_cast<int>(ix)});
                };
                auto isInCoarseBox2 = [&localCoarseFieldBox2, ix, &dim]() {
                    return localCoarseFieldBox2.contains(
                        SAMRAI::hier::Index{dim, static_cast<int>(ix)});
                };

                auto inBoxExpectedTest = [&field, ix, &affineFill, &position]() {
                    EXPECT_DOUBLE_EQ(field(ix), affineFill(position));
                };

                auto outBoxExpectedTest = [&field, ix]() { EXPECT_DOUBLE_EQ(field(ix), 0.); };



                // so we want to know if both coarseBox1 and coarseBox2 intersect this patch
                // in this case localCoarseFieldBox1 and localCoarseFieldBox2 have correct value
                // for indexing with the field
                if (coarseBox1OverlapWithFine && coarseBox2OverlapWithFine)
                {
                    if (isInCoarseBox1())
                    {
                        inBoxExpectedTest();
                    }
                    else if (isInCoarseBox2())
                    {
                        inBoxExpectedTest();
                    }
                    else
                    {
                        outBoxExpectedTest();
                    }
                }
                // if only one box intersect this patch
                // we have to check the index with this particular box
                else if (isOverlapWithFine)
                {
                    if (coarseBox1OverlapWithFine)
                    {
                        if (isInCoarseBox1())
                        {
                            inBoxExpectedTest();
                        }
                        else
                        {
                            outBoxExpectedTest();
                        }
                    }
                    else
                    {
                        if (isInCoarseBox2())
                        {
                            inBoxExpectedTest();
                        }
                        else
                        {
                            outBoxExpectedTest();
                        }
                    }
                }
                // finnaly none of the box intersect this patch
                // so it must be 0.
                else
                {
                    outBoxExpectedTest();
                }
            }
        }
    }
}




INSTANTIATE_TEST_SUITE_P(WithRatioFrom2To10TestThat, ALinearFieldCoarsen1DO1,
                         ::testing::Values(2, 3, 4, 5, 6, 7, 8, 9, 10));
