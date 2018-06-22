#include "test_basic_hierarchy.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "data/field/field.h"
#include "data/field/field_geometry.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutdefs.h"
#include "tools/amr_utils.h"

using GridYee1DO1 = GridLayoutImplYee<1, 1>;
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


    basicHierarchy.coarseIt();



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

            auto box1 = currentBox;
            auto box2 = currentBox;


            // Here the box are the one of the refine boxes tags
            box1.setLower(dirX, 4);
            box1.setUpper(dirX, 15);


            box2.setLower(dirX, 30);
            box2.setUpper(dirX, 50);

            bool const withGhost{true};

            auto box1Layout
                = FieldGeometry<GridYee1DO1, decltype(qty)>::layoutFromBox(box1, layout);
            auto box2Layout
                = FieldGeometry<GridYee1DO1, decltype(qty)>::layoutFromBox(box2, layout);


            // we have to consider them in "FieldGeometry space "

            box1 = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(box1, qty, box1Layout,
                                                                         !withGhost);

            box2 = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(box2, qty, box2Layout,
                                                                         !withGhost);




            auto fieldBoxWithoutGhost = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(
                fieldData->getBox(), qty, layout, !withGhost);

            auto fieldBox = FieldGeometry<GridYee1DO1, decltype(qty)>::toFieldBox(
                fieldData->getBox(), qty, layout, withGhost);


            auto const box1Restrict = fieldBox * box1;
            auto const box2Restrict = fieldBox * box2;

            // We coarse value on the interior of the patch

            bool const box1OverlapWithFine{box1.intersects(fieldBoxWithoutGhost)};
            bool const box2OverlapWithFine{box2.intersects(fieldBoxWithoutGhost)};

            bool const isOverlapWithFine = box1OverlapWithFine || box2OverlapWithFine;


            auto localFieldBox1 = AMRToLocal(box1Restrict, fieldBox);

            auto localFieldBox2 = AMRToLocal(box2Restrict, fieldBox);


            SAMRAI::tbox::Dimension dim{1};


            for (uint32 ix = iStartX; ix <= iEndX; ++ix)
            {
                auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);


                auto isInBox1 = [&localFieldBox1, ix, &dim]() {
                    return localFieldBox1.contains(SAMRAI::hier::Index{dim, static_cast<int>(ix)});
                };
                auto isInBox2 = [&localFieldBox2, ix, &dim]() {
                    return localFieldBox2.contains(SAMRAI::hier::Index{dim, static_cast<int>(ix)});
                };

                auto inBoxExpectedTest = [&field, ix, &affineFill, &position]() {
                    EXPECT_DOUBLE_EQ(field(ix), affineFill(position));
                };

                auto outBoxExpectedTest = [&field, ix]() { EXPECT_DOUBLE_EQ(field(ix), 0.); };



                // so we want to know if both box1 and box2 intersect this patch
                // in this case localFieldBox1 and localFieldBox2 have correct value
                // for indexing with the field
                if (box1OverlapWithFine && box2OverlapWithFine)
                {
                    if (isInBox1())
                    {
                        inBoxExpectedTest();
                    }
                    else if (isInBox2())
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
                    if (box1OverlapWithFine)
                    {
                        if (isInBox1())
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
                        if (isInBox2())
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




INSTANTIATE_TEST_CASE_P(WithRatioFrom2To10TestThat, ALinearFieldCoarsen1DO1,
                        ::testing::Values(2, 3, 4, 5, 6, 7, 8, 9, 10));
