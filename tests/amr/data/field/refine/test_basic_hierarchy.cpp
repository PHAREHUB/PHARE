#include "test_basic_hierarchy.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "data/field/field.h"
#include "data/field/field_geometry.h"
#include "data/grid/gridlayout_impl.h"
#include "data/grid/gridlayoutdefs.h"
#include "tools/amr_utils.h"

using GridYee1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using Field1D     = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;


template<typename GridLayoutT, typename FieldT>
class ALinearFieldRefine : public testing::TestWithParam<int>
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




using ALinearFieldRefine1DO1 = ALinearFieldRefine<GridYee1DO1, Field1D>;




TEST_P(ALinearFieldRefine1DO1, conserveLinearFunction)
{
    auto& basicHierarchy = *basicHierarchy_;

    auto& hierarchy = basicHierarchy.getHierarchy();


    //  Value is initialized with a affine function : ax + b
    // where a = 0.5, b = 2.0
    // see TagStrategy::affineFill implementation

    auto& affineFill = TagStrategy<GridYee1DO1, Field1D>::affineFill;


    // test coarse operation
    auto level = hierarchy.getPatchLevel(1);

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

                EXPECT_DOUBLE_EQ(field(ix), affineFill(position));
            }
        }
    }
}




INSTANTIATE_TEST_CASE_P(WithRatioFrom2To10TestThat, ALinearFieldRefine1DO1,
                        ::testing::Values(2, 3, 4, 5, 6, 7, 8, 9, 10));
