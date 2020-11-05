#include "test_field_refinement_on_hierarchy.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "core/data/field/field.h"
#include "amr/data/field/field_geometry.h"
#include "core/data/grid/gridlayout_impl.h"
#include "core/data/grid/gridlayoutdefs.h"
#include "amr/resources_manager/amr_utils.h"



template<typename TypeInfo /*= std::pair<DimConst<1>, InterpConst<1>>*/>
struct ALinearFieldRefineTest : public ::testing::Test
{
    static constexpr auto dim    = typename TypeInfo::first_type{}();
    static constexpr auto interp = typename TypeInfo::second_type{}();
    static constexpr auto refine = 2;

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
    using FieldND = Field<NdArrayVector<dim>, HybridQuantity::Scalar>;

public:
    void SetUp() override
    {
        // create a BasicHierarchy with a refinement factor equal 2
        basicHierarchy_ = std::make_shared<BasicHierarchy<GridYee, FieldND>>(refine);
    }

    void TearDown() override { basicHierarchy_->TearDown(); }

protected:
    //
    std::shared_ptr<BasicHierarchy<GridYee, FieldND>> basicHierarchy_;
};


// https://stackoverflow.com/questions/56115790/gtest-parametrized-tests-for-different-types

using LinearFieldRefineTupleInfos = testing::Types<std::pair<DimConst<1>, InterpConst<1>>>;


TYPED_TEST_SUITE(ALinearFieldRefineTest, LinearFieldRefineTupleInfos);


TYPED_TEST(ALinearFieldRefineTest, ConserveLinearFunction)
{
    TypeParam pair;
    auto constexpr dim    = pair.first();
    auto constexpr interp = pair.second();

    using GridYee = GridLayout<GridLayoutImplYee<dim, interp>>;
    using FieldND = Field<NdArrayVector<dim>, HybridQuantity::Scalar>;


    auto& basicHierarchy = this->basicHierarchy_;
    auto& hierarchy      = basicHierarchy->getHierarchy();

    // Value is initialized with a affine function : ax + by + cz + d
    // where a, b, c & d are given in TagStrategy::affineFill implementation
    auto& affineFill = TagStrategy<GridYee, FieldND>::affineFill;

    // test coarse operation
    auto level = hierarchy.getPatchLevel(1);

    for (auto& patch : *level)
    {
        for (auto const& variablesId : basicHierarchy->getVariables())
        {
            auto const& dataId = variablesId.second;
            auto fieldData     = std::dynamic_pointer_cast<FieldData<GridYee, FieldND>>(
                patch->getPatchData(dataId));

            auto& layout = fieldData->gridLayout;
            auto& field  = fieldData->field;

            if constexpr (dim == 1)
            {
                std::uint32_t iStartX = layout.ghostStartIndex(field, Direction::X);
                std::uint32_t iEndX   = layout.ghostEndIndex(field, Direction::X);

                for (std::uint32_t ix = iStartX; ix <= iEndX; ++ix)
                {
                    auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);

                    EXPECT_DOUBLE_EQ(field(ix), affineFill(position));
                }
            }
            if constexpr (dim == 2) {}
            if constexpr (dim == 3) {}
        }
    }
}




// _________________________________________________________________________



/*
using GridYee1DO1 = GridLayout<GridLayoutImplYee<1, 1>>;
using Field1D     = Field<NdArrayVector<1>, HybridQuantity::Scalar>;


template<typename GridLayoutT, typename FieldT>
class ALinearFieldRefineZZZ : public testing::TestWithParam<int>
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




using ALinearFieldRefine1DO1 = ALinearFieldRefineZZZ<GridYee1DO1, Field1D>;




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

            std::uint32_t iStartX = layout.ghostStartIndex(field, Direction::X);
            std::uint32_t iEndX   = layout.ghostEndIndex(field, Direction::X);



            for (std::uint32_t ix = iStartX; ix <= iEndX; ++ix)
            {
                auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);

                EXPECT_DOUBLE_EQ(field(ix), affineFill(position));
            }
        }
    }
}




INSTANTIATE_TEST_SUITE_P(WithRatioFrom2To10TestThat, ALinearFieldRefine1DO1,
                         ::testing::Values(2, 3, 4, 5, 6, 7, 8, 9, 10));
*/
