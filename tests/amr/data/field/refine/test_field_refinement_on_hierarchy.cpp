#include "test_field_refinement_on_hierarchy.hpp"

#include "gmock/gmock.h"
#include "gtest/gtest.h"


#include "core/data/field/field.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "amr/resources_manager/amr_utils.hpp"



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
    std::shared_ptr<BasicHierarchy<GridYee, FieldND>> basicHierarchy_;
};


using LinearFieldRefineTupleInfos
    = testing::Types<std::pair<DimConst<1>, InterpConst<1>>, std::pair<DimConst<2>, InterpConst<1>>,
                     std::pair<DimConst<3>, InterpConst<1>>>;


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
                std::uint32_t gsi_X = layout.ghostStartIndex(field, Direction::X);
                std::uint32_t gei_X = layout.ghostEndIndex(field, Direction::X);

                for (std::uint32_t ix = gsi_X; ix <= gei_X; ++ix)
                {
                    auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix);

                    EXPECT_DOUBLE_EQ(field(ix), affineFill(position, dataId));
                }
            }
            if constexpr (dim == 2)
            {
                std::uint32_t gsi_X = layout.ghostStartIndex(field, Direction::X);
                std::uint32_t gei_X = layout.ghostEndIndex(field, Direction::X);
                std::uint32_t gsi_Y = layout.ghostStartIndex(field, Direction::Y);
                std::uint32_t gei_Y = layout.ghostEndIndex(field, Direction::Y);

                for (std::uint32_t ix = gsi_X; ix <= gei_X; ++ix)
                {
                    for (std::uint32_t iy = gsi_Y; iy <= gei_Y; ++iy)
                    {
                        auto position = layout.fieldNodeCoordinates(field, layout.origin(), ix, iy);

                        EXPECT_DOUBLE_EQ(field(ix, iy), affineFill(position, dataId));
                    }
                }
            }
            if constexpr (dim == 3)
            {
                std::uint32_t gsi_X = layout.ghostStartIndex(field, Direction::X);
                std::uint32_t gei_X = layout.ghostEndIndex(field, Direction::X);
                std::uint32_t gsi_Y = layout.ghostStartIndex(field, Direction::Y);
                std::uint32_t gei_Y = layout.ghostEndIndex(field, Direction::Y);
                std::uint32_t gsi_Z = layout.ghostStartIndex(field, Direction::Z);
                std::uint32_t gei_Z = layout.ghostEndIndex(field, Direction::Z);

                for (std::uint32_t ix = gsi_X; ix <= gei_X; ++ix)
                {
                    for (std::uint32_t iy = gsi_Y; iy <= gei_Y; ++iy)
                    {
                        for (std::uint32_t iz = gsi_Z; iz <= gei_Z; ++iz)
                        {
                            auto position
                                = layout.fieldNodeCoordinates(field, layout.origin(), ix, iy, iz);

                            EXPECT_DOUBLE_EQ(field(ix, iy, iz), affineFill(position, dataId));
                        }
                    }
                }
            }
        }
    }
}




INSTANTIATE_TEST_SUITE_P(WithRatioFrom2To10TestThat, ALinearFieldRefine1DO1,
                         ::testing::Values(2, 3, 4, 5, 6, 7, 8, 9, 10));



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
