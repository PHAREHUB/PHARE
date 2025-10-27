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
    using GridND  = Grid<NdArrayVector<dim>, HybridQuantity::Scalar>;

public:
    void SetUp() override
    {
        // create a BasicHierarchy with a refinement factor equal 2
        basicHierarchy_ = std::make_shared<BasicHierarchy<GridYee, GridND>>(refine);
    }

    void TearDown() override { basicHierarchy_->TearDown(); }

protected:
    std::shared_ptr<BasicHierarchy<GridYee, GridND>> basicHierarchy_;
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

    using GridYee = typename TestFixture::GridYee;
    using GridND  = typename TestFixture::GridND;

    auto& basicHierarchy = this->basicHierarchy_;
    auto& hierarchy      = basicHierarchy->getHierarchy();

    // Value is initialized with a affine function : ax + by + cz + d
    // where a, b, c & d are given in TagStrategy::affineFill implementation
    auto& affineFill = TagStrategy<GridYee, GridND>::affineFill;

    auto level = hierarchy.getPatchLevel(1);

    for (auto& patch : *level)
    {
        for (auto const& variablesId : basicHierarchy->getVariables())
        {
            auto const& dataId = variablesId.second;
            auto fieldData     = std::dynamic_pointer_cast<FieldData<GridYee, GridND>>(
                patch->getPatchData(dataId));

            auto& layout = fieldData->gridLayout;
            auto& field  = fieldData->field;

            for (auto const amr_idx : layout.AMRGhostBoxFor(field))
            {
                auto const position = layout.fieldNodeCoordinates(field, amr_idx);
                auto const lcl_idx  = layout.AMRToLocal(amr_idx);
                EXPECT_DOUBLE_EQ(field(lcl_idx), affineFill(position, dataId));
            }
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
