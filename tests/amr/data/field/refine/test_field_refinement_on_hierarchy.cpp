/**
 * In this test, we check the linear property of the DEFAULT field refinement operator
 * Note this operator is used on all fields, including electric and magnetic fields
 * which, in production, use special operators.
 */

#include "phare_core.hpp"
#include "core/data/grid/grid.hpp"
#include "core/data/field/field.hpp"
#include "core/data/grid/gridlayout_impl.hpp"
#include "core/data/grid/gridlayoutdefs.hpp"
#include "core/data/ndarray/ndarray_vector.hpp"
#include "core/data/particles/particle_array_def.hpp"
#include "phare_simulator_options.hpp"

#include "amr/data/field/field_geometry.hpp"
#include "amr/data/field/field_geometry.hpp"
#include "amr/resources_manager/amr_utils.hpp"

#include "test_field_refinement_on_hierarchy.hpp"

#include "gtest/gtest.h"

using SimOpts = PHARE::SimOpts;

template<auto opts_>
struct TestParam
{
    auto static constexpr opts = opts_;
};


template<typename TestParam_t>
struct ALinearFieldRefineTest : public ::testing::Test
{
    auto static constexpr opts = TestParam_t::opts;

    static constexpr auto dim    = opts.dimension;
    static constexpr auto interp = opts.interp_order;
    static constexpr auto refine = 2;

    using PHARE_Types = PHARE::core::PHARE_Types<opts>;
    using GridYee     = PHARE_Types::GridLayout_t;
    using GridND      = PHARE_Types::Grid_t;

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


// clang-format off
using LinearFieldRefineTupleInfos = testing::Types<
    TestParam<SimOpts{1}>
   ,TestParam<SimOpts{2}>
   ,TestParam<SimOpts{3}>
PHARE_WITH_MKN_GPU(
   ,TestParam<SimOpts{1, 1, LayoutMode::AoSTS}>
   // ,TestParam<SimOpts{2, 1, LayoutMode::AoSTS}> // todo
   // ,TestParam<SimOpts{3, 1, LayoutMode::AoSTS}>
)

>;
// clang-format on

TYPED_TEST_SUITE(ALinearFieldRefineTest, LinearFieldRefineTupleInfos);


TYPED_TEST(ALinearFieldRefineTest, ConserveLinearFunction)
{
    auto constexpr dim    = TestFixture::dim;
    auto constexpr interp = TestFixture::interp;

    using GridYee = TestFixture::GridYee;
    using GridND  = TestFixture::GridND;

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

            if constexpr (dim == 1)
            {
                if constexpr (PHARE::core::is_field_tile_set_v<GridND>)
                {
                    for (auto const& tile : field())
                        for (auto const& bix : tile.layout().ghostBoxFor(field.physicalQuantity()))
                        {
                            auto position = tile.layout().fieldNodeCoordinates(
                                tile(), tile.layout().origin(), bix.as_signed());

                            EXPECT_DOUBLE_EQ(tile()(bix), affineFill(position, dataId))
                                << variablesId.first << " failed";
                        }
                }
                else
                {
                    for (auto const& bix : layout.ghostBoxFor(field.physicalQuantity()))
                    {
                        auto position
                            = layout.fieldNodeCoordinates(field, layout.origin(), bix.as_signed());
                        EXPECT_DOUBLE_EQ(field(bix), affineFill(position, dataId));
                    }
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
