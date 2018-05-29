#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/CellDataFactory.h>
#include <SAMRAI/pdat/NodeDataFactory.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/field/field_data.h"
#include "data/grid/gridlayout.h"
#include "data/grid/gridlayout_impl.h"
#include "field_data_test_param.h"
#include "utilities/point/point.h"


using testing::Eq;

namespace PHARE
{
using Field1D = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

using FieldDataTest1DOrder1 = FieldDataTestParam<GridLayoutImplYee<1, 1>, Field1D>;
using FieldDataTest1DOrder2 = FieldDataTestParam<GridLayoutImplYee<1, 2>, Field1D>;
using FieldDataTest1DOrder3 = FieldDataTestParam<GridLayoutImplYee<1, 3>, Field1D>;

using FieldDataTestList
    = ::testing::Types<FieldDataTest1DOrder1, FieldDataTest1DOrder2, FieldDataTest1DOrder3>;




TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx);

TYPED_TEST_P(AFieldData1DCenteredOnEx, CopyLikeACellData)
{
    this->cell0Data->fillAll(0.0);
    this->cell1Data->fillAll(1.0);

    auto& field0 = this->param.field0Data->field;


    this->param.field0Data->copy(*(this->param.field1Data));

    this->cell0Data->copy(*(this->cell1Data));

    auto iStart = this->param.field0Data->gridLayout.ghostStartIndex(field0, Direction::X);
    auto iEnd   = this->param.field0Data->gridLayout.ghostEndIndex(field0, Direction::X);


    double const* cellDataStart = this->cell0Data->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(field0(ix), Eq(cellDataStart[ix]));
    }
}

REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx, CopyLikeACellData);




INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEx,
                              FieldDataTestList);




TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy);

TYPED_TEST_P(AFieldData1DCenteredOnEy, CopyLikeACellData)
{
    this->node0Data->fillAll(0.0);
    this->node1Data->fillAll(1.0);

    auto& field0 = this->param.field0Data->field;


    this->param.field0Data->copy(*this->param.field1Data);

    this->node0Data->copy(*(this->node1Data));

    auto iStart = this->param.field0Data->gridLayout.ghostStartIndex(field0, Direction::X);
    auto iEnd   = this->param.field0Data->gridLayout.ghostEndIndex(field0, Direction::X);


    double const* nodeDataStart = this->node0Data->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(field0(ix), Eq(nodeDataStart[ix]));
    }
}

REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy, CopyLikeACellData);


INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEy,
                              FieldDataTestList);


} // namespace PHARE

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
