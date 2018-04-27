#include <SAMRAI/geom/CartesianPatchGeometry.h>
#include <SAMRAI/pdat/CellData.h>
#include <SAMRAI/pdat/CellDataFactory.h>
#include <SAMRAI/pdat/NodeDataFactory.h>
#include <SAMRAI/pdat/NodeGeometry.h>
#include <SAMRAI/tbox/MessageStream.h>
#include <SAMRAI/tbox/SAMRAIManager.h>
#include <SAMRAI/tbox/SAMRAI_MPI.h>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "data/field/field.h"
#include "data/field/field_data.h"
#include "data/field/field_variable.h"
#include "data/grid/gridlayout.h"
#include "data/ndarray/ndarray_vector.h"
#include "field_data_test_param.h"
#include "utilities/point/point.h"


using testing::Eq;
namespace PHARE
{
using Field1D = Field<NdArrayVector1D<>, HybridQuantity::Scalar>;

using FieldDataTest1DOrder1 = FieldDataTestParam<Layout::Yee, 1, 1, Field1D>;
using FieldDataTest1DOrder2 = FieldDataTestParam<Layout::Yee, 1, 2, Field1D>;
using FieldDataTest1DOrder3 = FieldDataTestParam<Layout::Yee, 1, 3, Field1D>;

using FieldDataTestList
    = ::testing::Types<FieldDataTest1DOrder1, FieldDataTest1DOrder2, FieldDataTest1DOrder3>;



TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx);


TYPED_TEST_P(AFieldData1DCenteredOnEx, PackStreamLikeACellData)
{
    auto& patch0 = this->patch1d.patch0;
    auto& patch1 = this->patch1d.patch1;


    auto& layout0 = this->param.field0Data->gridLayout;

    auto const& ghosts = this->ghosts;

    auto& param = this->param;



    std::shared_ptr<SAMRAI::hier::BoxGeometry> cell0Geom
        = std::make_shared<SAMRAI::pdat::CellGeometry>(patch0.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> cell1Geom
        = std::make_shared<SAMRAI::pdat::CellGeometry>(patch1.getBox(), ghosts);

    int lower = 6;
    int upper = 9;

    if (layout0.order() == 2)
    {
        upper = 15;
    }
    else if (layout0.order() == 3)
    {
        upper = 20;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};
    SAMRAI::hier::Box fillBox{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};

    bool overwriteInterior{true};

    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(this->dim)};

    auto fieldOverlap
        = std::dynamic_pointer_cast<FieldOverlap<1>>(param.field0Geom->calculateOverlap(
            *param.field1Geom, srcMask, fillBox, overwriteInterior, transformation));

    auto cellOverlap
        = std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(cell0Geom->calculateOverlap(
            *cell1Geom, srcMask, fillBox, overwriteInterior, transformation));


    this->cell0Data->fillAll(0.0);
    this->cell1Data->fillAll(1.0);

    auto& field0 = param.field0Data->field;


    // put field 1 data on a stream
    // and unpack it in field 0

    SAMRAI::tbox::MessageStream fieldStream;
    param.field1Data->packStream(fieldStream, *fieldOverlap);

    SAMRAI::tbox::MessageStream fieldReadStream{fieldStream.getCurrentSize(),
                                                SAMRAI::tbox::MessageStream::Read,
                                                fieldStream.getBufferStart()};

    param.field0Data->unpackStream(fieldReadStream, *fieldOverlap);


    SAMRAI::tbox::MessageStream cellStream;
    this->cell1Data->packStream(cellStream, *cellOverlap);

    SAMRAI::tbox::MessageStream cellReadStream{cellStream.getCurrentSize(),
                                               SAMRAI::tbox::MessageStream::Read,
                                               cellStream.getBufferStart()};

    this->cell0Data->unpackStream(cellReadStream, *cellOverlap);


    auto iStart = param.field0Data->gridLayout.ghostStartIndex(field0, Direction::X);
    auto iEnd   = param.field0Data->gridLayout.ghostEndIndex(field0, Direction::X);


    double const* cellDataStart = this->cell0Data->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(field0(ix), Eq(cellDataStart[ix]));
    }
}


TYPED_TEST_P(AFieldData1DCenteredOnEx, PackStreamWithPeriodicsLikeACellData)
{
    auto& patch0 = this->patch1d.patch0;
    auto& patch1 = this->patch1d.patch1;

    auto& layout0 = this->param.field0Data->gridLayout;


    auto const& ghosts = this->ghosts;


    std::shared_ptr<SAMRAI::hier::BoxGeometry> cell0Geom
        = std::make_shared<SAMRAI::pdat::CellGeometry>(patch0.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> cell1Geom
        = std::make_shared<SAMRAI::pdat::CellGeometry>(patch1.getBox(), ghosts);
    bool overwriteInterior{true};

    SAMRAI::hier::Transformation transformation{patch0.getBox().lower() - patch1.getBox().upper()};

    SAMRAI::hier::Box srcMask{this->cell1Data->getBox()};
    SAMRAI::hier::Box fillMask{this->cell0Data->getGhostBox()};

    auto fieldOverlap
        = std::dynamic_pointer_cast<FieldOverlap<1>>(this->param.field0Geom->calculateOverlap(
            *this->param.field1Geom, srcMask, fillMask, overwriteInterior, transformation));

    auto cellOverlap
        = std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(cell0Geom->calculateOverlap(
            *cell1Geom, srcMask, fillMask, overwriteInterior, transformation));


    this->cell0Data->fillAll(0.0);
    this->cell1Data->fillAll(1.0);

    auto& field0 = this->param.field0Data->field;


    // put field 1 data on a stream
    // and unpack it in field 0

    SAMRAI::tbox::MessageStream fieldStream;
    this->param.field1Data->packStream(fieldStream, *fieldOverlap);

    SAMRAI::tbox::MessageStream fieldReadStream{fieldStream.getCurrentSize(),
                                                SAMRAI::tbox::MessageStream::Read,
                                                fieldStream.getBufferStart()};

    this->param.field0Data->unpackStream(fieldReadStream, *fieldOverlap);


    SAMRAI::tbox::MessageStream cellStream;
    this->cell1Data->packStream(cellStream, *cellOverlap);

    SAMRAI::tbox::MessageStream cellReadStream{cellStream.getCurrentSize(),
                                               SAMRAI::tbox::MessageStream::Read,
                                               cellStream.getBufferStart()};

    this->cell0Data->unpackStream(cellReadStream, *cellOverlap);


    auto iStart = layout0.ghostStartIndex(field0, Direction::X);
    auto iEnd   = layout0.ghostEndIndex(field0, Direction::X);


    double const* cellDataStart = this->cell0Data->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(field0(ix), Eq(cellDataStart[ix]));
    }
}


REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx, PackStreamLikeACellData,
                           PackStreamWithPeriodicsLikeACellData);


INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEx,
                              FieldDataTestList);



TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy);


TYPED_TEST_P(AFieldData1DCenteredOnEy, PackStreamLikeANodeData)
{
    auto& patch0 = this->patch1d.patch0;
    auto& patch1 = this->patch1d.patch1;

    auto& param = this->param;

    auto& layout0 = param.field0Data->gridLayout;



    auto const& ghosts = this->ghosts;


    std::shared_ptr<SAMRAI::hier::BoxGeometry> node0Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(patch0.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> node1Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(patch1.getBox(), ghosts);

    int lower = 6;
    int upper = 9;

    if (layout0.order() == 2)
    {
        upper = 15;
    }
    else if (layout0.order() == 3)
    {
        upper = 20;
    }

    SAMRAI::hier::Box srcMask{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};
    SAMRAI::hier::Box fillBox{SAMRAI::hier::Index{this->dim, lower},
                              SAMRAI::hier::Index{this->dim, upper}, this->blockId};

    bool overwriteInterior{true};

    SAMRAI::hier::Transformation transformation{SAMRAI::hier::IntVector::getZero(this->dim)};

    auto fieldOverlap
        = std::dynamic_pointer_cast<FieldOverlap<1>>(param.field0Geom->calculateOverlap(
            *param.field1Geom, srcMask, fillBox, overwriteInterior, transformation));

    auto nodeOverlap
        = std::dynamic_pointer_cast<SAMRAI::pdat::NodeOverlap>(node0Geom->calculateOverlap(
            *node1Geom, srcMask, fillBox, overwriteInterior, transformation));


    this->node0Data->fillAll(0.0);
    this->node1Data->fillAll(1.0);

    auto& field0 = param.field0Data->field;


    // put field 1 data on a stream
    // and unpack it in field 0

    SAMRAI::tbox::MessageStream fieldStream;
    param.field1Data->packStream(fieldStream, *fieldOverlap);

    SAMRAI::tbox::MessageStream fieldReadStream{fieldStream.getCurrentSize(),
                                                SAMRAI::tbox::MessageStream::Read,
                                                fieldStream.getBufferStart()};

    param.field0Data->unpackStream(fieldReadStream, *fieldOverlap);


    SAMRAI::tbox::MessageStream nodeStream;
    this->node1Data->packStream(nodeStream, *nodeOverlap);

    SAMRAI::tbox::MessageStream nodeReadStream{nodeStream.getCurrentSize(),
                                               SAMRAI::tbox::MessageStream::Read,
                                               nodeStream.getBufferStart()};

    this->node0Data->unpackStream(nodeReadStream, *nodeOverlap);


    auto iStart = layout0.ghostStartIndex(field0, Direction::X);
    auto iEnd   = layout0.ghostEndIndex(field0, Direction::X);


    double const* nodeDataStart = this->node0Data->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(field0(ix), Eq(nodeDataStart[ix]));
    }
}

TYPED_TEST_P(AFieldData1DCenteredOnEy, PackStreamWithPeriodicsLikeANodeData)
{
    auto& patch0 = this->patch1d.patch0;
    auto& patch1 = this->patch1d.patch1;

    auto& param = this->param;

    auto& layout0 = param.field0Data->gridLayout;



    auto const& ghosts = this->ghosts;


    std::shared_ptr<SAMRAI::hier::BoxGeometry> node0Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(patch0.getBox(), ghosts);

    std::shared_ptr<SAMRAI::hier::BoxGeometry> node1Geom
        = std::make_shared<SAMRAI::pdat::NodeGeometry>(patch1.getBox(), ghosts);


    bool overwriteInterior{true};



    SAMRAI::hier::IntVector shift{param.patch0.getBox().lower() - param.patch1.getBox().upper()};
    SAMRAI::hier::Transformation transformation{shift};


    SAMRAI::hier::Box srcMask{this->node1Data->getBox()};
    SAMRAI::hier::Box fillMask{this->node0Data->getGhostBox()};


    auto fieldOverlap
        = std::dynamic_pointer_cast<FieldOverlap<1>>(param.field0Geom->calculateOverlap(
            *param.field1Geom, srcMask, fillMask, overwriteInterior, transformation));

    auto nodeOverlap
        = std::dynamic_pointer_cast<SAMRAI::pdat::NodeOverlap>(node0Geom->calculateOverlap(
            *node1Geom, srcMask, fillMask, overwriteInterior, transformation));


    this->node0Data->fillAll(0.0);
    this->node1Data->fillAll(1.0);

    auto& field0 = param.field0Data->field;


    // put field 1 data on a stream
    // and unpack it in field 0

    SAMRAI::tbox::MessageStream fieldStream;
    param.field1Data->packStream(fieldStream, *fieldOverlap);

    SAMRAI::tbox::MessageStream fieldReadStream{fieldStream.getCurrentSize(),
                                                SAMRAI::tbox::MessageStream::Read,
                                                fieldStream.getBufferStart()};

    param.field0Data->unpackStream(fieldReadStream, *fieldOverlap);


    SAMRAI::tbox::MessageStream nodeStream;
    this->node1Data->packStream(nodeStream, *nodeOverlap);

    SAMRAI::tbox::MessageStream nodeReadStream{nodeStream.getCurrentSize(),
                                               SAMRAI::tbox::MessageStream::Read,
                                               nodeStream.getBufferStart()};

    this->node0Data->unpackStream(nodeReadStream, *nodeOverlap);


    auto iStart = layout0.ghostStartIndex(field0, Direction::X);
    auto iEnd   = layout0.ghostEndIndex(field0, Direction::X);


    double const* nodeDataStart = this->node0Data->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(field0(ix), Eq(nodeDataStart[ix]));
    }
}
REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEy, PackStreamLikeANodeData,
                           PackStreamWithPeriodicsLikeANodeData);


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
