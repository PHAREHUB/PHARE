#ifndef PHARE_TEST_STREAM_PACK_CENTERED_EY_H
#define PHARE_TEST_STREAM_PACK_CENTERED_EY_H


#include "field_data_test_param.h"

#include <SAMRAI/tbox/MessageStream.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"

using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr_interface;



template<typename T>
struct Setup1DCenteredOnEy
{
    AFieldData1DCenteredOnEy<T>& testReference;

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destinationNodeGeometry;
    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceNodeGeometry;

    SAMRAI::hier::Box srcMask;
    SAMRAI::hier::Box fillMask;

    SAMRAI::hier::Transformation transformation;


    Setup1DCenteredOnEy(AFieldData1DCenteredOnEy<T>& instance, SAMRAI::hier::Box const& sourceMask,
                        SAMRAI::hier::Box const& destinationMask,
                        SAMRAI::hier::Transformation const& transform)
        : testReference{instance}
        , destinationNodeGeometry{std::make_shared<SAMRAI::pdat::NodeGeometry>(
              testReference.param.destinationPatch.getBox(), testReference.ghosts)}
        , sourceNodeGeometry{std::make_shared<SAMRAI::pdat::NodeGeometry>(
              testReference.param.sourcePatch.getBox(), testReference.ghosts)}
        , srcMask{sourceMask}

        , fillMask{destinationMask}
        , transformation{transform}
    {
        std::array<bool, 2> overwritePossibilites{{true, false}};


        // For overwriteInterior in (true,false):
        // fill source data with cos , destination with sin
        // compute the overlap from srcMask,destinationMask, transformation
        // that will be used for the transfert
        // then transmit the source patch to the destination patch with both nodeData and fieldData
        // finnaly compare nodeData and fieldData

        for (auto overwriteInterior : overwritePossibilites)
        {
            auto fieldOverlap = std::dynamic_pointer_cast<FieldOverlap<1>>(
                testReference.param.destinationFieldGeometry->calculateOverlap(
                    *testReference.param.sourceFieldGeometry, srcMask, fillMask, overwriteInterior,
                    transformation));

            auto nodeOverlap = std::dynamic_pointer_cast<SAMRAI::pdat::NodeOverlap>(
                destinationNodeGeometry->calculateOverlap(*sourceNodeGeometry, srcMask, fillMask,
                                                          overwriteInterior, transformation));


            auto& destinationField = testReference.param.destinationFieldData->field;
            auto& sourceField      = testReference.param.sourceFieldData->field;


            // As usual we will fill the destination and source with the help of
            // two different function (here: cos , and sin)
            double* destinationNodeStart = testReference.destinationNodeData->getPointer();


            // Since our data match a NodeData we can use our gridlayout to get correct
            // boundary
            auto iStart = testReference.param.destinationFieldData->gridLayout.ghostStartIndex(
                destinationField, Direction::X);
            auto iEnd = testReference.param.destinationFieldData->gridLayout.ghostEndIndex(
                destinationField, Direction::X);

            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                destinationNodeStart[ix] = testReference.param.destinationFill(ix);
            }

            double* sourceNodeStart = testReference.sourceNodeData->getPointer();

            iStart = testReference.param.sourceFieldData->gridLayout.ghostStartIndex(sourceField,
                                                                                     Direction::X);
            iEnd   = testReference.param.sourceFieldData->gridLayout.ghostEndIndex(sourceField,
                                                                                 Direction::X);


            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                sourceNodeStart[ix] = testReference.param.sourceFill(ix);
            }

            // We have set our data, now is time to packStream into a messageStream
            // for both FieldData and NodeData, and read from it to fill another
            // data

            SAMRAI::tbox::MessageStream fieldStream;
            testReference.param.sourceFieldData->packStream(fieldStream, *fieldOverlap);

            SAMRAI::tbox::MessageStream fieldReadStream{fieldStream.getCurrentSize(),
                                                        SAMRAI::tbox::MessageStream::Read,
                                                        fieldStream.getBufferStart()};

            testReference.param.destinationFieldData->unpackStream(fieldReadStream, *fieldOverlap);


            SAMRAI::tbox::MessageStream nodeStream;
            testReference.sourceNodeData->packStream(nodeStream, *nodeOverlap);

            SAMRAI::tbox::MessageStream nodeReadStream{nodeStream.getCurrentSize(),
                                                       SAMRAI::tbox::MessageStream::Read,
                                                       nodeStream.getBufferStart()};

            testReference.destinationNodeData->unpackStream(nodeReadStream, *nodeOverlap);


            iStart = testReference.param.destinationFieldData->gridLayout.ghostStartIndex(
                destinationField, Direction::X);
            iEnd = testReference.param.destinationFieldData->gridLayout.ghostEndIndex(
                destinationField, Direction::X);


            // Data has been transfered, now is time to check that we have the same values as
            // the NodeData
            double const* nodeDataStart = testReference.destinationNodeData->getPointer();
            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                EXPECT_THAT(destinationField(ix), Eq(nodeDataStart[ix]));
            }
            testReference.param.resetValues();
        }
    }
};


#endif
