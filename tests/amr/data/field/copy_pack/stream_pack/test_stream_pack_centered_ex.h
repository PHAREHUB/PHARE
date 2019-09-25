#ifndef PHARE_TEST_STREAM_PACK_CENTERED_EX_H
#define PHARE_TEST_STREAM_PACK_CENTERED_EX_H


#include "field_data_test_param.h"

#include <SAMRAI/tbox/MessageStream.h>

#include "gmock/gmock.h"
#include "gtest/gtest.h"


using testing::Eq;

using namespace PHARE::core;
using namespace PHARE::amr;




template<typename T>
struct Setup1DCenteredOnEx
{
    AFieldData1DCenteredOnEx<T>& testReference;

    std::shared_ptr<SAMRAI::hier::BoxGeometry> destinationCellGeometry;
    std::shared_ptr<SAMRAI::hier::BoxGeometry> sourceCellGeometry;

    SAMRAI::hier::Box srcMask;
    SAMRAI::hier::Box fillMask;

    SAMRAI::hier::Transformation transformation;


    Setup1DCenteredOnEx(AFieldData1DCenteredOnEx<T>& instance, SAMRAI::hier::Box const& sourceMask,
                        SAMRAI::hier::Box const& destinationMask,
                        SAMRAI::hier::Transformation const& transform)
        : testReference{instance}
        , destinationCellGeometry{std::make_shared<SAMRAI::pdat::CellGeometry>(
              testReference.param.destinationPatch.getBox(), testReference.ghosts)}
        , sourceCellGeometry{std::make_shared<SAMRAI::pdat::CellGeometry>(
              testReference.param.sourcePatch.getBox(), testReference.ghosts)}
        , srcMask{sourceMask}

        , fillMask{destinationMask}
        , transformation{transform}
    {
        std::array<bool, 2> overwritePossibilites{{true, false}};

        // For overwriteInterior in (true,false):
        // fill source data with a cos, destination with a sin
        // compute the overlap from srcMask,destinationMask, transformation
        // that will be used for the transfert
        // then transmit the source patch to the destination patch with both cellData and fieldData
        // finnaly compare cellData and fieldData

        for (auto overwriteInterior : overwritePossibilites)
        {
            auto fieldOverlap = std::dynamic_pointer_cast<FieldOverlap<1>>(
                testReference.param.destinationFieldGeometry->calculateOverlap(
                    *testReference.param.sourceFieldGeometry, srcMask, fillMask, overwriteInterior,
                    transformation));

            auto cellOverlap = std::dynamic_pointer_cast<SAMRAI::pdat::CellOverlap>(
                destinationCellGeometry->calculateOverlap(*sourceCellGeometry, srcMask, fillMask,
                                                          overwriteInterior, transformation));


            auto& destinationField = testReference.param.destinationFieldData->field;
            auto& sourceField      = testReference.param.sourceFieldData->field;



            // As usual we will fill the destination and source with the help of
            // two different function (here: cos , and sin)
            double* destinationCellStart = testReference.destinationCellData->getPointer();


            // Since our data match a CellData we can use our gridlayout to get correct
            // boundary
            auto iStart = testReference.param.destinationFieldData->gridLayout.ghostStartIndex(
                destinationField, Direction::X);
            auto iEnd = testReference.param.destinationFieldData->gridLayout.ghostEndIndex(
                destinationField, Direction::X);

            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                destinationCellStart[ix] = testReference.param.destinationFill(ix);
            }

            double* sourceCellStart = testReference.sourceCellData->getPointer();

            iStart = testReference.param.sourceFieldData->gridLayout.ghostStartIndex(sourceField,
                                                                                     Direction::X);
            iEnd   = testReference.param.sourceFieldData->gridLayout.ghostEndIndex(sourceField,
                                                                                 Direction::X);


            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                sourceCellStart[ix] = testReference.param.sourceFill(ix);
            }

            // We have set our data, now is time to packStream into a messageStream
            // for both FieldData and CellData, and read from it to fill another
            // data


            SAMRAI::tbox::MessageStream fieldStream;
            testReference.param.sourceFieldData->packStream(fieldStream, *fieldOverlap);

            SAMRAI::tbox::MessageStream fieldReadStream{fieldStream.getCurrentSize(),
                                                        SAMRAI::tbox::MessageStream::Read,
                                                        fieldStream.getBufferStart()};

            testReference.param.destinationFieldData->unpackStream(fieldReadStream, *fieldOverlap);


            SAMRAI::tbox::MessageStream cellStream;
            testReference.sourceCellData->packStream(cellStream, *cellOverlap);

            SAMRAI::tbox::MessageStream cellReadStream{cellStream.getCurrentSize(),
                                                       SAMRAI::tbox::MessageStream::Read,
                                                       cellStream.getBufferStart()};

            testReference.destinationCellData->unpackStream(cellReadStream, *cellOverlap);

            iStart = testReference.param.destinationFieldData->gridLayout.ghostStartIndex(
                destinationField, Direction::X);
            iEnd = testReference.param.destinationFieldData->gridLayout.ghostEndIndex(
                destinationField, Direction::X);


            // Data has been transfered, now is time to check that we have the same values as
            // the CellData
            double const* cellDataStart = testReference.destinationCellData->getPointer();
            for (auto ix = iStart; ix <= iEnd; ++ix)
            {
                EXPECT_THAT(destinationField(ix), Eq(cellDataStart[ix]));
            }
            testReference.param.resetValues();
        }
    }
};


#endif
