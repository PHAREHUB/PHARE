#include "field_data_test_param.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

using testing::Eq;

using namespace PHARE;

TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx);

TYPED_TEST_P(AFieldData1DCenteredOnEx, CopyLikeACellData)
{
    auto& destinationField = this->param.destinationFieldData->field;
    auto& sourceField      = this->param.sourceFieldData->field;

    double* destinationCellStart = this->destinationCellData->getPointer();


    auto iStart = this->param.destinationFieldData->gridLayout.ghostStartIndex(destinationField,
                                                                               Direction::X);
    auto iEnd   = this->param.destinationFieldData->gridLayout.ghostEndIndex(destinationField,
                                                                           Direction::X);

    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        destinationCellStart[ix] = this->param.destinationFill(ix);
    }

    double* sourceCellStart = this->sourceCellData->getPointer();

    iStart = this->param.sourceFieldData->gridLayout.ghostStartIndex(sourceField, Direction::X);
    iEnd   = this->param.sourceFieldData->gridLayout.ghostEndIndex(sourceField, Direction::X);


    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        sourceCellStart[ix] = this->param.sourceFill(ix);
    }


    this->param.destinationFieldData->copy(*(this->param.sourceFieldData));

    this->destinationCellData->copy(*(this->sourceCellData));

    iStart = this->param.destinationFieldData->gridLayout.ghostStartIndex(destinationField,
                                                                          Direction::X);
    iEnd   = this->param.destinationFieldData->gridLayout.ghostEndIndex(destinationField,
                                                                      Direction::X);


    double const* cellDataStart = this->destinationCellData->getPointer();
    for (auto ix = iStart; ix <= iEnd; ++ix)
    {
        EXPECT_THAT(destinationField(ix), Eq(cellDataStart[ix]));
    }
}

REGISTER_TYPED_TEST_CASE_P(AFieldData1DCenteredOnEx, CopyLikeACellData);




INSTANTIATE_TYPED_TEST_CASE_P(TestWithOrderFrom1To3That, AFieldData1DCenteredOnEx,
                              FieldDataTestList);
